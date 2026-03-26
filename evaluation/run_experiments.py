import sys
import yaml
import subprocess
import itertools
from pathlib import Path
from tqdm import tqdm
import argparse

def load_yaml(path):
    if not Path(path).exists():
        print(f"❌ Config not found: {path}")
        sys.exit(1)
    with open(path, 'r') as f:
        return yaml.safe_load(f)

def get_files_from_dataset(conf_global, dataset_conf):
    """
    Returns a list of dicts: {'path': Path, 'subdir': str}
    'subdir' corresponds to the list name (e.g., "testset_no_noise") or "" if using a pattern.
    """
    base_dir = Path(conf_global['dataset_dir']) / dataset_conf['dir']
    ext = dataset_conf.get('extension', '.xyz')
    
    raw_pattern = dataset_conf.get('pattern', '*')
    if isinstance(raw_pattern, str):
        patterns = [raw_pattern]
    else:
        patterns = raw_pattern 

    files_data = []

    if not base_dir.exists():
        print(f"⚠️  Dataset directory missing: {base_dir}")
        return []

    # CASE 1: Using LISTS (Creates subdirectories)
    if 'lists' in dataset_conf:
        for list_name in dataset_conf['lists']:
            # Search for the list file
            list_path = base_dir / list_name
            if not list_path.exists(): list_path = base_dir / f"{list_name}.txt"
            if not list_path.exists(): list_path = Path(conf_global['dataset_dir']) / list_name
            if not list_path.exists(): list_path = Path(conf_global['dataset_dir']) / f"{list_name}.txt"

            if list_path.exists():
                subdir_name = Path(list_name).stem 

                with open(list_path, 'r') as f:
                    lines = [l.strip() for l in f.readlines() if l.strip()]
                    for l in lines:
                        fname = l if l.endswith(ext) else l + ext
                        fpath = base_dir / fname
                        if fpath.exists(): 
                            files_data.append({'path': fpath, 'subdir': subdir_name})
            else:
                print(f"⚠️  List not found: {list_name}")

    # CASE 2: Using PATTERNS (No subdirectory)
    else:
        found_files = []
        for pat in patterns:
            full_glob = f"{pat}{ext}"
            found = list(base_dir.glob(full_glob))
            if not found:
                print(f"⚠️  No files found for '{full_glob}' in {base_dir}")
            found_files.extend(found)
        
        found_files = list(set(found_files))
        
        for f in found_files:
            files_data.append({'path': f, 'subdir': ''})
    
    return files_data

def get_method_name(job, params):
    if 'method_label' in job: return job['method_label']
    if '-e' in params: return params['-e']
    if '--estimator' in params: return params['--estimator']
    return "default"

def get_clean_shape_name(file_stem):
    """
    Transforms 'shape_nbpoints' to 'shape' if the suffix is numeric.
    """
    if "_" in file_stem:
        parts = file_stem.rsplit("_", 1)
        if parts[1].isdigit():
            return parts[0]
    return file_stem

def generate_commands(job, file_info, conf_global, verbose=False):
    input_file = file_info['path']
    subdir = file_info['subdir']

    build_dir = Path(conf_global['build_dir'])
    exec_path = build_dir / job['executable']
    
    if not exec_path.exists() and exec_path.with_suffix(".exe").exists():
        exec_path = exec_path.with_suffix(".exe")
    if not exec_path.exists():
        raise FileNotFoundError(f"Executable not found: {exec_path}")

    save_config = job.get('save', conf_global.get('save', ["stats"]))

    full_grid = {}
    if 'grid' in job: full_grid.update(job['grid'])
    if 'noise_grid' in job: full_grid.update(job['noise_grid'])

    keys = list(full_grid.keys())
    values = list(full_grid.values())
    combinations = list(itertools.product(*values))

    tasks = []
    
    file_unique_name = input_file.stem 
    shape_group_name = get_clean_shape_name(file_unique_name)
    
    root_out = Path(conf_global['output_dir'])

    # --- HIERARCHY ---
    base_job_dir = root_out / job['name'] 
    
    # Stats Directory
    dir_stats = base_job_dir / "stats" / subdir
    
    if "stats" in save_config:
        dir_stats.mkdir(parents=True, exist_ok=True)
        stats_file = dir_stats / f"allErrorStats-{shape_group_name}.dat"
    else:
        stats_file = None

    for combo in combinations:
        current_params = dict(zip(keys, combo))
        method_name = get_method_name(job, current_params)

        cmd_args = [str(exec_path)]
        cmd_args.extend(["-i", str(input_file)])

        if 'static_args' in job:
            for k, v in job['static_args'].items():
                cmd_args.append(k)
                if v: cmd_args.append(str(v))

        estimator_params = False
        suffix_parts = []
        for k, v in current_params.items():
            if k in ['-e', '--estimator']:
                estimator_params = True
            cmd_args.extend([k, str(v)])
            clean_k = k.replace('--', '').replace('-', '')
            suffix_parts.append(f"{clean_k}{v}")
        
        if not estimator_params:
            suffix_parts.insert(0, "e" + method_name)
        
        suffix = "_" + "_".join(suffix_parts)

        if stats_file:
            cmd_args.extend(["--output-stats", str(stats_file)])

        if "estimation" in save_config:
            dest_dir = base_job_dir / "estimation" / method_name / subdir
            dest_dir.mkdir(parents=True, exist_ok=True)
            outfile = dest_dir / f"{file_unique_name}{suffix}.pts"
            cmd_args.extend(["-o", str(outfile)])

        if "error" in save_config:
            dest_dir = base_job_dir / "error" / method_name / subdir
            dest_dir.mkdir(parents=True, exist_ok=True)
            outfile = dest_dir / f"{file_unique_name}{suffix}.pts"
            cmd_args.extend(["--output-error", str(outfile)])

        log_file = None
        if verbose:
            log_dir = base_job_dir / "logs" / method_name / subdir
            log_dir.mkdir(parents=True, exist_ok=True)
            log_file = log_dir / f"{file_unique_name}{suffix}.log"
            with open(log_file, "w") as f_log:
                f_log.write(" ".join([str(x) for x in cmd_args]) + "\n\n")

        tasks.append({
            "cmd": cmd_args,
            "log_file": log_file,
            "stats_file": stats_file,
            "verbose": verbose,
            "id": f"{job['name']} | {subdir} | {file_unique_name} | {suffix}"
        })

    return tasks

def run_task(task):
    cmd_str = [str(x) for x in task['cmd']]
    
    if task['verbose'] and task['log_file']:
        with open(task['log_file'], "w") as f_log:
            res = subprocess.run(cmd_str, stdout=f_log, stderr=subprocess.STDOUT, text=True)
    else:
        res = subprocess.run(cmd_str, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
    return res.returncode == 0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True, help='Path to configuration YAML file')
    parser.add_argument('--verbose', action='store_true', help='Enable log files generation')
    args = parser.parse_args()

    config = load_yaml(args.config)
    global_conf = config['global']
    all_tasks = []
    
    print("🛠️  Generating tasks...")
    for job in config['jobs']:
        for ds in config['datasets']:
            files_info = get_files_from_dataset(global_conf, ds)
            if not files_info: continue
            
            for f_info in files_info:
                all_tasks.extend(generate_commands(job, f_info, global_conf, args.verbose))

    if not all_tasks:
        print("Nothing to do.")
        sys.exit(0)

    print("🧹 Cleaning stats files (unique targets)...")
    stats_files_to_clean = set()
    for t in all_tasks:
        if t['stats_file']:
            stats_files_to_clean.add(t['stats_file'])
    
    for f in stats_files_to_clean:
        if f.exists():
            f.unlink()
            
    print(f"🚀 Launching {len(all_tasks)} simulations sequentially...")

    for task in tqdm(all_tasks, unit="sim"):
        success = run_task(task)
        if not success:
            tqdm.write(f"❌ Error on {task['id']} (Enable --verbose to see logs)")

    print("\n✅ Done.")

if __name__ == "__main__":
    main()