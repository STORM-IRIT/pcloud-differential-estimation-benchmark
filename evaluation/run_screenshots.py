import sys
import yaml
import subprocess
import argparse
import difflib
from pathlib import Path
from tqdm import tqdm

def load_yaml(path):
    if not Path(path).exists():
        print(f"❌ Config not found: {path}")
        sys.exit(1)
    with open(path, 'r') as f:
        return yaml.safe_load(f)

def run_render(cmd, verbose=False):
    cmd_str = [str(x) for x in cmd]
    if verbose:
        print(" ".join(cmd_str))
        subprocess.run(cmd_str, check=True)
    else:
        subprocess.run(cmd_str, stdout=subprocess.DEVNULL, check=True)

def get_clean_shape_name(file_stem):
    """
    Extraction of the "flat" shape name from a complex filename.
    Logic: Take everything before the first numeric part.

    Examples:
    - sphere_25000_ASO      -> sphere
    - goursat_hole_10k_MLS  -> goursat_hole
    - bunny                 -> bunny
    """
    parts = file_stem.split('_')
    clean_parts = []

    for part in parts:
        if part and part[0].isdigit():
            clean_parts.append(part)
            break
        clean_parts.append(part)

    if not clean_parts:
        return file_stem

    return "_".join(clean_parts)

def find_best_camera_match(camera_dir, file_stem):
    """
    Go through all available cameras and find the one whose name is contained in the file name (file_stem).
    We keep the longest match (to distinguish Area_1 from Area_10).
    """
    available_cameras = list(camera_dir.glob("*.json")) + list(camera_dir.glob("*.txt"))
    
    matches = []

    for cam_path in available_cameras:
        cam_name = cam_path.stem
        clean_cam_name = cam_name.replace("cameraSettings_", "").replace("camera_", "")
        
        if clean_cam_name in file_stem:
            matches.append((len(clean_cam_name), cam_path))

    matches.sort(key=lambda x: x[0], reverse=True)

    if matches:
        return matches[0][1]
    
    return None

def resolve_camera_path(camera_dir, cam_config, file_stem):
    """
    Solve the camera path based on the configuration:
    """
    if not cam_config:
        return None

    if cam_config == "auto":
        best_match = find_best_camera_match(camera_dir, file_stem)
        if best_match:
            return best_match
        return None

    candidates = [
        cam_config,
        f"{cam_config}.json",
        f"{cam_config}.txt"
    ]
    
    for cand in candidates:
        p = camera_dir / cand
        if p.exists():
            return p
            
    return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True, help='Path to render_config.yaml')
    parser.add_argument('--verbose', action='store_true', help='Print full commands')
    args = parser.parse_args()

    config = load_yaml(args.config)

    bin_path = Path(config['global']['bin_path'])
    if not bin_path.exists() and bin_path.with_suffix(".exe").exists():
        bin_path = bin_path.with_suffix(".exe")
    if not bin_path.exists():
        print(f"❌ Executable not found: {bin_path}")
        sys.exit(1)

    input_root = Path(config['global']['input_root'])
    output_root = Path(config['global']['output_root'])
    camera_dir = Path(config['global']['camera_dir'])

    jobs = config['jobs']

    print("🚀 Starting Batch Rendering...")

    for job in jobs:
        job_name = job['name']
        job_input_base = input_root / job['input_subdir']

        raw_pattern = job.get('pattern', '*.pts')
        if isinstance(raw_pattern, str):
            patterns = [raw_pattern]
        else:
            patterns = raw_pattern 

        files = []
        for p in patterns:
            print (f"🔍 Searching for pattern '{p}' in {job_input_base}...")
            files.extend(list(job_input_base.glob(p)))

        files = sorted(list(set(files)))

        if not files:
            print(f"⚠️  No files found for job '{job_name}' in {job_input_base} with pattern '{patterns}'")
            continue

        job_output_root = output_root / job_name
        job_output_root.mkdir(parents=True, exist_ok=True)

        is_error = job.get('is_error', False)
        point_radius = job.get('point_radius', 0.005)
        properties = job.get('properties', ["kMean"])

        cameras = job.get('cameras', [])
        if not cameras:
            cameras = [None]

        bounds = job.get('bounds', [-1.0, 1.0])
        diag_crop = job.get('diag_crop', False)
        colorbar = job.get('colorbar', False)
        use_abs = job.get('use_abs', False)
        vector_as_color = job.get('vector_as_color', False)
        fast_render = job.get('fast', False)

        for f in tqdm(files, desc=f"Rendering {job_name}", unit="file"):
            try:
                rel_path = f.relative_to(job_input_base)
                sub_folder = rel_path.parent
            except ValueError:
                sub_folder = Path(".")

            target_dir = job_output_root / sub_folder
            target_dir.mkdir(parents=True, exist_ok=True)

            file_stem = f.stem

            for cam_config in cameras:
                cam_file = resolve_camera_path(camera_dir, cam_config, file_stem)

                if cam_config is not None and cam_file is None:
                    if args.verbose:
                        print(f"⚠️  Camera not found for config '{cam_config}' on shape '{file_stem}'. Skipping.")
                    continue

                for prop in properties:
                    if cam_config == "auto":
                        cam_suffix = "auto"
                    elif cam_config is None:
                        cam_suffix = "default"
                    else:
                        cam_suffix = cam_config

                    out_name = f"{file_stem}_{prop}.png"
                    out_path = target_dir / out_name

                    cmd = [
                        bin_path,
                        "-i", f,
                        "--screenshot", out_path,
                        "--property", prop,
                        "--pointRadius", point_radius,
                        "--minBound", bounds[0],
                        "--maxBound", bounds[1]
                    ]

                    if cam_file:
                        cmd.extend(["--camera", cam_file])

                    if is_error:
                        cmd.append("--errorFile")
                        if use_abs: cmd.append("--abs")

                    if diag_crop: cmd.append("--diagCrop")
                    if colorbar: cmd.append("--colorbar")
                    if vector_as_color: cmd.append("--vectorAsColor")
                    if fast_render: cmd.append("--fast")

                    try:
                        run_render(cmd, args.verbose)
                    except subprocess.CalledProcessError as e:
                        print(f"❌ Error on {f.name}: {e}")

    print("✅ All renders complete.")

if __name__ == "__main__":
    main()