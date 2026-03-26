import os
import argparse

def process_dat_files(file1_path, file2_path, method_name, output_path):
    """
    Processes two individual files and saves the result.
    """
    kept_lines_1 = []
    with open(file1_path, 'r', encoding='utf-8') as f1:
        for line in f1:
            if not line.strip():
                continue
            if line.strip().split()[-1] != method_name:
                kept_lines_1.append(line)

    kept_lines_2 = []
    with open(file2_path, 'r', encoding='utf-8') as f2:
        for line in f2:
            if not line.strip():
                continue
            if line.strip().split()[-1] == method_name:
                kept_lines_2.append(line)

    combined_lines = kept_lines_1 + kept_lines_2

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, 'w', encoding='utf-8') as out_f:
        out_f.writelines(combined_lines)

def process_directories(dir1, dir2, method_name, out_dir):
    """
    Iterates through dir1, finds matches in dir2, and saves the output to out_dir.
    """
    processed_count = 0

    for root, dirs, files in os.walk(dir1):
        for file in files:
            if file.endswith(".dat"):
                file1_path = os.path.join(root, file)

                relative_path = os.path.relpath(file1_path, dir1)

                file2_path = os.path.join(dir2, relative_path)
                output_path = os.path.join(out_dir, relative_path)

                if not os.path.exists(file2_path):
                    print(f"Warning: File {file2_path} was not found. File skipped.")
                    continue

                process_dat_files(file1_path, file2_path, method_name, output_path)
                processed_count += 1
                print(f"Processed: {relative_path}")

    print(f"\nDone! {processed_count} files were successfully merged into '{out_dir}'.")

def main():
    parser = argparse.ArgumentParser(description="Recursively merges method data across directories of .dat files.")
    parser.add_argument('dir1', type=str, help="Directory 1 (contains files from which the method will be removed)")
    parser.add_argument('dir2', type=str, help="Directory 2 (contains files from which the method will be extracted)")
    parser.add_argument('method_name', type=str, help="Name of the method to replace (e.g., JetFitting)")
    parser.add_argument('out_dir', type=str, help="Output directory (will contain the new directory tree)")

    args = parser.parse_args()

    process_directories(args.dir1, args.dir2, args.method_name, args.out_dir)

if __name__ == "__main__":
    main()