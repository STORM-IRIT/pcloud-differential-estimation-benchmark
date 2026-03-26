# Evaluation & Automation Pipeline
This directory contains the automation backbone of the PointCloudDiffBenchmark project. It includes Python scripts and Bash wrappers designed to batch-execute experiments, automate Polyscope screenshot generation, parse results.

Additionally, this directory contains the data_analysis subfolder, which contains all the necessary tools to process the generated statistics and plot comparison graphs.

## Main Execution Scripts (Bash Wrappers)
These bash scripts serve as the primary entry points for running batch operations. They automatically iterate through YAML configuration files stored in specific subdirectories.

- `run_all_experiments.sh`: Scans the `configs_experiments/` directory for YAML files and sequentially executes the benchmarking pipeline.

- `run_all_screenshots.sh`: Scans the `configs_screenshots/` directory for YAML files and automates the batch rendering of point cloud visualizations.

- `run_data_website.sh`: Scans the `configs_data_website/` directory to run a specific set of experiments, and upon success, automatically triggers the conversion of .pts outputs into our web JSON files.

## Core Python Engines
The bash scripts rely on two Python YAML configuration parsers that dynamically build the command-line arguments for the C++ executables.

- `run_experiments.py`
This is the core benchmarking engine.

Parses global settings, datasets, and job definitions from a YAML file.

Dynamically resolves dataset files either by using direct file patterns (e.g., *.pts) or by reading lists of files from text documents.

Generates all combinations of parameters (e.g., different noise levels or radius values).

Executes the C++ estimators and organizes the output .pts files (estimations and errors) and .dat statistics files into structured directories.

- `run_screenshots.py`
This script automates the Polyscope-based viewer (`displayPTS2`) to generate consistent visualizations.

## Utility Scripts
- `convert_pts_to_web.py`: Parses the .pts output files (both estimations and errors) and converts data into optimized JSON arrays. It maps specific properties (like Gaussian or Mean curvature) into a format strictly required for the accompanying website.

- `merge_dat_files.py`: A utility tool that recursively merges two directories of .dat statistic files. It is specifically designed to extract a specific method's data from one directory and cleanly inject/replace it in another, preventing the need to rerun entire benchmarks for a single method update.

- `png_to_jpg.sh:` A quick ImageMagick wrapper (convert) that finds all .png screenshots, flattens them onto a white background, and compresses them into .jpg files at 50% quality to save space.

## Data Analysis
- `data_analysis/`
(Refer to the internal scripts of this folder for more details).
This subfolder contains the scripts to analyse the .dat files and to generate the web dataframes.