from math import sqrt
import math
import os
import json
import argparse
import re
from pathlib import Path
from collections import defaultdict

SIGNED_CURVATURE = False

METHODS_INDICES = {
    'GroundTruth': 3, 'Mean': 4, 'Cov2D': 5, 'NormCov2D': 6,
    'NormCov3D': 7, 'ShapeOperator': 8, 'PCA': 9, '2-Monge': 10,
    'PC-MLS': 11, 'JetFitting': 12, 'WaveJets': 13, 'Sphere': 14,
    'APSS': 15, 'UnorientedSphere': 16, 'ASO': 17, '3DQuadric': 18,
    'Varifolds': 19, 'AvgHexagram': 20,
}

PROPERTIES_TO_EXPORT = ["kGauss", "kMean"]

EXPORT_TO_HEADER_MAP = {
    "kGauss": "Gauss",
    "kMean": "Mean",
}

RADIUS_REGEX = re.compile(r"[_|-]r(?:adius)?([0-9]+\.[0-9]+)")

radii = ["0.075", "0.1", "0.2"]
output_names = ["CAD", "CAD_helios", "DGtal", "DGtal_helios", "PCPNet"] 
gt_names = ["dataset_CAD", "dataset_helios", "dataset_implicit", "dataset_implicit_helios", "PCPNet_rescaled"]

estimation_header = ["x", "y", "z", "Gauss", "Mean", "nx", "ny", "nz", "k1", "k2", "d1x", "d1y", "d1z", "d2x", "d2y", "d2z"]
error_header = ["x", "y", "z", "iShape", "Gauss", "Mean", "k1", "k2", "d1", "d2", "normal", "pos"]


def parse_file(filepath, type_dir):
    """
    Parses files dynamically based on provided headers.
    """
    coords = []
    props_data = defaultdict(list)
    
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.split()
            if not parts or parts[0].startswith('#') or len(parts) < 3:
                continue
            
            try:
                x, y, z = map(float, parts[:3])
                coords.append((x, y, z))
                
                if type_dir in ["estims", "gt"]:
                    for idx, prop in enumerate(estimation_header[3:], start=3):
                        val = float(parts[idx]) if idx < len(parts) else 0.0
                        # Check if val is inf or NaN, if so set to 0.0
                        if math.isinf(val) or math.isnan(val):
                            val = 0.0
                        props_data[prop].append(val)
                        
                elif type_dir == "errors":
                    for idx, prop in enumerate(error_header[3:], start=3):
                        val = float(parts[idx]) if idx < len(parts) else 0.0
                        # Check if val is inf or NaN, if so set to 0.0
                        if math.isinf(val) or math.isnan(val):
                            val = 0.0
                        props_data[prop].append(sqrt(val))
                        
            except ValueError:
                continue

    if not SIGNED_CURVATURE and "k1" in props_data and "k2" in props_data:            
        props_data["Gauss"] = [abs(kgauss) for kgauss in props_data["Gauss"]]
        props_data["Mean"] = [abs(kmean) for kmean in props_data["Mean"]]
        props_data["k1"] = [min(abs(k1), abs(k2)) for k1, k2 in zip(props_data["k1"], props_data["k2"])]
        props_data["k2"] = [max(abs(k1), abs(k2)) for k1, k2 in zip(props_data["k1"], props_data["k2"])]

    return coords, props_data



def process_all_datasets(input_dir, output_dir, datasets_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    datasets_dir = Path(datasets_dir)
    
    max_col_index = max(METHODS_INDICES.values())
    files_created = 0

    for idx_dataset, output_name in enumerate(output_names):
        gt_dataset_name = gt_names[idx_dataset]
        gt_dir = datasets_dir / gt_dataset_name
        
        gt_index = {}
        if gt_dir.exists():
            for f in gt_dir.rglob("*"):
                if f.suffix in ['.pts', '.xyz']:
                    gt_index[f.stem.lower()] = f
        else:
            print(f"Warning: Ground Truth directory not found: {gt_dir}")
            
        for radius in radii:
            for raw_type_dir, out_type_dir in [("estimation", "estims"), ("error", "errors")]:
                
                base_res_dir = input_dir / output_name / radius / raw_type_dir
                if not base_res_dir.exists():
                    continue
                    
                print(f"Processing: {output_name} | {out_type_dir} | r={radius}")
                
                shapes_map = defaultdict(dict)
                shape_originals = {}
                
                for method_name in METHODS_INDICES.keys():
                    if method_name == 'GroundTruth':
                        continue
                        
                    method_dir = base_res_dir / method_name
                    if not method_dir.exists():
                        continue
                        
                    for filepath in method_dir.rglob("*.pts"):
                        file_stem = filepath.stem.replace("selle", "saddle")
                        
                        shape_original = file_stem
                        if f"_e{method_name}" in shape_original:
                            shape_original = shape_original.split(f"_e{method_name}")[0]
                        else:
                            match = RADIUS_REGEX.search(shape_original)
                            if match: shape_original = shape_original[:match.start()]
                                
                        shape_clean = shape_original
                        if "_" in shape_clean:
                            subparts = shape_clean.rsplit("_", 1)
                            if subparts[1].isdigit(): shape_clean = subparts[0]
                                
                        shapes_map[shape_clean][method_name] = filepath
                        shape_originals[shape_clean] = shape_original

                for shape_clean, methods_files in shapes_map.items():
                    shape_original = shape_originals[shape_clean]
                    shape_data = {}
                    
                    possible_names = [shape_original.lower(), shape_clean.lower()]
                    if "saddle" in shape_original.lower():
                        possible_names.extend([shape_original.lower().replace("saddle", "selle"), shape_clean.lower().replace("saddle", "selle")])
                        
                    gt_filepath = next((gt_index[name] for name in possible_names if name in gt_index), None)
                    if not gt_filepath: 
                        for name in possible_names:
                            for indexed_name, path in gt_index.items():
                                if name in indexed_name:
                                    gt_filepath = path
                                    break
                            if gt_filepath: break
                            
                    if gt_filepath:
                        shape_data['GroundTruth'] = parse_file(gt_filepath, "gt")
                    else:
                        print(f"  [!] Missing Ground Truth for '{shape_original}', skipped.")
                        continue
                        
                    for method, filepath in methods_files.items():
                        shape_data[method] = parse_file(filepath, out_type_dir)
                        
                    coords = shape_data['GroundTruth'][0]
                    num_points = len(coords)
                    
                    for prop_export in PROPERTIES_TO_EXPORT:
                        prop_internal = EXPORT_TO_HEADER_MAP.get(prop_export, prop_export)
                        out_data = []
                        
                        for i in range(num_points):
                            row = [0.0] * (max_col_index + 1)
                            
                            if i < len(coords):
                                row[0], row[1], row[2] = coords[i]
                            else:
                                continue
                            
                                
                            has_valid_data = False
                            
                            for method, method_idx in METHODS_INDICES.items():
                                if method in shape_data:
                                    try:
                                        # If estimating, position error is 0.0
                                        if out_type_dir == "estims" and prop_internal == "pos":
                                            val = 0.0
                                        else:
                                            val = shape_data[method][1][prop_internal][i]
                                            
                                        # Check if val is inf or NaN, if so set to 0.0
                                        if math.isinf(val) or math.isnan(val):
                                            val = 0.0
                                            
                                        if val != 0.0:
                                            has_valid_data = True
                                        row[method_idx] = val
                                    except IndexError:
                                        pass
                                        
                            if has_valid_data:
                                out_data.append(row)
                                
                        if out_data:
                            final_out_dir = output_dir / radius / out_type_dir / output_name / prop_export
                            final_out_dir.mkdir(parents=True, exist_ok=True)
                            
                            with open(final_out_dir / f"{shape_clean}.json", 'w', encoding='utf-8') as f:
                                json.dump(out_data, f, separators=(',', ':'))
                            files_created += 1

    print(f"\nDone! {files_created} JSON files generated.")

def main():
    parser = argparse.ArgumentParser(description="Converts .pts experiment files into structured JSON files.")
    parser.add_argument('dir_in', type=str, help="Input directory (results folder)")
    parser.add_argument('dir_out', type=str, help="Output directory (public/data folder)")
    parser.add_argument('dir_dataset', type=str, help="Datasets directory (for Ground Truth)")

    args = parser.parse_args()
    process_all_datasets(args.dir_in, args.dir_out, args.dir_dataset)

if __name__ == "__main__":
    main()
