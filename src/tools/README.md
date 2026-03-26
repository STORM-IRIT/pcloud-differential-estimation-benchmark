# Point Cloud Utilities & Visualization Tools
This repository contains a collection of utility scripts for manipulating, evaluating, cleaning, and visualizing 3D point clouds and meshes. These tools are designed to work seamlessly with the custom .pts file format used by the estimators.

# Visualization Tools (Polyscope)
These visualizers rely on the polyscope library to render 3D point clouds, with their differential quantities or errors.

## displayPTS2
The primary viewer for .pts point clouds with curvature or error data. 

### Features:

- Renders scalar quantities like kMin, kMax, kMean, kGauss, and Shape Index.
- Displays vectors for Normal, kMin Direction, and kMax Direction.
- Supports generating screenshots (`--screenshot`) with options for automated diagonal cropping (`--diagCrop`).
- Visualizes error using the viridis colormap when the `--errorFile` flag is active (only working on error point clouds).
- Saves and loads custom camera positions and point sizes via JSON (`--camera`).

## displayPTS_mesh
A variant viewer built to load both meshes (.obj, .ply) and point clouds simultaneously.

**Note:** This file contains specific hard-coded translations and rotations applied automatically to point clouds upon loading. It was tailored for a specific, isolated use case.

# Data Manipulation Tools

## `addNormalsFromFile`
This utility merges separate coordinate, normal, and curvature files into a single unified .pts file.

### Parameters:

- `-i, --input`: The input coordinates as a .xyz file.
- `-n, --normal`: The input normals as a .normals file.
- `-c, --curvs`: The input curvatures as a .curvs file.
- `-o, --output`: The target output .pts filename.
- `-r, --ratio`: A rescaling ratio applied to adjust the input curvatures during the merge.

## `rescale_pointCloud`
Scales and centers a point cloud, adjusting all associated differential quantities (Gaussian curvature, mean curvature, principal curvatures) to remain mathematically consistent with the new scale.

### Parameters:

- `-r, --radius` & `-k, --kNN`: Rescales the point cloud so that a given radius contains an average targeted number of kNN neighbors.
- `-R, --ratio`: Forces rescaling using a specific provided ratio.
- `--normalize`: Rescales the point cloud bounding box to fit within [-1, 1]^3.
- `--centering`: Translates the point cloud to center its bounding box at the origin.
- `--transformation`: Applies a 4x4 affine transformation matrix from a text file (compatible with CloudCompare matrix exports).
- `--XYZ`: Restricts the rescaling operation exclusively to the positions, leaving differential quantities untouched.

## `comparePTS`
Checks for differences between two point cloud files, which is useful for validating deterministic outputs across different pipeline runs.

### Parameters:

- `-i, --input`: Takes exactly two point cloud files to compare.
- `-e, --epsilon`: The tolerance threshold for the error comparison (default is 0.0001).

**Behavior:** It compares all 16 columns of the custom format point by point and tallies how many points differ beyond the specified epsilon.

# CAD & Mesh Processing

## `CADprocessing`
Aligns estimated point cloud normals with a ground-truth CAD mesh.

### Parameters:

- `--mesh-input`: The reference mesh file (.ply or .obj).
- `--real-normals`: If flagged, overwrites the point cloud normals directly with the closest normals found on the mesh. If not flagged, it merely reorients the existing point cloud normals to face the same hemisphere as the mesh.
- `--scale, --rotations, --translationX/Y/Z`: Applies geometric transformations to the reference mesh prior to normal extraction.

## meshCleaner
Filters out sparse vertices from a mesh and can overlay ground-truth point clouds.

### Parameters:

- `--mesh-input`: The target mesh to clean.
- `--kmin`: The minimum number of neighbors a vertex must have to be preserved (default is 2).
- `--gui`: Opens a Polyscope visualization window to review the cleaned mesh.

# Evaluation Tools

## `measure_normal_noise`
Evaluates the deviation between ground-truth normals and noisy estimations.

**Behavior:** Computes the angular difference (in degrees) between normal vectors and logs the minimum, maximum, and mean angular error across the dataset.

**Parameters:** Requires the `--noise-normal` standard deviation value to properly log the context of the error.