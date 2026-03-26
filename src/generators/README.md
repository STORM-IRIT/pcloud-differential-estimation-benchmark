# Point Cloud Generation and Conversion Tools
This repository includes a suite of tools designed to generate synthetic point clouds from implicit polynomial surfaces or convert existing 3D models into a specialized ASCII format. These tools are particularly useful for generating ground-truth data with exact differential quantities (curvatures and normals) for testing estimators.

## Implicit Surface Generators
These tools use the DGtal library to generate point clouds from mathematical implicit surfaces (such as `goursat`, `sphere1`, `ellipsoid`, `cylinder`, `torus`, etc.). They compute and export the exact Gaussian curvature, mean curvature, principal curvatures, principal directions, and normal vectors for each generated point.

## `generateFromImplicit`
This is the standard tool for creating point clouds with ground-truth.

### Core Parameters:

- `-n or --nbPts`: The number of samples to generate.
- `-s or --surface`: The type of implicit polynomial surface to sample (default: goursat).
- `--nbSteps`: The number of steps for the Lloyd relaxation optimization.
- `-r or --rescale`: Rescales the surface bounding box to the range [-1, 1].

### Noise and Outliers:

- `-P or --noise_pos`: Adds Gaussian noise to the point positions.
- `-N or --noise_normal`: Adds Gaussian noise to the exact normal vectors.
- `-F or --noise_flip`: Randomly flips the orientation of the normal vectors based on a given ratio.
- `--outlier_rate` & `--outlier-noise`: Injects a specified ratio of outlier points scattered with a given noise variance.
- `--nb_gen_noise`: Generates multiple noisy variations of the same ground-truth point cloud in a batch.

## `generate_adaptator_pcpnet`
This specialized tool shares the same underlying generation logic but is specifically structured to output data for PCPNet training/evaluation.

- It automatically creates a structured output directory containing a groundtruth/ subfolder with .pts files.
- It exports the noisy point cloud variations strictly as .xyz files containing only coordinates.
- It automatically compiles a dataset.txt file listing all generated .xyz files without their extensions, ensuring compatibility with PCPNet's expected dataset format.
