# Point Cloud Differential Quantities Estimators

These files provide tools to compute differential quantities (such as principal curvatures, mean/Gaussian curvatures, and normal vectors) on point clouds. The primary estimators available are based on CGAL's Jet Fitting and the Ponca library.

## Common Base Parameters

Both main tools share a core set of command-line arguments to handle inputs, outputs, neighborhood definitions, and noise simulation:

### I/O Parameters:

- `-i or --input`: The input point cloud file (required).
- `-o or --output`: The output file to save the computed curvature values.

### Neighborhood Parameters:

- `-r`: The radius to use for neighborhood research.
- `-k or --kNN`: The exact number of nearest neighbors to use.

### Noise & Outlier Simulation:

- `--noise-position`: Standard deviation of Gaussian noise added to point positions.
- `--noise-normal`: Standard deviation of Gaussian noise added to normal vectors.
- `--noise-flip`: Ratio of points where the normal vector orientation is flipped.
- `--outlier-rate`: Ratio of outlier points to inject.
- `--outlier-noise`: Standard deviation of the noise applied to outliers.
- `--seed`: Seed for the random number generator.

### Statistics & Error Metrics:

- `--stats`: Flag to compute error statistics against ground truth data.
- `--output-stats`: Appends the computed statistics to a specified text file.
- `--output-error`: Outputs the point cloud along with curvature error values to a file.
- `--abs`: Use absolute errors instead of the default MSE (Mean Squared Error) for individual points, and RMSE for global stats.
- `--stats-position`: Flag to compute statistics specifically on point positions.

### CGAL Jet Fitting (`compute_jetfitting`)
This tool estimates differential quantities by fitting a Monge patch using jet fitting over the neighborhood of each point, powered by CGAL.

**Usage Note:** It automatically computes a Monge precondition based on the fitting degree (defaulted to d_fitting = 4 and d_monge = 4). Points with insufficient neighbors to satisfy this precondition are marked as unstable.

### Ponca Estimators (`compute_ponca_estimators`)
This tool leverages the Ponca library to estimate differential properties under various methods.

### Specific Parameters
- `-e or --estimator`: Selects the estimator primitive to use. The available options are:
    `ASO`, `PCA`, `MeanPLANE`, `APSS`, `UnorientedSphere`, `Ellipsoid`, `PC-MLS`, `3DQuadric`, `AvgHexagram`, `Hexagram`, `Uniform`, `Independent`, `WaveJets`, `Varifold`, `2-Monge`.
    (Default: `ASO`)

- `--kernel`: Selects the weight function for MLS. Available options:
    `smooth`, `const`, `wendland`, `singular`, `expo`.
    (Default: `smooth`).

- `--mls-iter`: The number of iterations to perform for MLS convergence. (Default: `3`).
- `--knn-graph`: Flag to construct and use Ponca's kNN graph, which is particularly useful for preserving features on thin surfaces.
- `--unsigned-curvature`: Flag to compute curvatures as unsigned absolute values.
- `--idx`: Process and output data for a single point at the given index, useful for debugging a specific neighborhood.
