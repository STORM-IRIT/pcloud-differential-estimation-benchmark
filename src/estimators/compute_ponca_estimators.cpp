#include <iostream>
#include <vector>
#include <chrono>
#include <string>

#include <CLI11.hpp>

#include <DGtal/math/Statistic.h>

#include "../IO/PointCloudDiff.h"
#include "ponca_estimators/estimators.h"


// KdTree_ponca ponca_kdtree;

std::string estimator = "ASO";
std::string ponca_weightFunc = "smooth";

int main( int argc, char ** argv )
{
  CLI::App app{ "Differential quantities from point cloud using Ponca" };

  // app.add_option( "-e,--estimator", estimator, "Estimator to use (ASO, PLANE, APSS, BCylinder, NOCylinder, BOCylinder, FOCylinder, BO2D, FO2D, ellipsoid, ellipsoidV2, CNC, waveJets, orientedWaveJets, varifold; Unoriented_Sphere)" );
  
  estimator = "ASO";
  app.add_option( "-e,--estimator", 
                  estimator, 
                  "Estimator to use (ASO, PCA, MeanPLANE, APSS, UnorientedSphere, Ellipsoid, PC-MLS, OrientedPC-MLS, Direct2-Monge, 3DQuadric, AvgHexagram, Hexagram, Uniform, Independent, WaveJets, OrientedWaveJets, Varifold, VarifoldsMeanPlane, 2-Monge, Oriented2-Monge)" );

  kNN = -1;
  app.add_option( "-k, --kNN", kNN, "Number of neighbors" );

  std::string input;
  app.add_option( "-i,--input", input, "Input pointcloud (pts)" )->required()->check( CLI::ExistingFile );

  std::string outputFilename = "";
  app.add_option( "-o,--output", outputFilename, "output curvatures values" );

  radius = 1.0;
  app.add_option( "-r", radius, "Radius for the neighboring points" );

  Scalar noise_position = 0.0;
  app.add_option( "--noise-position", noise_position, "Standard deviation of the noise on the positions (Gaussian)" );
  
  Scalar noise_normal = 0.0;
  app.add_option( "--noise-normal", noise_normal, "Standard deviation of the noise on the normal vectors (Gaussian)" );

  Scalar flip_ratio = 0.0;
  app.add_option( "--noise-flip", flip_ratio, "Ratio of points to flip the normal vector" );

  Scalar outlier_rate = 0.0;
  app.add_option( "--outlier-rate", outlier_rate, "Ratio of outliers" );

  Scalar outlier_noise = 0.0;
  app.add_option( "--outlier-noise", outlier_noise, "Standard deviation of the noise on the outliers (Gaussian)" );

  mls_iter = 3; 
  app.add_option( "--mls-iter", mls_iter, "Number of iterations for MLS." );
  
  bool computeStats = false;
  app.add_flag( "--stats", computeStats, "Compute error statistics." );

  std::string statsFilename = "";
  app.add_option( "--output-stats", statsFilename, "Output the stats to a txt file (append mode)" );

  std::string outputErrorFilename = "";
  app.add_option( "--output-error", outputErrorFilename, "output errors of curvature values" );

  ponca_weightFunc = "smooth";
  app.add_option( "--kernel", ponca_weightFunc, "Weight function for MLS (smooth, const, wendland, singular, expo)" );

  seed = 0;
  app.add_option( "--seed", seed, "Seed for the random generator" );

  bool useABS = false;
  app.add_flag( "--abs", useABS, "Use absolute errors, default MSE for each points, and RMSE for global stats." );

  use_kNNGraph = false;
  app.add_flag( "--knn-graph", use_kNNGraph, "Use the kNN graph from ponca, usefull for thin surfaces." );

  int uniqueIdx = -1;
  app.add_option( "--idx", uniqueIdx, "Process only one point at the given idx." );

  bool statsPosition = false;
  app.add_flag( "--stats-position", statsPosition, "Compute the stats on the position" );

  bool unsignedCurvature = false;
  app.add_flag( "--unsigned-curvature", unsignedCurvature, "Compute the curvature as unsigned values" );

  bool statsAsEstimations = false;
  app.add_flag( "--stats-as-estimations", statsAsEstimations, "Compute the stats as estimations" );

  CLI11_PARSE( app, argc, argv );

  spdlog::info( "Using " + estimator + " estimator");

  if (ponca_weightFunc != "smooth" && 
        ponca_weightFunc != "const" && 
        ponca_weightFunc != "wendland" && 
        ponca_weightFunc != "expo" &&
        ponca_weightFunc != "singular")
  {
    ponca_weightFunc = "smooth";
    spdlog::warn( "Weight function not recognized, using smooth weight function" );
  }
  spdlog::info( "Using " + ponca_weightFunc + " weight function");

  // If type of Scalar is double, print std::cout "double" else "float"
  std::string scalar_type = std::is_same<Scalar, float>::value ? "Scalar type float" : "Scalar type double";
  spdlog::info ( scalar_type );

  spdlog::info( "Reading pts and groundtruth data..." );

  PointCloudDiff<Scalar> groundtruth;
  groundtruth.loadPointCloud( input );

  // TODO : Check if the "Default" copy the size of Matricies
  PointCloudDiff<Scalar> estimation( groundtruth );

  // Add noise if needed
  if ( noise_normal != 0.0 )
    estimation.addNoiseNormal( noise_normal );
  if ( noise_position != 0.0 )
    estimation.addNoisePosition( noise_position );
  if ( flip_ratio != 0.0 )
    estimation.flipNormal( flip_ratio );
  if ( outlier_rate != 0.0 )
    estimation.addOutliers( outlier_rate, outlier_noise );

  // PONCA KdTree
  spdlog::info( "Building the ponca KdTree..." );

  buildKdTree<PPAdapter, VectorType>(estimation.points, estimation.normals, ponca_kdtree);

  if (use_kNNGraph)
  {
    spdlog::info( "Building the ponca kNN graph..." );
    // delete ponca_knnGraph; // optional as it is set to nullptr by default
    if (kNN == -1){
      ponca_knnGraph = new KnnGraph(*ponca_kdtree);
    }
    else {
      ponca_knnGraph = new KnnGraph(*ponca_kdtree);
    }
  }

  // if kNN is set, estimate the radius from the kNN
  if ( kNN > 0 )
  {
    spdlog::info( "Estimating the max radius corresponding to kNN..." );
    radius = estimateRadiusFromKNN( kNN );
  }

  /////////////////////////////////////////////////////////////////
  // Main Loop
  /////////////////////////////////////////////////////////////////

  spdlog::info( "Computing the differential quantities..." );
  
  DifferentialQuantities<Scalar> diff_quantities;

  if ( ponca_weightFunc == "const" )
    diff_quantities = estimateDifferentialQuantities<ConstWeightFunc>(estimator, uniqueIdx);
  else if ( ponca_weightFunc == "wendland")
    diff_quantities = estimateDifferentialQuantities<WendlandWeightFunc>(estimator, uniqueIdx);
  else if ( ponca_weightFunc == "singular")
    diff_quantities = estimateDifferentialQuantities<SingularWeightFunc>(estimator, uniqueIdx);
  // else if ( ponca_weightFunc == "expo")
    // diff_quantities = estimateDifferentialQuantities<ExponentialWeightFunc>(estimator, uniqueIdx);
  else
    diff_quantities = estimateDifferentialQuantities<SmoothWeightFunc>(estimator, uniqueIdx);

  estimation.k1    = diff_quantities.k1();
  estimation.k2    = diff_quantities.k2();
  estimation.mean  = diff_quantities.mean();
  estimation.gauss = diff_quantities.gauss();

  estimation.v1    = diff_quantities.d1();
  estimation.v2    = diff_quantities.d2();
  estimation.normals = diff_quantities.normal();

  estimation.statNeighbors = diff_quantities.statNeighbors();
  estimation.statTimings   = diff_quantities.statTimings();

  if ( statsPosition )
    estimation.points = diff_quantities.position();

  Scalar non_stable_ratio = diff_quantities.getNonStableRatio();
  estimation.non_stable_idx = diff_quantities.getNonStableVector();

  if ( unsignedCurvature || ! diff_quantities.isOriented() ){
    spdlog::info( "Switching to absolute values..." );
    estimation.switch_to_abs();
    groundtruth.switch_to_abs();
    estimation.oriented = false;
    groundtruth.oriented = false;
  }
  
  spdlog::info( "Checking if k1 < k2 for all the data..." );
  estimation.check_k1_k2();

  spdlog::info( "done." );

  /////////////////////////////////////////////////////////////////
  // Export
  /////////////////////////////////////////////////////////////////
  
  
  // if kNN is set, make the radius being the kNN to facilitate the comparison
  if ( kNN > 0 )
  {
    spdlog::info( "Puting to radius the value of kNN..." );
    radius = kNN;
  }

  // estimator += "_PONCA";

  if ( uniqueIdx > -1 ) {
    spdlog::info( "Exporting the differential quantities for the unique point..." );
    estimation.compute_unique_error( groundtruth, uniqueIdx, useABS );
    std::string stats = estimation.getStats( estimator, radius, noise_position, noise_normal, flip_ratio, useABS );
    std::cout << stats;

    if ( statsFilename != "" )
    {
      std::ofstream ofs( statsFilename, std::ios::out | std::ios::app );
      ofs << stats;
      ofs.close();
    }
    return 0;
  }

  if ( outputFilename != "" )
  {
    spdlog::info( "Exporting..." );
    estimation.savePointCloud( outputFilename );
  }

  if ( computeStats )
  {
    spdlog::info( "Computing error statistics..." );
    groundtruth.compare ( estimation, useABS, statsAsEstimations );
    estimation.compare ( groundtruth, useABS, statsAsEstimations );

    if ( outputErrorFilename != "" )
    {
      spdlog::info("Exporting point cloud with errors...");
      groundtruth.savePointCloudAsErrors( outputErrorFilename );
    }

    // estimation.compare( groundtruth );
    spdlog::info("Average number of neighbors: {}", estimation.statNeighbors.mean());

    std::string stats = estimation.getStats ( estimator, radius, noise_position, noise_normal, flip_ratio, useABS );
    std::cout << stats;

    if ( statsFilename != "" )
    {
      std::ofstream ofs( statsFilename, std::ios::out | std::ios::app );
      ofs << stats;
      ofs.close();
    }

  }

  return 0;
}
