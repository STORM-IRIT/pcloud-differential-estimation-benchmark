#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <chrono>
#include <utility>

#include <CLI11.hpp>

#include <DGtal/math/Statistic.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/vcm_estimate_edges.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_points.h>

#include <Eigen/Dense>
#include <CGAL/Eigen_diagonalize_traits.h>
#include "../IO/PointCloudDiff.h"
#include "KDTree.h"

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

// Point with normal vector stored in a std::pair.
typedef std::pair<Point, Vector> PointVectorPair;
typedef std::vector<PointVectorPair> PointList;

typedef std::array<double,6> Covariance;

int main( int argc, char ** argv )
{
  CLI::App app{ "Differential quantities from point cloud" };
  std::string estimator = "VCM";

  std::string input;
  app.add_option( "-i,--input", input, "Input pointcloud (pts)" )->required()->check( CLI::ExistingFile );

  std::string outputFilename = "";
  app.add_option( "-o,--output", outputFilename, "output curvatures values" );

  Scalar radius = 1.0;
  app.add_option( "-r", radius, "Radius for the neighboring points (not yet used)" );

  Scalar noise_position = 0.0;
  app.add_option( "--noise-position", noise_position, "Standard deviation of the noise on the positions (Gaussian)" );
  Scalar noise_normal = 0.0;
  app.add_option( "--noise-normal", noise_normal, "Standard deviation of the noise on the normal vectors (Gaussian)" );
  Scalar flip_ratio = 0.0;
  app.add_option( "--noise-flip", flip_ratio, "Ratio of points to flip the normal vector" );

  bool computeStats = false;
  app.add_flag( "--stats", computeStats, "Compute error statistics." );

  std::string statsFilename = "";
  app.add_option( "--output-stats", statsFilename, "Output the stats to a txt file (append mode)" );

  std::string outputErrorFilename = "";
  app.add_option( "--output-error", outputErrorFilename, "output errors of curvature values" );

  seed = 0;
  app.add_option( "--seed", seed, "Seed for the random generator" );

  bool unsignedCurvature = false;
  app.add_flag( "--unsigned-curvature", unsignedCurvature, "Compute the curvature as unsigned values" );

  bool useABS = false;
  app.add_flag( "--abs", useABS, "Use absolute errors, default MSE for each points, and RMSE for global stats." );

  bool statsAsEstimations = false;
  app.add_flag( "--stats-as-estimations", statsAsEstimations, "Compute the stats as estimations" );

  Scalar convolution_radius = Scalar(0);
  app.add_option ( "--convolution_radius", convolution_radius, "Convolution radius" );

  CLI11_PARSE( app, argc, argv );

  Scalar offset_radius = radius;
  spdlog::info( "Using " + estimator + " estimator");

  std::string scalar_type = std::is_same<Scalar, float>::value ? "Scalar type float" : "Scalar type double";
  spdlog::info ( scalar_type );

  spdlog::info( "Reading pts and groundtruth data..." );

  PointCloudDiff<Scalar> groundtruth;
  groundtruth.loadPointCloud( input );
  // KDTree theTree( groundtruth.points );

  PointCloudDiff<Scalar> estimation( groundtruth );

  // Add noise if needed
  if ( noise_normal != 0.0 )
    estimation.addNoiseNormal( noise_normal );
  if ( noise_position != 0.0 )
    estimation.addNoisePosition( noise_position );
  if ( flip_ratio != 0.0 )
    estimation.flipNormal( flip_ratio );

  std::vector<PointVectorPair> points( estimation.points.size() );

  // #pragma omp parallel for
  for ( auto id = 0; id < estimation.points.size(); ++id )
  {
    PointVectorPair p;
    p.first = { estimation.points[ id ].x(), estimation.points[ id ].y(), estimation.points[ id ].z() };
    p.second = { estimation.normals[ id ].x(), estimation.normals[ id ].y(), estimation.normals[ id ].z() };
    points[ id ] = p;
  }

  std::vector<Covariance> cov;
  CGAL::First_of_pair_property_map<PointVectorPair> point_map;
  CGAL::compute_vcm(points, cov, offset_radius, convolution_radius,
                    CGAL::parameters::point_map (point_map).geom_traits (Kernel()));

  /////////////////////////////////////////////////////////////////
  // Compute VCM
  spdlog::info( "Computing the differential quantities..." );

  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  
  // #pragma omp parallel for
  for ( auto id = 0; id < estimation.points.size(); ++id ){
    Covariance c = cov[ id ];
    // Do the eigen decomposition of the covariance matrix
    Eigen::Matrix<Scalar, 3, 3> cov_mat;
    cov_mat << c[0], c[1], c[2],
               c[1], c[3], c[4],
               c[2], c[4], c[5];

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, 3, 3>> eigensolver( cov_mat );
    if ( eigensolver.info() != Eigen::Success )
      std::cout << "Eigen decomposition failed" << std::endl;

    Scalar eig_mid = eigensolver.eigenvalues()[ 0 ];
    Scalar eig_min = eigensolver.eigenvalues()[ 1 ];

    // Based on : 
    // Quentin Mérigot, Maks Ovsjanikov, Leonidas J. Guibas. Voronoi-Based Curvature and Feature Estimation from Point Clouds. IEEE Transactions on Visualization and Computer Graphics, 2011, 17 (6), pp.743 - 756. ff10.1109/TVCG.2010.261ff. ffinria-00406575v2

    estimation.k1[ id ]      = std::sqrt( eig_min * (Scalar(4) / ( convolution_radius * convolution_radius ) ) ) ;
    estimation.k2[ id ]      = std::sqrt( eig_mid * (Scalar(4) / ( convolution_radius * convolution_radius ) ) ) ;
    
    estimation.v1[ id ]      = eigensolver.eigenvectors().col( 1 );
    estimation.v2[ id ]      = eigensolver.eigenvectors().col( 0 );
    estimation.normals[ id ] = eigensolver.eigenvectors().col( 2 );
    estimation.mean[ id ]  = Scalar(0.5) * ( estimation.k1[ id ] + estimation.k2[ id ] );
    estimation.gauss[ id ] = estimation.k1[ id ] * estimation.k2[ id ];
  }

  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  auto total_time = ( std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count() );

  // #pragma omp parallel for
  for ( auto id = 0; id < estimation.points.size(); ++id ){
    estimation.statTimings.addValue( total_time / estimation.points.size());
  }

  spdlog::info( "done." );

  spdlog::info( "Switching to absolute values for the groundtruth..." );
  estimation.switch_to_abs();
  groundtruth.switch_to_abs();
  estimation.oriented = false;
  groundtruth.oriented = false;

  spdlog::info( "Checking if k1 < k2 for all the data..." );
  estimation.check_k1_k2();

  if ( outputFilename != "" )
  {
    spdlog::info( "Exporting..." );
    estimation.savePointCloud( outputFilename );
  }

  if ( computeStats )
  {
    spdlog::info( "Computing error statistics..." );

    // Matched with Ponca logic
    groundtruth.compare ( estimation, useABS, statsAsEstimations );
    estimation.compare ( groundtruth, useABS, statsAsEstimations );

    if ( outputErrorFilename != "" )
    {
      spdlog::info("Exporting point cloud with errors...");
      groundtruth.savePointCloudAsErrors( outputErrorFilename );
    }

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
