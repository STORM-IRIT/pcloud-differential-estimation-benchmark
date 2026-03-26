#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <chrono>

#include <CLI11.hpp>

#include <DGtal/math/Statistic.h>

#include "../IO/PointCloudDiff.h"
#include "KDTree.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>

typedef Scalar DFT;
typedef CGAL::Simple_cartesian<DFT> Data_Kernel;
typedef Data_Kernel::Point_3 DPoint;
typedef CGAL::Monge_via_jet_fitting<Data_Kernel> My_Monge_via_jet_fitting;
typedef My_Monge_via_jet_fitting::Monge_form My_Monge_form;

int main( int argc, char ** argv )
{
  CLI::App app{ "Differential quantities from point cloud" };

  std::string input;
  app.add_option( "-i,--input", input, "Input pointcloud (pts)" )->required()->check( CLI::ExistingFile );

  std::string outputFilename = "";
  app.add_option( "-o,--output", outputFilename, "output curvatures values" );

  Scalar radius = 1.0;
  app.add_option( "-r", radius, "Radius for the neighboring points (not yet used)" );

  size_t kNN = 0; 
  app.add_option( "-k, --kNN", kNN, "Number of neighbors" );

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
  seed = 0;
  app.add_option( "--seed", seed, "Seed for the random generator" );

  bool computeStats = false;
  app.add_flag( "--stats", computeStats, "Compute error statistics." );

  std::string statsFilename = "";
  app.add_option( "--output-stats", statsFilename, "Output the stats to a txt file (append mode)" );

  std::string outputErrorFilename = "";
  app.add_option( "--output-error", outputErrorFilename, "output errors of curvature values" );
  
  bool useABS = false;
  app.add_flag( "--abs", useABS, "Use absolute errors, default MSE for each points, and RMSE for global stats." );

  bool statsPosition = false;
  app.add_flag( "--stats-position", statsPosition, "Compute the stats on the position" );

  CLI11_PARSE( app, argc, argv );

  // If type of Scalar is double, print std::cout "double" else "float"
  std::string scalar_type = std::is_same<Scalar, float>::value ? "Scalar type float" : "Scalar type double";
  spdlog::info ( scalar_type );

  spdlog::info( "Reading pts and groundtruth data..." );

  PointCloudDiff<Scalar> groundtruth;
  groundtruth.loadPointCloud( input );

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

  KDTree theTree( estimation.points );

  // Jet-fitting settings
  const size_t d_fitting = 4;
  const ssize_t d_monge  = 4;
  const size_t zero = 0;

  Scalar monge_precondition = ( ( d_fitting + 1 ) * ( d_fitting + 2 ) ) / 2;
  spdlog::info( "Monge precondition: {}", monge_precondition );
  spdlog::info( "Computing the differential quantities..." );

  estimation.non_stable_idx = std::vector<unsigned int>( estimation.points.size(), 0 );
  std::vector<Scalar> statNeighbors_values = std::vector<Scalar>( estimation.points.size(), 0 );
  std::vector<Scalar> statTimings_values = std::vector<Scalar>( estimation.points.size(), 0 );


#pragma omp parallel
  {
    My_Monge_form monge_form_thread;
    My_Monge_via_jet_fitting monge_fit_thread;
    std::vector<DPoint> in_points_thread; // CGAL points
    std::vector<std::pair<size_t, Scalar>> indices_thread;

#pragma omp for schedule(dynamic)
    for ( auto id = 0; id < estimation.points.size(); ++id )
    {
      if ( omp_get_thread_num() == 0 && id % 100 == 0 ) {
        DGtal::trace.progressBar( id, estimation.points.size() );
      }

      // Get the neighbors
      auto query = estimation.points[ id ];

      if ( kNN > zero ) {
        indices_thread = theTree.knnSearch( query, kNN );
      } else {
        indices_thread = theTree.radiusSearch( query, radius * radius );
      }

      // Init timer
      std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();

      statNeighbors_values[ id ] = indices_thread.size();

      if ( indices_thread.size() <= monge_precondition ){
        estimation.non_stable_idx[ id ] = 1;
        statTimings_values[ id ] = 0;
        continue;
      }

      in_points_thread.clear();
      in_points_thread.reserve( indices_thread.size() + 1 );
      in_points_thread.push_back( { query[0], query[1], query[2] } );

      bool found_central = false;

      for ( const auto& pair : indices_thread )
      {
        size_t nn = pair.first;
        if ( nn == (size_t)id ) {
          found_central = true;
          continue; // On ne l'ajoute pas ici, il est déjà en position 0
        }

        in_points_thread.push_back( { estimation.points[nn][0], estimation.points[nn][1], estimation.points[nn][2] } );
      }

      monge_form_thread = monge_fit_thread( in_points_thread.begin(), in_points_thread.end(), d_fitting, d_monge );

      Scalar ka  = monge_form_thread.principal_curvatures( 0 );
      Scalar kb  = monge_form_thread.principal_curvatures( 1 );
      auto dmax  = monge_form_thread.maximal_principal_direction();
      auto dmin  = monge_form_thread.minimal_principal_direction();
      auto normal_dir = monge_form_thread.normal_direction();

      // Stop timer

      if ( statsPosition )
        estimation.points[ id ] = { monge_form_thread.origin()[0], monge_form_thread.origin()[1], monge_form_thread.origin()[2] };

      estimation.normals[ id ] = { normal_dir.x(), normal_dir.y(), normal_dir.z() };

      if ( fabs( ka ) < fabs( kb ) )
      {
        estimation.k1[ id ] = fabs( ka );
        estimation.k2[ id ] = fabs( kb );
        estimation.v2[ id ] = { dmin.x(), dmin.y(), dmin.z() };
        estimation.v1[ id ] = { dmax.x(), dmax.y(), dmax.z() };
      }
      else
      {
        estimation.k1[ id ] = fabs( kb );
        estimation.k2[ id ] = fabs( ka );
        estimation.v1[ id ] = { dmin.x(), dmin.y(), dmin.z() };
        estimation.v2[ id ] = { dmax.x(), dmax.y(), dmax.z() };
      }
      estimation.mean[ id ]  = fabs( 0.5 * ( ka + kb ) );
      estimation.gauss[ id ] = fabs( ka * kb );

      std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
      std::chrono::nanoseconds time_span = std::chrono::duration_cast<std::chrono::nanoseconds>( t1 - t0 );
      statTimings_values[ id ] = time_span.count();
    }
  }

  for ( int i = 0; i < statNeighbors_values.size(); i++ )
  {
    Scalar v = statNeighbors_values[ i ];
    Scalar t = statTimings_values[ i ];
    estimation.statNeighbors.addValue( v );
    estimation.statTimings.addValue( t );
  }

  spdlog::info( "done." );

  spdlog::info( "Switching to absolute values for the groundtruth..." );
  groundtruth.switch_to_abs();
  groundtruth.oriented = false;
  estimation.oriented = false;

  if ( outputFilename != "" )
  {
    spdlog::info( "Exporting..." );
    estimation.savePointCloud( outputFilename );
  }

  if ( computeStats )
  {
    spdlog::info( "Computing error statistics..." );
    groundtruth.compare( estimation, useABS );
    estimation.compare( groundtruth, useABS );
    if (outputErrorFilename != "")
    {
      spdlog::info("Exporting point cloud with errors...");
      groundtruth.savePointCloudAsErrors( outputErrorFilename );
    }

    spdlog::info("Average number of neighbors: {}", estimation.statNeighbors.mean());
    
    std::string variant = "JetFitting";

    std::string stats = estimation.getStats ( variant, radius, noise_position, noise_normal, flip_ratio, useABS );
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
