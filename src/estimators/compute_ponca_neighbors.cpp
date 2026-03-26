#include <iostream>
#include <vector>
#include <chrono>
#include <string>

#include <CLI11.hpp>

#include <DGtal/math/Statistic.h>

#include "../IO/PointCloudDiff.h"
#include "ponca_estimators/estimators.h"


// KdTree_ponca ponca_kdtree;

int main( int argc, char ** argv )
{
  CLI::App app{ "Neighborhood request from point cloud using Ponca" };

  std::string input;
  app.add_option( "-i,--input", input, "Input pointcloud (pts)" )->required()->check( CLI::ExistingFile );

  std::string outputFilename = "";
  app.add_option( "-o,--output", outputFilename, "output curvatures values" )->required();

  radius = 0.5;
  app.add_option( "-r", radius, "Radius for the neighboring points" );

  Scalar noise_position = 0.0;
  app.add_option( "--noise-position", noise_position, "Standard deviation of the noise on the positions (Gaussian)" );
  
  Scalar noise_normal = 0.0;
  app.add_option( "--noise-normal", noise_normal, "Standard deviation of the noise on the normal vectors (Gaussian)" );

  kNN = -1;
  app.add_option( "-k, --kNN", kNN, "Number of neighbors" );

  std::string ponca_weightFunc = "smooth";
  app.add_option( "--kernel", ponca_weightFunc, "Weight function for MLS (smooth, const, wendland, singular)" );

  seed = 0;
  app.add_option( "--seed", seed, "Seed for the random generator" );

  int index = 0;
  app.add_option( "--index", index, "Index of the point to compute the neighborhood request, default 0." );

  std::vector<int> find_radius_for_kNNs;
  app.add_option( "--find-radius", find_radius_for_kNNs, "Find the radius for a given kNN, or list of kNNs" );

  use_kNNGraph = false;
  app.add_flag( "--knn-graph", use_kNNGraph, "Use the kNN graph from ponca, usefull for thin surfaces." );

  mls_iter = 1;

  CLI11_PARSE( app, argc, argv );

  std::string estimator = "dry";

  spdlog::info( "Using " + estimator + " estimator");

  
  if (ponca_weightFunc != "smooth" && 
        ponca_weightFunc != "const" && 
        ponca_weightFunc != "wendland" && 
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

  /////////////////////////////////////////////////////////////////
  // Main Loop
  /////////////////////////////////////////////////////////////////


      /////////////////////////////////////////////////////////////
      // Find a list of radius for a given kNN or list of kNNs if needed
      /////////////////////////////////////////////////////////////


  if ( find_radius_for_kNNs.size() > 0 )
  {
    spdlog::info( "Computing the radius for the given kNNs..." );
    for ( auto k : find_radius_for_kNNs )
    {
      radius = estimateRadiusFromKNN( k );
      radius = 0;
      spdlog::info( "kNN = " + std::to_string( k ) + " -> radius = " + std::to_string( radius ) );
    }
    return 0;
  }

      /////////////////////////////////////////////////////////////
      // Compute the neighborhood request
      /////////////////////////////////////////////////////////////

  spdlog::info( "Computing the neighborhood request using " + std::string(( (kNN == -1)? "radius" : "kNN")) + " search..." );
  
  std::vector<Scalar> neighborhood_value;

  if ( ponca_weightFunc == "const" )
    neighborhood_value = neiRequest<ConstWeightFunc>(index);
  else if ( ponca_weightFunc == "wendland")
    neighborhood_value = neiRequest<WendlandWeightFunc>(index);
  else if ( ponca_weightFunc == "singular")
    neighborhood_value = neiRequest<SingularWeightFunc>(index);
  else
    neighborhood_value = neiRequest<SmoothWeightFunc>(index);

  spdlog::info( "done." );

  /////////////////////////////////////////////////////////////////
  // Export
  /////////////////////////////////////////////////////////////////

  estimator += "_PONCA";

  spdlog::info( "Exporting..." );
  
  std::ofstream ofs( outputFilename, std::ofstream::out );
  ofs << "# x y z nx ny nz nei_val\n";
  std::cout << "Number of points : " << groundtruth.points.size() << std::endl;
  for ( auto i = 0u; i < groundtruth.points.size(); ++i )
  {
    ofs << groundtruth.points[ i ][ 0 ] << " " << groundtruth.points[ i ][ 1 ] << " " << groundtruth.points[ i ][ 2 ];
    ofs << " " << groundtruth.normals[ i ][ 0 ] << " " << groundtruth.normals[ i ][ 1 ] << " " << groundtruth.normals[ i ][ 2 ];
    ofs << " " << neighborhood_value[ i ];
    ofs << std::endl;
  }
  ofs.close();


  return 0;
}
