#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <set>

#include <Eigen/Dense>
#include "../IO/PointCloudDiff.h"

#include <CLI11.hpp>

void load_normals( PointCloudDiff<Scalar> &pointcloud, std::string input_normals ){
  std::ifstream file( input_normals );
  if( !file.is_open() ){
    spdlog::error( "Could not open file {}", input_normals );
    exit( EXIT_FAILURE );
  }

  int idx = 0;
  std::string line;
  while( std::getline( file, line ) ){
    std::istringstream iss( line );
    Scalar nx, ny, nz;
    iss >> nx >> ny >> nz;
    pointcloud.normals[ idx ] = Eigen::Vector<Scalar, 3>( nx, ny, nz );
    idx++;
  }
  file.close();
}

void load_curvs( PointCloudDiff<Scalar> &pointcloud, std::string input_curvs, Scalar rescale_ratio ){
  std::ifstream file( input_curvs );
  if( !file.is_open() ){
    spdlog::error( "Could not open file {}", input_curvs );
    exit( EXIT_FAILURE );
  }

  int idx = 0;
  std::string line;
  while( std::getline( file, line ) ){
    std::istringstream iss( line );
    Scalar maxCurv, minCurv;
    iss >> maxCurv >> minCurv;
    pointcloud.k1[ idx ] = minCurv / rescale_ratio;
    pointcloud.k2[ idx ] = maxCurv / rescale_ratio;
    pointcloud.gauss[ idx ] = ( minCurv * maxCurv ) / ( rescale_ratio * rescale_ratio );
    pointcloud.mean[ idx ] = ( pointcloud.k1[ idx ] + pointcloud.k2[ idx ] ) / 2.0 ;
    idx++;
  }
  file.close();
}

int main( int argc, char ** argv )
{
  CLI::App app{ "displayPTS" };
  std::string input_xyz;
  app.add_option( "-i,--input", input_xyz, "Input point clouds as .xyz" )->required()->check(CLI::ExistingFile);

  std::string input_normals;
  app.add_option( "-n,--normal", input_normals, "Input point clouds as .normals" )->required()->check(CLI::ExistingFile);

  std::string input_curvs;
  app.add_option( "-c,--curvs", input_curvs, "Input point clouds as .curvs" )->required()->check(CLI::ExistingFile);

  std::string outputFilename = "";
  app.add_option( "-o,--output", outputFilename, "output pointCloud .pts" )->required();

  Scalar ratio = 1.0;
  app.add_option( "-r,--ratio", ratio, "rescale ratio" );

  CLI11_PARSE( app, argc, argv );

  spdlog::info( "Reading pts data..." );

  PointCloudDiff<Scalar> pointCloud;
  pointCloud.loadPointCloud( input_xyz );

  spdlog::info( "Reading normals data..." );
  load_normals( pointCloud, input_normals );

  spdlog::info( "Reading curvs data..." );
  load_curvs( pointCloud, input_curvs, ratio );

  spdlog::info( "Saving the pointCloud..." );

  pointCloud.savePointCloud( outputFilename );

  return EXIT_SUCCESS;
}
