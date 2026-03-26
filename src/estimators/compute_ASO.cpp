#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <string>

#include <CLI11.hpp>

#include <DGtal/math/Statistic.h>

#include "../IO/PointCloudDiff.h"
#include "KDTree.h"

#include "ASO/AlgebraicShapeOperator.h"

int main( int argc, char ** argv )
{
  CLI::App app{ "Differential quantities from point cloud" };

  std::string input;
  app.add_option( "-i,--input", input, "Input pointcloud (pts)" )->required()->check( CLI::ExistingFile );

  std::string outputFilename = "";
  app.add_option( "-o,--output", outputFilename, "output curvatures values" );

  Scalar radius = 1.0;
  app.add_option( "-r", radius, "Radius for the neighboring points" );

  Scalar noise_position = 0.0;
  app.add_option( "--noise-position", noise_position, "Standard deviation of the noise on the positions (Gaussian)" );

  Scalar noise_normal = 0.0;
  app.add_option( "--noise-normal", noise_normal, "Standard deviation of the noise on the normal vectors (Gaussian)" );

  bool computeStats = false;
  app.add_flag( "--stats", computeStats, "Compute error statistics." );

  std::string statsFilename = "";
  app.add_option( "--output-stats", statsFilename, "Output the stats to a txt file (append mode)" );

  std::string outputErrorFilename = "";
  app.add_option( "--output-error", outputErrorFilename, "output errors of curvature values" );

  seed = 0;
  app.add_option( "--seed", seed, "Seed for the random generator" );
  
  CLI11_PARSE( app, argc, argv );

  spdlog::info( "Reading pts and groundtruth data..." );

  PointCloudDiff<Scalar> groundtruth;
  groundtruth.loadPointCloud( input );

  PointCloudDiff<Scalar> estimation( groundtruth );

  // Add noise if needed
  if ( noise_normal != 0.0 )
    estimation.addNoiseNormal( noise_normal );
  if ( noise_position != 0.0 )
    estimation.addNoisePosition( noise_position );

  KDTree<Scalar> theTree( estimation.points );
  
  /////////////////////////////////////////////////////////////////
  // Main Loop
  spdlog::info( "Computing the differential quantities..." );
  for ( auto id = 0; id < estimation.points.size(); ++id )
  {
    DGtal::trace.progressBar( id, estimation.points.size() );

    // Get the knn
    auto query   = estimation.points[ id ];
    auto indices = theTree.radiusSearch( query, radius * radius );
    // auto indices = theTree.knnSearch( query, radius );

    estimation.statNeighbors.addValue( indices.size() );

    std::vector<typename PointCloudDiff<Scalar>::Point> neigPoints( indices.size() ), neigNormals( indices.size() );
    for ( auto j = 0; j < indices.size(); ++j )
    {
      const auto nid   = indices[ j ].first; // Especially for the radius search
      // const auto nid   = indices[ j ]; // Especially for the kNN search
      neigPoints[ j ]  = estimation.points[ nid ];
      neigNormals[ j ] = estimation.normals[ nid ];
    }

    // ASO Main
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    auto resAlgebraic = aso::compute( query, radius, neigPoints.begin(), neigPoints.end(), neigNormals.begin(), neigNormals.end() );

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    estimation.statTimings.addValue( std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count() );

    estimation.mean[ id ]  = resAlgebraic.H();
    estimation.gauss[ id ] = resAlgebraic.K();
    Scalar asok1 = resAlgebraic.k1(), asok2 = resAlgebraic.k2();
    auto dd1 = resAlgebraic.d1(), dd2 = resAlgebraic.d2();
    if ( asok1 >= asok2 )
    {
      std::swap( asok1, asok2 );
      std::swap( dd1, dd2 );
    }
    estimation.k1[ id ] = asok1;
    estimation.k2[ id ] = asok2;
    estimation.v1[ id ] = { dd1( 0 ), dd1( 1 ), dd1( 2 ) };
    estimation.v2[ id ] = { dd2( 0 ), dd2( 1 ), dd2( 2 ) };
  }
  spdlog::info( "done." );

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
    estimation.compare( groundtruth, -1 );

    if (outputErrorFilename != "")
    {
      spdlog::info("Exporting point cloud with errors...");
      estimation.savePointCloudAsErrors( outputErrorFilename );
    }
    
    spdlog::info("Average number of neighbors: {}", estimation.statNeighbors.mean());

    std::cout << "# nbPoints radius noise-position noise-normal flip-normal {avg,max,variance} for [nbNeighbors,mean,gauss,k1,k2,d1,d2,timings,pos,idxShape,normal] non_stable_ratio and variant-name\n";
    std::string variant = "ASO";
    std::cout << groundtruth.points.size() << " " << radius << " " << noise_position << " " << noise_normal << " " << 0.0 << " " 
              << estimation.statNeighbors.mean() << " " << estimation.statNeighbors.max() << " " << estimation.statNeighbors.unbiasedVariance() << " " 
              << estimation.errorMean.mean() << " " << estimation.errorMean.max() << " " << estimation.errorMean.unbiasedVariance() << " " 
              << estimation.errorGauss.mean() << " " << estimation.errorGauss.max() << " " << estimation.errorGauss.unbiasedVariance() << " " 
              << estimation.errorK1.mean() << " " << estimation.errorK1.max() << " " << estimation.errorK1.unbiasedVariance() << " " 
              << estimation.errorK2.mean() << " " << estimation.errorK2.max() << " " << estimation.errorK2.unbiasedVariance() << " " 
              << estimation.errorD1.mean() << " " << estimation.errorD1.max() << " " << estimation.errorD1.unbiasedVariance() << " " 
              << estimation.errorD2.mean() << " " << estimation.errorD2.max() << " " << estimation.errorD2.unbiasedVariance() << " " 
              << estimation.statTimings.mean() << " " << estimation.statTimings.max() << " " << estimation.statTimings.unbiasedVariance() << " "
              << estimation.errorPos.mean() << " " << estimation.errorPos.max() << " " << estimation.errorPos.unbiasedVariance() << " " 
              << estimation.errorShapeIndex.mean() << " " << estimation.errorShapeIndex.max() << " " << estimation.errorShapeIndex.unbiasedVariance() << " " 
              << estimation.errorNormal.mean() << " " << estimation.errorNormal.max() << " " << estimation.errorNormal.unbiasedVariance() << " " 
              << 0 << " " << variant << std::endl;

    if ( statsFilename != "" )
    {
      std::ofstream ofs( statsFilename, std::ios::out | std::ios::app );

      ofs << groundtruth.points.size() << " " << radius << " " << noise_position << " " << noise_normal << " " << 0.0 << " "
          << estimation.statNeighbors.mean() << " " << estimation.statNeighbors.max() << " " << estimation.statNeighbors.unbiasedVariance() << " " 
          << estimation.errorMean.mean() << " " << estimation.errorMean.max() << " " << estimation.errorMean.unbiasedVariance() << " " 
          << estimation.errorGauss.mean() << " " << estimation.errorGauss.max() << " " << estimation.errorGauss.unbiasedVariance() << " " 
          << estimation.errorK1.mean() << " " << estimation.errorK1.max() << " " << estimation.errorK1.unbiasedVariance() << " " 
          << estimation.errorK2.mean() << " " << estimation.errorK2.max() << " " << estimation.errorK2.unbiasedVariance() << " " 
          << estimation.errorD1.mean() << " " << estimation.errorD1.max() << " " << estimation.errorD1.unbiasedVariance() << " " 
          << estimation.errorD2.mean() << " " << estimation.errorD2.max() << " " << estimation.errorD2.unbiasedVariance() << " " 
          << estimation.statTimings.mean() << " " << estimation.statTimings.max() << " " << estimation.statTimings.unbiasedVariance() << " "
          << estimation.errorPos.mean() << " " << estimation.errorPos.max() << " " << estimation.errorPos.unbiasedVariance() << " "
          << estimation.errorShapeIndex.mean() << " " << estimation.errorShapeIndex.max() << " " << estimation.errorShapeIndex.unbiasedVariance() << " " 
          << estimation.errorNormal.mean() << " " << estimation.errorNormal.max() << " " << estimation.errorNormal.unbiasedVariance() << " " 
          << 0 << " " << variant << std::endl;

      ofs.close();
    }

  }
  return 0;
}
