#include <iostream>
#include <vector>
#include <chrono>
#include <string>

#include <CLI11.hpp>

#include <DGtal/math/Statistic.h>

#include "../IO/PointCloudDiff.h"

Scalar angularDifference( const Eigen::Vector<Scalar, 3> & n1, const Eigen::Vector<Scalar, 3> & n2 )
{
  Scalar dot_product = n1.dot( n2 );
  dot_product = std::min( Scalar(1.0), std::max( Scalar( -1.0 ), dot_product ) );

  Scalar angle = acos( dot_product ) * Scalar( 180.0 ) / M_PI;
  return std::min( angle, Scalar( 180.0 ) - angle );
}

int main( int argc, char ** argv )
{
  CLI::App app{ "AngularDifference" };

  std::vector<std::string> inputs;
  app.add_option( "-i,--input, 1", inputs, "Input pointcloud (pts)" )->required();

  Scalar noise_normal;
  app.add_option( "--noise-normal", noise_normal, "Standard deviation of the noise on the normal vectors (Gaussian)" );
  
  std::string outputFilename = "";
  app.add_option( "-o,--output", outputFilename, "output curvatures values" );

  seed = 1;
  app.add_option( "--seed", seed, "Seed for the random generator" );

  CLI11_PARSE( app, argc, argv );

  // If type of Scalar is double, print std::cout "double" else "float"
  std::string scalar_type = std::is_same<Scalar, float>::value ? "Scalar type float" : "Scalar type double";
  spdlog::info ( scalar_type );

  spdlog::info( "Measuring the max angular added with the noise..." );

  std::string output_content;

  Scalar errorNormal_Max = 0.0;
  Scalar errorNormal_Min = 0.0;
  Scalar errorNormal_Mean = 0.0;

  for ( auto input : inputs )
  {
    
    PointCloudDiff<Scalar> groundtruth;
    groundtruth.loadPointCloud( input );
    PointCloudDiff<Scalar> estimation( groundtruth );
    estimation.addNoiseNormal( noise_normal );
    Scalar current_errorNormal_Mean = 0.0;
    // eval the angular difference between the normal vectors
    for ( auto i = 0u; i < groundtruth.normals.size(); ++i )
    {
      Scalar err = angularDifference( groundtruth.normals[ i ].normalized(), estimation.normals[ i ].normalized() );
      errorNormal_Max = std::max( err, errorNormal_Max );
      errorNormal_Min = std::min( err, errorNormal_Max );
      current_errorNormal_Mean += err;
    }
    current_errorNormal_Mean /= groundtruth.normals.size();
    errorNormal_Mean += current_errorNormal_Mean;
  }
  errorNormal_Mean /= inputs.size();


  output_content = "Noise normal = " + std::to_string( noise_normal ) + " -> errorNormal : Min = " + std::to_string( errorNormal_Min ) + " Max = " + std::to_string( errorNormal_Max ) + " Mean = " + std::to_string( errorNormal_Mean ) + "\n";
  
  spdlog::info( output_content );


  if ( !outputFilename.empty() )
  {
    std::ofstream ofs( outputFilename, std::ofstream::out | std::ios::app );
    ofs << output_content;
    ofs.close();
  }

  return 0;
}
