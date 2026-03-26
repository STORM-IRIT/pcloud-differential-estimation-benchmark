#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <set>

#include <Eigen/Dense>
#include <spdlog/spdlog.h>
#include <CLI11.hpp>

/***
 * @brief Process the two point clouds to compare them.
 * It compares two pts point clouds, to check if they are the same up to a certain epsilon.
 * Usefull to test the output of a point cloud processing, to check two noisy point clouds at a same seed.
 * @param filenames : the two point clouds to compare.
 * @param epsilon : the epsilon for the error.
 * @return void
*/
void processFiles( const std::vector<std::string> & filenames , double epsilon )
{ 

  std::array< std::vector < std::array< double, 16 > >, 2 > pointClouds;
  
  pointClouds[1] = std::vector< std::array< double, 16 > >();

  for (int i = 0; i < filenames.size(); i++) {
    
    std::ifstream ifs( filenames[i] );
    std::string line;
    std::vector < std::array< double , 16 > > pc;
    while ( std::getline( ifs, line ) )
    {
      if ( line.find( "#" ) != std::string::npos )
      {
        spdlog::info( "Skipping # line" );
      }
      else
      {
        std::stringstream linestream( line );
        std::vector<double> data;
        auto i = 0;
        while ( linestream.good() )
        {
          double v;
          linestream >> v;
          data.push_back( v );
          ++i;
        }
        std::array<double, 16> point = {data[0], data[1], data[2], data[3], data[4], data[5], data[6],
                                        data[7], data[8], data[9], data[10], data[11], data[12], data[13],
                                        data[14], data[15]};
        pc.push_back(point);
      }
      pointClouds[i] = pc;
    }
  }
    spdlog::info("pointClouds[0].size() = {}", pointClouds[0].size());
    spdlog::info("pointClouds[1].size() = {}", pointClouds[1].size());

    if (pointClouds[0].size() != pointClouds[1].size()) {
      spdlog::error("The two point clouds must have the same number of points.");
      exit( 1 );
    }

    int nb_diff = 0;
    for (int i = 0; i < pointClouds[0].size(); i++){
        for (int j = 0; j < 16; j++){
            // check if the two points are equal up to epsilon
            if ( std::abs (std::abs(pointClouds[0][i][j]) - std::abs(pointClouds[1][i][j])) > epsilon)
             {
              spdlog::warn("Difference found at point {} coord {} : value {} != {}.", i, j, pointClouds[0][i][j], pointClouds[1][i][j]);
              // exit( 1 );
              nb_diff++;
            }
        }
    }
    spdlog::warn ( "The two point clouds have {} differences.", nb_diff );

}

int main( int argc, char ** argv )
{
  CLI::App app{ "comparePTS" };
  std::vector<std::string> filenames;
  app.add_option( "-i,--input,1", filenames, "2 point clouds to compare." )->required();
  
  double epsilon = 0.0001;
  app.add_option("-e, --epsilon", epsilon, "Epsilon for the error.");

  CLI11_PARSE( app, argc, argv );

  if ( filenames.size() != 2 )
  {
    spdlog::warn("Input files must be at number of 2.");
    exit( 1 );
  }

  // Process comparaison
  processFiles( filenames , epsilon );

  return EXIT_SUCCESS;
}
