#include <iostream>
#include <vector>
#include <chrono>
#include <string>

#include <CLI11.hpp>

#include <DGtal/math/Statistic.h>

#include "../IO/PointCloudDiff.h"

/***
 * @brief Read the estimations from PCPNet for the curv estimation
 * It is stored as a .curv file
 * The format is the following : 
 * For each line corresponding to the same line in the point cloud file : K1 K2 
 */
void read_estimations_PCPNet_curv( const std::string & filename, PointCloudDiff<Scalar> & estimation )
{
  std::ifstream ifs( filename );
  std::string line;
  Scalar k1, k2;
  int idx = 0;
  while ( std::getline( ifs, line ) )
  {
    std::stringstream linestream( line );
    linestream >> k1 >> k2;
    if ( k1 > k2 )
    {
      estimation.k1[ idx ]  = k2;
      estimation.k2[ idx ]  = k1;
    }
    else
    {
      estimation.k1[ idx ]  = k1;
      estimation.k2[ idx ]  = k2;
    }
    estimation.mean[ idx ]  = ( estimation.k1[ idx ] + estimation.k2[ idx ] ) / 2.0;
    estimation.gauss[ idx ] = estimation.k1[ idx ] * estimation.k2[ idx ];
    idx++;
  }
  ifs.close();
  spdlog::info( "Read: {} points with the .curv estimation file", idx);
  if ( idx != estimation.points.size() )
  {
    spdlog::error( "Number of points in the estimation file does not match the number of points in the point cloud" );
    exit( 1 );
  }
}

/***
 * @brief Read the estimations from PCPNet for the normal estimation
 * It is stored as a .normal file
 * The format is the following : 
 * For each line corresponding to the same line in the point cloud file : nx ny nz 
 */
void read_estimations_PCPNet_normal( const std::string & filename, PointCloudDiff<Scalar> & estimation )
{
  std::ifstream ifs( filename );
  std::string line;
  Scalar nx, ny, nz;
  int idx = 0;
  while ( std::getline( ifs, line ) )
  {
    std::stringstream linestream( line );
    linestream >> nx >> ny >> nz;
    Eigen::Vector<Scalar, 3> normal( nx, ny, nz );
    estimation.normals[ idx ] = normal;
    idx++;
  }
  ifs.close();
  spdlog::info( "Read: {} points with the .normals estimation file", idx);
  if ( idx != estimation.points.size() )
  {
    spdlog::error( "Number of points in the estimation file does not match the number of points in the point cloud" );
    exit( 1 );
  }
}

/**
 * @brief Read the noise parameters from the filename
 * formated as [surfaceName]_[nbPts]_[noisePos]_[noiseNorm]_[noiseFlip].[extension] 
 * or [surfaceName]_[nbPts].[extension] for 0 noise
 */
void read_noise_parameters (const std::string & filename, Scalar & n_position, Scalar & n_normal, Scalar & f_ratio)
{
  std::string delimiter = "_";
  std::vector <std::string> tokens;
  std::string noise_position, noise_normal, flip_ratio;
  size_t pos = 0;
  std::string token;
  std::string filename_copy = filename;
  
  while ((pos = filename_copy.find(delimiter)) != std::string::npos) {
    token = filename_copy.substr(0, pos);
    tokens.push_back(token);
    filename_copy.erase(0, pos + delimiter.length());
  }
  tokens.push_back(filename_copy);
  
  if (tokens.size() == 2){
    // No noise parameters
    n_position = Scalar( 0 );
    n_normal = Scalar( 0 );
    f_ratio = Scalar( 0 );
    return;
  }

  for (auto & token : tokens){
    spdlog::info(token);
  }

  noise_position = tokens[2];
  noise_normal = tokens[3];
  flip_ratio = tokens[4];

  // noise_position / noise_normal and flip_ratio are formated like "0.0 or 0,0", we need to convert coma to dot
  std::replace(noise_position.begin(), noise_position.end(), ',', '.');
  std::replace(noise_normal.begin(), noise_normal.end(), ',', '.');
  std::replace(flip_ratio.begin(), flip_ratio.end(), ',', '.');

  if ( std::is_same<Scalar, double>::value ) {
    n_position = std::stod(noise_position);
    n_normal = std::stod(noise_normal);
    f_ratio = std::stod(flip_ratio);
  }
  else {
    n_position = std::stof(noise_position);
    n_normal = std::stof(noise_normal);
    f_ratio = std::stof(flip_ratio);
  }
}

int main( int argc, char ** argv )
{
  CLI::App app{ "Differential quantities from point cloud by reading PCPNet outputs" };

  std::string estimator = "PCPNet";
  app.add_option( "--estimator", estimator, "Estimator to use (PCPNet)" );

  std::string input;
  app.add_option( "-i,--input", input, "Input pointcloud (pts)" )->required()->check( CLI::ExistingFile );

  std::string inputEstimation;
  app.add_option( "-e,--estimation-file", inputEstimation, "Input PCPNet estimation (formated as name_nbPts_[noisePos]_[noiseNorm]_[noiseFlip][.extension])" )->required();

  std::string outputFilename = "";
  app.add_option( "-o,--output", outputFilename, "output curvatures values" );

  bool computeStats = false;
  app.add_flag( "--stats", computeStats, "Compute error statistics." );

  std::string statsFilename = "";
  app.add_option( "--output-stats", statsFilename, "Output the stats to a txt file (append mode)" );

  std::string outputErrorFilename = "";
  app.add_option( "--output-error", outputErrorFilename, "output errors of curvature values" );

  bool unsignedCurvature = false;
  app.add_flag( "--unsigned-curvature", unsignedCurvature, "Compute the curvature as unsigned values" );

  bool useABS = false;
  app.add_flag( "--abs", useABS, "Use absolute errors, default MSE for each points, and RMSE for global stats." );

  bool statsAsEstimations = false;
  app.add_flag( "--stats-as-estimations", statsAsEstimations, "Compute the stats as estimations" );
  CLI11_PARSE( app, argc, argv );

  if ( estimator != "PCPNet" )
  {
    spdlog::warn( "Default estimator : PCPNet" );
    estimator = "PCPNet";
  }

  /////////////////////////////////////////////////////////////////
  // Check the inputEstimation PCPNet estimation
  /////////////////////////////////////////////////////////////////

  // Separate the path and the filename
  std::string path, filename;
  std::string delimiter = "/";
  size_t pos = 0;
  std::string token;
  std::string inputEstimation_copy = inputEstimation;
  while ((pos = inputEstimation_copy.find(delimiter)) != std::string::npos) {
    token = inputEstimation_copy.substr(0, pos);
    path = path + token + "/";
    inputEstimation_copy.erase(0, pos + delimiter.length());
  }
  filename = inputEstimation_copy.substr(0, inputEstimation_copy.find_last_of("."));
  inputEstimation = path + filename;

  std::ifstream file(inputEstimation + ".curv");
  if ( ! file.good() ){
    spdlog::error("File {} does not exist", inputEstimation);
    return 1;    
  }

  std::ifstream file2(inputEstimation + ".normals");
  if ( ! file2.good() ){
    spdlog::error("File {} does not exist", inputEstimation);
    return 1;
  }

  spdlog::info( "Reading noise parameters..." );

  Scalar noise_position, noise_normal, flip_ratio;
  read_noise_parameters(filename, noise_position, noise_normal, flip_ratio);

  spdlog::info( "Noise position: {}", noise_position );
  spdlog::info( "Noise normal: {}", noise_normal );
  spdlog::info( "Flip ratio: {}", flip_ratio );

  spdlog::info( "Checking " + estimator + " estimation" );

  // If type of Scalar is double, print std::cout "double" else "float"
  std::string scalar_type = std::is_same<Scalar, float>::value ? "Scalar type float" : "Scalar type double";
  
  spdlog::info ( scalar_type );

  spdlog::info( "Reading pts and groundtruth data..." );

  PointCloudDiff<Scalar> groundtruth;
  groundtruth.loadPointCloud( input );

  // TODO : Check if the "Default" copy the size of Matricies
  PointCloudDiff<Scalar> estimation( groundtruth );

  /////////////////////////////////////////////////////////////////
  // Main Loop
  /////////////////////////////////////////////////////////////////

  spdlog::info( "Reading the differential quantities estimation..." );

  read_estimations_PCPNet_curv( inputEstimation + ".curv", estimation );

  read_estimations_PCPNet_normal( inputEstimation + ".normals", estimation );

  std::vector<Eigen::Vector<Scalar, 3>> zero(groundtruth.points.size(), Eigen::Vector<Scalar, 3>(0, 0, 0));

  estimation.v1    = zero;
  estimation.v2    = zero;

  // estimation.statNeighbors = diff_quantities.statNeighbors();
  // estimation.statTimings   = diff_quantities.statTimings();
  // estimation.points = diff_quantities.position();

  spdlog::info( "done." );

  /////////////////////////////////////////////////////////////////
  // Export
  /////////////////////////////////////////////////////////////////

  spdlog::info( "Switching to absolute values..." );
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
    groundtruth.compare ( estimation, useABS, statsAsEstimations );
    estimation.compare ( groundtruth, useABS, statsAsEstimations );

    if (outputErrorFilename != "")
    {
      spdlog::info("Exporting point cloud with errors...");
      groundtruth.savePointCloudAsErrors( outputErrorFilename );
    }

    spdlog::info("Average number of neighbors: {}", estimation.statNeighbors.mean());

    Scalar radius = 0.0;

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
