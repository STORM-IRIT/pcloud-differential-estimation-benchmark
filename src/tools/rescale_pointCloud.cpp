#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <set>

#include <Eigen/Dense>
#include "../IO/PointCloudDiff.h"
#include "../estimators/ponca_estimators/estimators.h"

#include <CLI11.hpp>

bool only_positions=false;


struct TransformationComponents {
    double scale;
    std::vector<double> rotations;
    std::vector<double> translations;
};

TransformationComponents decomposeMatrix(const Eigen::Matrix4d& matrix) {
    TransformationComponents result;

    // Extraire la translation
    result.translations = {matrix(0,3), matrix(1,3), matrix(2,3)};

    // Extraire la matrice de rotation 3x3
    Eigen::Matrix3d rotationMatrix = matrix.block<3,3>(0,0);

    // Calculer le facteur d'échelle
    double scale = rotationMatrix.col(0).norm();
    result.scale = scale;

    // Normaliser la matrice de rotation
    rotationMatrix /= scale;

    // Convertir la matrice de rotation en angles d'Euler
    Eigen::Vector3d eulerAngles = rotationMatrix.eulerAngles(0, 1, 2);

    // Convertir les angles en degrés
    result.rotations = {
        eulerAngles[0] * 180.0 / M_PI,
        eulerAngles[1] * 180.0 / M_PI,
        eulerAngles[2] * 180.0 / M_PI
    };

    return result;
}

Eigen::Matrix4d readMatrixFromFile(const std::string& filename) {
    std::ifstream file(filename);
    Eigen::Matrix4d matrix;

    if (!file.is_open()) {
        throw std::runtime_error("Impossible d'ouvrir le fichier: " + filename);
    }

    std::string line;
    int row = 0;
    while (std::getline(file, line) && row < 4) {
        std::istringstream iss(line);
        for (int col = 0; col < 4; ++col) {
            if (!(iss >> matrix(row, col))) {
                throw std::runtime_error("Erreur de lecture de la matrice dans le fichier");
            }
        }
        ++row;
    }

    file.close();

    if (row != 4) {
        throw std::runtime_error("Le fichier ne contient pas une matrice 4x4 complète");
    }

    return matrix;
}

void applyRotation ( PointCloudDiff<Scalar>& pointCloud, const double& scale, const std::vector <double> &rotations, const std::vector <double> &translations ){

      int nb_points = pointCloud.points.size();
      VectorType rotation    = VectorType (rotations[0], rotations[1], rotations[2]);
      VectorType translation = VectorType(translations[0], translations[1], translations[2]);

      // Rotation (degree) to radian
      rotation = rotation * M_PI / 180.0;

      Eigen::Matrix<Scalar, 3, 3> rot, rotX, rotY, rotZ;
      rotX << 1, 0, 0, 
              0, cos(rotation[0]), -sin(rotation[0]), 
              0, sin(rotation[0]), cos(rotation[0]);
      rotY << cos(rotation[1]), 0, sin(rotation[1]), 
              0, 1, 0, 
              -sin(rotation[1]), 0, cos(rotation[1]);
      rotZ << cos(rotation[2]), -sin(rotation[2]), 0, 
              sin(rotation[2]),  cos(rotation[2]), 0, 
              0, 0, 1;
      rot = rotX.transpose() * rotY.transpose() * rotZ.transpose();

      #pragma omp for
      for (int i = 0; i < nb_points; ++i){
        if (scale != 1.0){
          pointCloud.points[i] = scale * pointCloud.points[i];
        }
        pointCloud.points[i] = rot.transpose() * pointCloud.points[i] + translation;
        pointCloud.normals[i] = rot.transpose() * pointCloud.normals[i];
        pointCloud.normals[i].normalize();
      }
    }

void applyCloudCompareTransformation( PointCloudDiff<Scalar>& pointCloud, const std::string& filename) {
    try {
        Eigen::Matrix4d ccMatrix = readMatrixFromFile(filename);
        TransformationComponents components = decomposeMatrix(ccMatrix);
        
        // Afficher les composants pour vérification
        std::cout << "Scale: " << components.scale << std::endl;
        std::cout << "Rotations: " << components.rotations[0] << ", " 
                  << components.rotations[1] << ", " << components.rotations[2] << std::endl;
        std::cout << "Translations: " << components.translations[0] << ", " 
                  << components.translations[1] << ", " << components.translations[2] << std::endl;
        
        // Appeler votre fonction existante avec les composants décomposés
        applyRotation(pointCloud, components.scale, components.rotations, components.translations);
    } catch (const std::exception& e) {
        std::cerr << "Erreur: " << e.what() << std::endl;
    }
}


Eigen::Vector<Scalar, 3> compute_barycenter ( PointCloudDiff<Scalar>& pointCloud )
{
  Eigen::Vector<Scalar, 3> barycenter = Eigen::Vector<Scalar, 3>::Zero();
  for ( auto i = 0u; i < pointCloud.points.size(); ++i )
  {
    barycenter += pointCloud.points[ i ];
  }
  barycenter /= pointCloud.points.size();

  return barycenter;
}

Scalar compute_max_dist ( PointCloudDiff<Scalar>& pointCloud, const Eigen::Vector<Scalar, 3>& barycenter)
{
  Scalar max_dist = 0.0;
  for ( auto i = 0u; i < pointCloud.points.size(); ++i )
  {
    max_dist = std::max( max_dist, ( pointCloud.points[ i ] - barycenter ).norm() );
  }

  return max_dist;
}

void centeringPointCloud ( PointCloudDiff<Scalar>& pointCloud ) {
  VectorType min = pointCloud.points[0];
  VectorType max = pointCloud.points[0];

  for ( int i = 1; i < pointCloud.points.size(); i++ )
  {
    min = min.cwiseMin( pointCloud.points[i] );
    max = max.cwiseMax( pointCloud.points[i] );
  }

  VectorType center = ( min + max ) / 2.0;
  #pragma omp for
  for ( int i = 0; i < pointCloud.points.size(); i++ )
  {
    pointCloud.points[i] -= center;
  }
}


Scalar compute_rescaling( PointCloudDiff<Scalar>& pointCloud, bool centering, Scalar current_radius, Scalar aimed_radius )
{
  Scalar ratio = aimed_radius / current_radius;
  Scalar max_dist = 0;
  VectorType barycenter = VectorType::Zero();
  if (centering) {
    centeringPointCloud( pointCloud );
  }
  else {
    barycenter = compute_barycenter( pointCloud );
  }
  max_dist = compute_max_dist( pointCloud, barycenter );
  
  for ( int i = 0 ; i < pointCloud.points.size(); ++i)
  {
    pointCloud.points[ i ] = ratio * ( pointCloud.points[ i ] - barycenter );
    if (!only_positions) {
      pointCloud.gauss [ i ] = pointCloud.gauss [ i ] / (ratio * ratio);
      pointCloud.mean  [ i ] = pointCloud.mean [ i ] / ratio;
      pointCloud.k1    [ i ] = pointCloud.k1 [ i ] / ratio;
      pointCloud.k2    [ i ] = pointCloud.k2 [ i ] / ratio;
    }
  }
  return ratio;
}

Scalar compute_rescaling( PointCloudDiff<Scalar>& pointCloud, bool centering, Scalar ratio )
{
  VectorType barycenter = VectorType::Zero();
  if ( centering ){
    centeringPointCloud( pointCloud );
  }
  else { 
    barycenter = compute_barycenter( pointCloud );
  }

  Scalar max_dist = compute_max_dist( pointCloud, barycenter );

  for ( int i = 0 ; i < pointCloud.points.size(); ++i)
  {
    pointCloud.points[ i ] = ratio * ( pointCloud.points[ i ] - barycenter );
    if (!only_positions) {
      pointCloud.gauss [ i ] = pointCloud.gauss [ i ] / (ratio * ratio);
      pointCloud.mean  [ i ] = pointCloud.mean [ i ] / ratio;
      pointCloud.k1    [ i ] = pointCloud.k1 [ i ] / ratio;
      pointCloud.k2    [ i ] = pointCloud.k2 [ i ] / ratio;
    }
  }
  return ratio;
}

Scalar compute_rescaling( PointCloudDiff<Scalar>& pointCloud, bool centering )
{
  VectorType barycenter = VectorType::Zero();
  if ( centering ){
    centeringPointCloud( pointCloud );
  }
  else { 
    barycenter = compute_barycenter( pointCloud );
  }

  Scalar max_dist = compute_max_dist( pointCloud, barycenter );

  Scalar ratio = 1.0 / max_dist;

  for ( int i = 0 ; i < pointCloud.points.size(); ++i)
  {
    pointCloud.points[ i ] = ratio * ( pointCloud.points[ i ] - barycenter );
    if (!only_positions) {
      pointCloud.gauss [ i ] = pointCloud.gauss [ i ] / (ratio * ratio);
      pointCloud.mean  [ i ] = pointCloud.mean [ i ] / ratio;
      pointCloud.k1    [ i ] = pointCloud.k1 [ i ] / ratio;
      pointCloud.k2    [ i ] = pointCloud.k2 [ i ] / ratio;
    }
  }
  return ratio;
}


int main( int argc, char ** argv )
{
  CLI::App app{ "displayPTS" };
  std::string input;
  app.add_option( "-i,--input", input, "Input point clouds" )->required()->check(CLI::ExistingFile);

  std::string outputFilename = "";
  app.add_option( "-o,--output", outputFilename, "output rescaled pointCloud pts" )->required();

  Scalar radius = 0.1;
  app.add_option( "-r, --radius", radius, "Radius for which we want an average of kNN neighbors for each points." );

  int kNN = 50;
  app.add_option( "-k, --kNN", kNN, "Aim number of kNN neighbors for the given radius." );

  Scalar ratio = -1;
  app.add_option( "-R, --ratio", ratio, "Ratio for the rescaling, if given, it uses it." );

  only_positions = false;
  app.add_flag( "--XYZ", only_positions, "Only rescale the positions." );

  std::string outputRatioFile = "";
  app.add_option( "--outputRatio", outputRatioFile, "Output the ratio used for the rescaling." );

  bool normalize = false;
  app.add_flag("--normalize", normalize, "Normalize the point cloud to [-1, 1]^3.");

  bool centering = false;
  app.add_flag("--centering", centering, "Center the point cloud.");

  std::string transformationFile = "";
  app.add_option( "--transformation", transformationFile, "Transformation file to apply to the point cloud." );

  CLI11_PARSE( app, argc, argv );

  spdlog::info( "Reading pts data..." );

  PointCloudDiff<Scalar> pointCloud;
  pointCloud.loadPointCloud( input );

  // PONCA KdTree
  spdlog::info( "Building the ponca KdTree..." );

  buildKdTree<PPAdapter, VectorType>(pointCloud.points, pointCloud.normals, ponca_kdtree);

  spdlog::info( "Find the current radius for the given kNN..." );
  Scalar current_radius = estimateRadiusFromKNN( kNN );
  spdlog::info ( "Current radius for the given kNN: " + std::to_string( current_radius ) );

  if (ratio != -1) {
    spdlog::info( "Rescaling the pointCloud using the ratio " + std::to_string( ratio ) + " ..." );  
    ratio = compute_rescaling( pointCloud, centering, ratio );
    
  }
  else {
    if (normalize){
      spdlog::info( "Rescaling the pointCloud to [-1, 1]^3..." );
      ratio = compute_rescaling( pointCloud, centering );
      spdlog::info ( "Rescaling ratio: " + std::to_string( ratio ) );
    }
    else {
      spdlog::info( "Rescaling the pointCloud using aimed radius..." );
      ratio = compute_rescaling( pointCloud, centering, current_radius, radius );
      spdlog::info ( "Rescaling ratio: " + std::to_string( ratio ) );
    }
  }

  if ( !transformationFile.empty() ) {
    spdlog::info( "Applying transformation from CloudCompare file..." );
    applyCloudCompareTransformation(pointCloud, transformationFile);
  }

  // Save the ratio used for the rescaling
  if ( ! outputRatioFile.empty() )
  {
    // Create / clear the file
    std::ofstream file( outputRatioFile );
    file << ratio;
    file.close();
  }

  spdlog::info( "Saving the rescaled pointCloud..." );
  if ( ! only_positions )
    pointCloud.savePointCloud( outputFilename );
  else {
    pointCloud.savePointCloudAsPositions( outputFilename );
  }

  return EXIT_SUCCESS;
}
