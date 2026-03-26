#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include <set>

#include <igl/read_triangle_mesh.h>
#include <igl/readPLY.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>

// Includes to process PCA for normal estimation (without orientation)
#include "../IO/PointCloudDiff.h"
#include "../estimators/ponca_estimators/estimators.h"

#include <CLI11.hpp>

typedef Eigen::Vector<Scalar, 3>    VectorType;
typedef Eigen::Matrix<Scalar, 3, 3> MatrixType;

typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>  SampleMatrixType;
typedef Eigen::Vector<Scalar, Eigen::Dynamic>                  SampleVectorType;

class Mesh {

  public :
    SampleMatrixType V;
    SampleMatrixType N;
    Eigen::MatrixXi F;

  public :

    SampleMatrixType getV() const { return V; }
    SampleMatrixType getN() const { return N; }
    std::vector<VectorType> getVVector() const {
      std::vector<VectorType> V_vector(V.rows());
      for (int i = 0; i < V.rows(); ++i)
        V_vector[i] = V.row(i).transpose();
      return V_vector;
    }
    std::vector<VectorType> getNVector() const {
      std::vector<VectorType> N_vector(N.rows());
      for (int i = 0; i < N.rows(); ++i)
        N_vector[i] = N.row(i).transpose();
      return N_vector;
    }

  public :

    Mesh () {}

    Mesh (const SampleMatrixType & _V, const SampleMatrixType & _N, const Eigen::MatrixXi & _F) : V(_V), N(_N), F(_F) {}

    Mesh (const std::vector<VectorType> &_V, const std::vector<VectorType> &_N){
      V.resize(_V.size(), 3);
      N.resize(_N.size(), 3);
      for (int i = 0; i < _V.size(); ++i){
          V.row(i) = _V[i].transpose();
          N.row(i) = _N[i].transpose();
      }
    }


    void applyRotation ( const double& scale, const std::vector <double> &rotations, const std::vector <double> &translations ){

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
      for (int i = 0; i < V.rows(); ++i){
        if (scale != 1.0){
          V.row(i) = scale * V.row(i);
        }
        V.row(i) = V.row(i) * rot + translation.transpose();
        N.row(i) = N.row(i) * rot;
        N.row(i).normalize();
      }
    }

    const VectorType closestNormal (const VectorType & point) const {
      VectorType normal_i = VectorType::Zero();
      double minDist = std::numeric_limits<double>::max();
      for (int j = 0; j < V.rows(); ++j) {
        VectorType v_point = V.row(j).transpose();
        double dist = (point - v_point).norm();
        if (dist < minDist){
            minDist = dist;
            normal_i = N.row(j);
        }
      }
      return normal_i;
    }
};


Mesh loadMesh( const std::string & filename, bool per_face=true )
{
  SampleMatrixType V;
  SampleMatrixType N;
  Eigen::MatrixXi F;
  

  if (filename.substr(filename.find_last_of(".") + 1) == "ply"){
      Eigen::MatrixXi cloudE;
      SampleMatrixType cloudUV;
      igl::readPLY(filename, V, F, cloudE, N, cloudUV);
  }
  else {
    igl::read_triangle_mesh(filename, V, F);
  }
  if (N.rows() == 0) {
      if (per_face) {
        igl::per_face_normals(V, F, N);
      } else { 
        igl::per_vertex_normals(V, F, N);
      }
  }
  if ( F.rows() != 0 && N.rows() == 0 ) {
      std::cerr << "[libIGL] An error occurred when computing the normal vectors from the mesh. Aborting..."
              << std::endl;
      exit (EXIT_FAILURE);
  }

  for (int i = 0; i < N.rows(); ++i)
      N.row(i).normalize();

  if ( per_face ){
    SampleMatrixType V_faces(F.rows(), 3);
    for (int i = 0; i < F.rows(); ++i){
        V_faces.row(i) = ( V.row( F( i, 0 ) ) + V.row( F( i, 1 ) ) + V.row( F( i, 2 ) ) ) / 3.0;
    }
    return Mesh(V_faces, N, F);
  }

  return Mesh(V, N, F);
}

/**
 * @brief It find the correct normal of each point, 
 * giving a mesh file corresponding to the captures point cloud.
 * @param[pointCloud] The point cloud to reorient.
 * @param[mesh_filename] The mesh file to use to reorient the normals / compute the normals.
 * @param[rotation_mesh_filename] txt file giving a rotation that need to be applied to the mesh.
 * @param[real_normals] It takes the normals directly from the mesh.
 */
Mesh reorientNormals( PointCloudDiff<Scalar> & pointCloud, const std::string & mesh_filename, const double& scale, const std::vector<double>& rotations, const std::vector<double>& translations, bool real_normals ){

  // Open the mesh file
  Mesh mesh = loadMesh( mesh_filename );

  mesh.applyRotation( scale, rotations, translations );

  // [TODO] Change it using a KdTree.

  VectorType point_i;
  VectorType normal_i = VectorType::Zero();

  #pragma omp parallel for
  for ( int i = 0 ; i < pointCloud.points.size() ; i++ ){
      point_i = pointCloud.points[ i ];
      normal_i = mesh.closestNormal( point_i );

      if ( real_normals ){
          pointCloud.normals[i] = normal_i;
      }
      else {
          // Check if the normal is in the same direction as the point_i
          if ( pointCloud.normals[ i ].dot(normal_i) < 0 )
              pointCloud.normals[ i ] *= -1;
      }
  }  

  return mesh;
}

int main( int argc, char ** argv )
{
  CLI::App app{ "CAD processing" };

  std::string input_filename;
  app.add_option( "-i,--input", input_filename, "Input point cloud" )->required();

  std::string mesh_filename;
  app.add_option( "--mesh-input", mesh_filename, "Use mesh file to reorient the estimated normals." );

  std::vector <double> rotations = {0.0, 0.0, 0.0};
  app.add_option( "--rotations", rotations, "Rotation to apply to the mesh." );

  std::vector <double> translations = {0.0, 0.0, 0.0};
  app.add_option( "--translationX", translations[0], "Translation to apply on X axis to the mesh." );

  app.add_option( "--translationY", translations[1], "Translation to apply on Y axis to the mesh." );

  app.add_option( "--translationZ", translations[2], "Translation to apply on Z axis to the mesh." );

  double scale = 1.0;
  app.add_option( "--scale", scale, "Scale to apply to the mesh." );

  bool real_normals = false;
  app.add_flag( "--real-normals", real_normals, "Compute the normals directly from the mesh." );

  std::string output_filename;
  app.add_option( "-o,--output", output_filename, "Output point cloud" );

  kNN = 10;
  app.add_option( "-k,--kNN", kNN, "Number of neighbors for the normal estimation." );

  radius = 20.0;

  CLI11_PARSE( app, argc, argv );

  if ( ! std::filesystem::exists( input_filename ) )
  {
    spdlog::error( "File {} does not exist", input_filename );
    return EXIT_FAILURE;
  }

  if ( mesh_filename != "" && ! std::filesystem::exists( mesh_filename ) )
  {
    spdlog::error( "File {} does not exist", mesh_filename );
    return EXIT_FAILURE;
  }

  // To force the PCA to be only a PCA and not MLS (for the normal estimation)
  mls_iter = 1;


  PointCloudDiff<Scalar> original_pc;
  original_pc.loadPointCloud( input_filename );

  Mesh mesh;

  // PONCA KdTree
  spdlog::info( "Building the ponca KdTree..." );

  buildKdTree<PPAdapter, VectorType>(original_pc.points, original_pc.normals, ponca_kdtree);

  if ( ! real_normals ) {
    DifferentialQuantities<Scalar> diff_quantities = estimateDifferentialQuantities<ConstWeightFunc>("PCA");
    spdlog::info( "Mean number of neighbors {}", diff_quantities.statNeighbors().mean() );
    original_pc.normals = diff_quantities.normal();
  }
  
  if ( mesh_filename != "" )
    mesh = reorientNormals ( original_pc, mesh_filename, scale, rotations, translations, real_normals );

  if ( output_filename != "" )
  {
    original_pc.savePointCloud( output_filename );
  }

  return EXIT_SUCCESS;
}
