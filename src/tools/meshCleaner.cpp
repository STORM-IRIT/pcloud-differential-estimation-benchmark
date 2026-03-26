#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include <set>

#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/readPLY.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>

#include <polyscope/polyscope.h>
#include <polyscope/point_cloud.h>
#include "polyscope/messages.h"
#include <polyscope/surface_mesh.h>

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
    Eigen::MatrixXi getF() const { return F; }
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
    std::vector<std::array<int, 3>> getFVector() const {
      std::vector<std::array<int, 3>> F_vector(F.rows());
      for (int i = 0; i < F.rows(); ++i)
        F_vector[i] = {F(i, 0), F(i, 1), F(i, 2)};
      return F_vector;
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

    void saveMesh ( const std::string & filename ) const {
      igl::write_triangle_mesh(filename, V, F);
    }

    void clean_mesh ( const int & min_knn ) {
      std::vector<int> bad_idx;
      int one_good = 0;
      for (int i = 0; i < V.rows(); ++i){
        VectorType point = V.row(i).transpose();
        int nb_neighbors = 0;
        processNeighbors(point, 0, [&nb_neighbors]( int ){
                nb_neighbors ++;
            });
        if (nb_neighbors < min_knn)
          bad_idx.push_back(i);
        else
          one_good = i;
      }
      
      std::cout << "Number of bad points: " << bad_idx.size() << std::endl;
      std::cout << "Selected good point: " << one_good << std::endl;
      std::cout << "Its point coordinates: " << V.row(one_good) << std::endl;
      clean_mesh(bad_idx, one_good);

    }

    private: 

        bool is_bad_face ( std::vector<int> & bad_idx, int i ) {
          bool bad = ( std::find(bad_idx.begin(), bad_idx.end(), F(i, 0)) != bad_idx.end() ||
                 std::find(bad_idx.begin(), bad_idx.end(), F(i, 1)) != bad_idx.end() ||
                 std::find(bad_idx.begin(), bad_idx.end(), F(i, 2)) != bad_idx.end() );
          
          // compute max dist between the 3 points of the face
          VectorType p0 = V.row(F(i, 0)).transpose();
          VectorType p1 = V.row(F(i, 1)).transpose();
          VectorType p2 = V.row(F(i, 2)).transpose();

          return bad;
        }

        void clean_mesh ( std::vector<int> bad_idx, int one_good ) {
          SampleMatrixType new_V(V.rows(), 3);
          SampleMatrixType new_N(N.rows(), 3);
          std::vector< std::array<int, 3> > new_F;

          VectorType one_good_point = V.row(one_good).transpose();
          VectorType one_good_normal = N.row(one_good).transpose();


          for (int i = 0 ; i < V.rows(); ++i){
            if (std::find(bad_idx.begin(), bad_idx.end(), i) == bad_idx.end()){
              new_V.row(i) = V.row(i);
              new_N.row(i) = N.row(i);
            }
            else {
              new_V.row(i) = one_good_point;
              new_N.row(i) = one_good_normal;
            }
          }

          for (int i = 0; i < F.rows(); ++i){
            if (!is_bad_face(bad_idx, i)){
              new_F.push_back({F(i, 0), F(i, 1), F(i, 2)});
            }
          }

          Eigen::MatrixXi new_F_matrix(new_F.size(), 3);
          for (int i = 0; i < new_F.size(); ++i){
            new_F_matrix.row(i) << new_F[i][0], new_F[i][1], new_F[i][2];
          }

          V = new_V;
          F = new_F_matrix;
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
        std::cout << "Per face normals" << std::endl;
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

  // if ( per_face ){
  //   SampleMatrixType V_faces(F.rows(), 3);
  //   for (int i = 0; i < F.rows(); ++i){
  //       V_faces.row(i) = ( V.row( F( i, 0 ) ) + V.row( F( i, 1 ) ) + V.row( F( i, 2 ) ) ) / 3.0;
  //   }
  //   return Mesh(V_faces, N, F);
  // }

  return Mesh(V, N, F);
}

void translate_other_pc(PointCloudDiff<Scalar>& other_pc) {
    VectorType rotation = VectorType(0.0, -45.0, -45.0);  // Inverse of the original rotation
    VectorType translation = VectorType(0.0, 0.0, -10.0);  // Note the negative Z for inverse translation
    Scalar scale = Scalar(1.0 / 9.0);  // Inverse of the original scale

    rotation = rotation * M_PI / 180.0;

    Eigen::Matrix<Scalar, 3, 3> rot, rotY, rotZ;
    rotY << cos(rotation[1]), 0, sin(rotation[1]), 
            0, 1, 0, 
            -sin(rotation[1]), 0, cos(rotation[1]);
    rotZ << cos(rotation[2]), -sin(rotation[2]), 0, 
            sin(rotation[2]),  cos(rotation[2]), 0, 
            0, 0, 1;
    rot = rotY * rotZ ;


    #pragma omp for
    for (int i = 0; i < other_pc.points.size(); ++i) {
        other_pc.points[i] = scale * ( rot * ( other_pc.points[i] + translation ) );
        if (!other_pc.normals.empty()) {
            other_pc.normals[i] = rot * other_pc.normals[i];
            other_pc.normals[i].normalize();
        }
    }
}

void add_ground_truth ( const PointCloudDiff<Scalar> & original_pc, PointCloudDiff<Scalar>& other_pc ){

  kNN = 1;
  int size = other_pc.points.size();
  std::vector<Scalar> k1(size), k2(size), kmean(size), kgauss(size);
  std::vector<VectorType> normals(size), v1(size), v2(size);

  // Translate the other point cloud to align with the original point cloud
  translate_other_pc(other_pc);

  for (int i = 0; i < size; ++i){
    VectorType point = other_pc.points[i];
    int idx = -1;
    processNeighbors(point, 0, [&idx]( int i ){
            idx = i;
        });
    if (idx == -1){
      std::cerr << "Error: no neighbors found for point " << i << std::endl;
      continue;
    }
    k1[i] = original_pc.k1[idx];
    k2[i] = original_pc.k2[idx];
    kmean[i] = original_pc.mean[idx];
    kgauss[i] = original_pc.gauss[idx];
    normals[i] = original_pc.normals[idx];
    v1[i] = original_pc.v1[idx];
    v2[i] = original_pc.v2[idx];
  }

  other_pc.k1 = k1;
  other_pc.k2 = k2;
  other_pc.mean = kmean;
  other_pc.gauss = kgauss;
  other_pc.normals = normals;
  other_pc.v1 = v1;
  other_pc.v2 = v2;
}

void callback()
{
  // Create a window
  ImGui::PushItemWidth( 100 );
  // The size of the window, the position is set by default
  ImGui::SetNextWindowSize( ImVec2( 300, 600 ), ImGuiCond_FirstUseEver );
  // Save the camera settings
  if ( ImGui::Button( "Save camera settings" ) ){
      std::string base_name = "cameraSettings";
      std::string num = "0";
      std::string extension = ".json";
      std::string view = polyscope::view::getViewAsJson();
      // While it exists a file with the same name, add a number at the end of the name
      while ( std::filesystem::exists( base_name + num + extension ) ) num = std::to_string (std::stoi( num ) + 1 );
      std::ofstream file( base_name + num + extension );
      file << view;
      file.close();
  }
}

int main( int argc, char ** argv )
{
  CLI::App app{ "Mesh Cleaner" };

  std::string input_filename;
  app.add_option( "-i,--input", input_filename, "Input point cloud" )->required();

  std::string mesh_filename;
  app.add_option( "--mesh-input", mesh_filename, "Mesh to clean from input point cloud." );

  std::string xyz_filename;
  app.add_option( "--xyz-input", xyz_filename, "Use xyz file to add groundtruth from input point cloud." );

  std::string output_filename;
  app.add_option( "-o,--output", output_filename, "Output point cloud" );

  kNN = -1;
  app.add_option( "-k,--kNN", kNN, "Number of neighbors for the normal estimation." );

  int k = 2;
  app.add_option( "--kmin", k, "Number of neighbors needed at minimum to conserve a point of the mesh." );

  radius = 0.05;
  app.add_option( "-r,--radius", radius, "Radius for the research." );

  bool show_gui = false;
  app.add_flag( "--gui", show_gui, "Show the GUI." );

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

  if ( xyz_filename != "" && ! std::filesystem::exists( xyz_filename ) )
  {
    spdlog::error( "File {} does not exist", xyz_filename );
    return EXIT_FAILURE;
  }

  spdlog::info( "Loading the pointCloud..." );
  PointCloudDiff<Scalar> original_pc;
  original_pc.loadPointCloud( input_filename );
  spdlog::info( "Building the ponca KdTree..." );
  buildKdTree<PPAdapter, VectorType>(original_pc.points, original_pc.normals, ponca_kdtree);

  if (mesh_filename != ""){
    spdlog::info( "Loading the mesh..." );
    Mesh mesh;
    mesh = loadMesh( mesh_filename );
    spdlog::info( "Cleaning the mesh..." );
    mesh.clean_mesh(k);
    if (output_filename != ""){
      spdlog::info( "Saving the new mesh..." );
      mesh.saveMesh(output_filename);
    }
    if (show_gui){
      polyscope::init();
      polyscope::state::userCallback = callback;
      polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
      polyscope::view::lookAt( { 28, 18, 28 }, { 0, 0, 0 } );

      // Open the mesh file
      polyscope::registerSurfaceMesh("Mesh", mesh.getVVector(), mesh.getFVector());
      polyscope::show();
    }
  }

  if (xyz_filename != "") {
    spdlog::info ("Loading the xyz file...");
    PointCloudDiff<Scalar> xyz_pc;
    xyz_pc.loadPointCloud(xyz_filename);
    spdlog::info ("Adding the ground truth...");
    add_ground_truth(original_pc, xyz_pc);

    if (output_filename != ""){
      spdlog::info( "Saving the new pts file..." );
      xyz_pc.savePointCloud(output_filename);
    }

    if (show_gui){
      // Open the 2 point clouds
      polyscope::init();
      polyscope::state::userCallback = callback;
      polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
      polyscope::view::lookAt( { 28, 18, 28 }, { 0, 0, 0 } );
      
      polyscope::registerPointCloud("Original", original_pc.points);
      polyscope::getPointCloud("Original")->addVectorQuantity("normals", original_pc.normals);
      polyscope::registerPointCloud("Ground Truth", xyz_pc.points);
      polyscope::getPointCloud("Ground Truth")->addVectorQuantity("normals", xyz_pc.normals);
      polyscope::show();
    }

  }

  return EXIT_SUCCESS;
}
