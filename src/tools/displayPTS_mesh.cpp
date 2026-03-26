#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include <set>

#include "../IO/PointCloudDiff.h"

#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/readPLY.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>

#include <polyscope/polyscope.h>
#include <polyscope/point_cloud.h>
#include "polyscope/messages.h"
#include <polyscope/surface_mesh.h>


#include <CLI11.hpp>

typedef Eigen::Vector<Scalar, 3>    Point;
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

void processMesh ( const std::string & filename, bool per_face=true ) {
  std::cout << "Process file called with " << filename << std::endl;
  Mesh mesh = loadMesh(filename, per_face);
  std::cout << "Successfully loaded " << filename << std::endl;

  auto ps = polyscope::registerSurfaceMesh( filename, mesh.getVVector(), mesh.getFVector() );

}

void translate_other_pc(PointCloudDiff<Scalar>& pc) {
  // // Implicit
  // VectorType rotation = VectorType(0.0, -45.0, -45.0);  // Inverse of the original rotation
  // VectorType translation = VectorType(0.0, 0.0, -10.0);  // Note the negative Z for inverse translation
  // Scalar scale = Scalar(1.0 / 9.0);  // Inverse of the original scale

  VectorType rotation = VectorType(0.0, -45.0, -45.0);  // Inverse of the original rotation
  VectorType translation = VectorType(-20.0, 10.0, 15.0);  // Note the negative Z for inverse translation
  Scalar scale = Scalar(2.0);  // Inverse of the original scale

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
  for (int i = 0; i < pc.points.size(); ++i) {
    pc.points[i] = scale * ( rot * ( pc.points[i] + translation ) );
      if (!pc.normals.empty()) {
        pc.normals[i] = rot * pc.normals[i];
        pc.normals[i].normalize();
      }
  }
}

void processFile( std::string & filename, double pointRadius, bool fast_render = false)
{
  PointCloudDiff<Scalar> pc;
  std::cout << "Process file called with " << filename << std::endl;
  pc.loadPointCloud( filename );

  translate_other_pc(pc);

  std::cout << "Successfully loaded " << filename << std::endl;

  auto ps = polyscope::registerPointCloud( filename, pc.points );

  if ( fast_render )
    ps->setPointRenderMode( polyscope::PointRenderMode::Quad );

  // Overriding point radius.
  if ( pointRadius != 0.0 )
    ps->setPointRadius( pointRadius );

  ps->addVectorQuantity( "Normal", pc.normals );

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
  CLI::App app{ "displayPTS" };
  std::vector<std::string> filenames;
  app.add_option( "-i,--input,1", filenames, "Input point clouds" );

  std::string camFilename = "";
  app.add_option( "-c, --camera", camFilename, "Camera filename)." );

  std::string screenshotFilename = "";
  app.add_option( "--screenshot", screenshotFilename, "Screenshot filename)." );

  double pointRadius = 0.0;
  app.add_option( "--pointRadius", pointRadius, "Override the default point radius." );

  bool fast_render = false;
  app.add_flag( "--fast", fast_render, "Fast render mode." );

  float shadow = 0.5;
  app.add_option( "--shadowDark", shadow, "Drop shadow darkness in [0, 1]. Default 0.5." );

  int blur = 10;
  app.add_option( "--shadowBlur", shadow, "Drop shadow blurring iterations. Default 10." );

  std::string upDir = "Z";
  app.add_option( "--upDir", upDir, "Up direction for the camera (Z, negZ, Y, negY, X, negX)." );

  bool per_face_normal=true;
  app.add_flag( "--faceNorm", per_face_normal, "Take the normals from the faces of the obj." );

  CLI11_PARSE( app, argc, argv );

  // Initialize polyscope
  polyscope::init();

  // Set the callback
  polyscope::state::userCallback = callback;

  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
  polyscope::options::shadowDarkness = shadow;
  polyscope::options::shadowBlurIters = blur;

  if ( upDir == "Z" )
    polyscope::view::setUpDir(polyscope::UpDir::ZUp);
  else if ( upDir == "negZ" )
    polyscope::view::setUpDir(polyscope::UpDir::NegZUp);
  else if ( upDir == "Y" )
    polyscope::view::setUpDir(polyscope::UpDir::YUp);
  else if ( upDir == "negY" )
    polyscope::view::setUpDir(polyscope::UpDir::NegYUp);
  else if ( upDir == "X" )
    polyscope::view::setUpDir(polyscope::UpDir::XUp);
  else if ( upDir == "negX" )
    polyscope::view::setUpDir(polyscope::UpDir::NegXUp);
  else {
    spdlog::error( "Error: unknown up direction {}", upDir );
  }

  // polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
  polyscope::view::lookAt( { 28, 18, 28 }, { 0, 0, 0 } );
  if ( camFilename != "" ){
    std::ifstream camFile( camFilename );
    if ( !camFile.is_open() ){
      spdlog::error( "Error: could not open file {}", camFilename );
      return EXIT_FAILURE;
    }
    std::string cam_json( ( std::istreambuf_iterator<char>( camFile ) ), std::istreambuf_iterator<char>() );
    spdlog::info( "Loaded camera :\n\t {}", cam_json );
    if ( cam_json.empty() )
      spdlog::info("The camera file is empty");
    else 
      polyscope::view::setViewFromJson(cam_json, false);

  }

  // Process all files
  for ( auto & filename : filenames ){
    // if the file ends with .txt or is a dir, continue 
    if ( filename.find( ".txt" ) != std::string::npos || std::filesystem::is_directory( filename ) )
      continue;
    if ( filename.find( ".obj" ) != std::string::npos ) {
      processMesh( filename, per_face_normal );
      continue;
    }
    processFile( filename, pointRadius, fast_render );

  }

  if ( screenshotFilename != "" ){
    polyscope::options::ssaaFactor      = 4;
    polyscope::screenshot( screenshotFilename, true );
  }
  else {
    // Render as quad if fast_render is set
    polyscope::show();
  }

  return EXIT_SUCCESS;
}
