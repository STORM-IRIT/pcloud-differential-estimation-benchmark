#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>
#include <set>
#include <polyscope/polyscope.h>
#include <polyscope/point_cloud.h>
#include "polyscope/messages.h"
#include <polyscope/surface_mesh.h>
#include "../IO/PointCloudDiff.h"
#include "../IO/PNGHandler.h"

#include <CLI11.hpp>

typedef Eigen::Vector<Scalar, 3>    Point;

float pointSize = 0.02f; // Rayon par défaut
std::vector<polyscope::PointCloud*> g_pointClouds; // Pointeurs vers tous les nuages
bool center = false;
double minBound = -1.0;
double maxBound = 1.0;
bool errorFile = false;

inline std::map<int, polyscope::view::UpDir> upDirections = {
  {0, polyscope::view::UpDir::XUp}, {1, polyscope::view::UpDir::YUp},
  {2, polyscope::view::UpDir::ZUp}, {3, polyscope::view::UpDir::NegXUp},
  {4, polyscope::view::UpDir::NegYUp}, {5, polyscope::view::UpDir::NegZUp}
};

std::map<polyscope::view::UpDir, int> upDirectionsToInt = {
  {polyscope::view::UpDir::XUp, 0}, {polyscope::view::UpDir::YUp, 1},
  {polyscope::view::UpDir::ZUp, 2}, {polyscope::view::UpDir::NegXUp, 3},
  {polyscope::view::UpDir::NegYUp, 4}, {polyscope::view::UpDir::NegZUp, 5}
};

void updatePointSizes() {
    for (auto& pc : g_pointClouds) {
        if (pc) pc->setPointRadius(pointSize, false);
    }
}

void loadCamera ( const std::string &filename )
{
  std::string view;
  std::ifstream file(filename);
  if (file.is_open()) {
    std::string line;
    std::getline(file, line); view = line;
    if (std::getline(file, line)) {
      pointSize = std::stof(line);
      updatePointSizes();
    }
    if (std::getline(file, line)) {
      polyscope::view::upDir = upDirections[std::stoi(line)];
    }
    file.close();
  }
  polyscope::view::setViewFromJson(view, false);
}

void saveCamera ( const std::string &camName )
{
  std::string view = polyscope::view::getViewAsJson();
  std::ofstream file(camName);
  file << view;
  file << "\n" << pointSize;
  file << "\n" << upDirectionsToInt[polyscope::view::getUpDir()];
  file.close();
}

void centering(std::vector<Point>& points, Point minBox, Point maxBox) {
    Point p_center = (minBox + maxBox) / 2.0;
    for (auto& p : points)
        p -= p_center;
}

void addColorQuantity ( polyscope::PointCloud * ps, std::string name, std::string propertyName, std::string & propertyToDisplay, std::vector<Point> & values)
{
  auto qt = ps->addColorQuantity( propertyName, values );
  if ( propertyName == propertyToDisplay )
    qt->setEnabled( true );
}

void addScalarQuantity ( polyscope::PointCloud * ps, std::string name, std::string propertyName, std::string & propertyToDisplay, std::vector<Scalar> & values)
{
  auto qt = ps->addScalarQuantity( propertyName, values );
  if ( propertyName == propertyToDisplay ){
    qt->setEnabled( true );
  }
  if ( ! errorFile )
    qt->setColorMap( "coolwarm" );
  else
    qt->setColorMap( "viridis" );
  qt->setMapRange( std::pair<Scalar, Scalar>( minBound, maxBound ) );
}

void addVectorQuantity ( polyscope::PointCloud * ps, std::string name, std::string propertyName, std::string & propertyToDisplay, std::vector<Point> & values)
{
  auto qt = ps->addVectorQuantity( propertyName, values );
  if ( propertyName == propertyToDisplay )
    qt->setEnabled( true );
}

std::pair<Point, Point> processBbox ( std::vector<Point> & points )
{
  Point min = { std::numeric_limits<Scalar>::max(), std::numeric_limits<Scalar>::max(), std::numeric_limits<Scalar>::max() };
  Point max = { std::numeric_limits<Scalar>::min(), std::numeric_limits<Scalar>::min(), std::numeric_limits<Scalar>::min() };

  for ( auto & p : points ){
    for ( auto i = 0u; i < 3; ++i ){
      if ( p[ i ] < min[ i ] )
        min[ i ] = p[ i ];
      if ( p[ i ] > max[ i ] )
        max[ i ] = p[ i ];
    }
  }

  spdlog::info( "Bounding box: min = ({}, {}, {}), max = ({}, {}, {})", min[ 0 ], min[ 1 ], min[ 2 ], max[ 0 ], max[ 1 ], max[ 2 ] );
  return std::make_pair( min, max );
}

void processErrorFile( std::string & filename, double pointRadius, std::string &property, bool useABS,  double minBound = -1.0, double maxBound = 1.0)
{
  // If the filename isn't a pts file, abort.
  if ( filename.find( ".pts" ) == std::string::npos ){
    std::cerr << "Error: " << filename << " is not a .pts file." << std::endl;
    return;
  }

  std::ifstream file( filename );
  if ( !file.is_open() ){
    std::cerr << "Error: could not open file " << filename << std::endl;
    return;
  }

  std::vector<Point> points;
  std::vector<Scalar> shapeIdxs, k1s, k2s, means, gausss, d1s, d2s, normals, poss;

  std::string x, y, z, shapeIdx, gauss, mean, k1, k2, d1, d2, normal, pos;

  std::string line;
  while ( std::getline( file, line ) ){
    if ( ( line.find( "#" ) ) != std::string::npos )
      spdlog::info( "Comment: {}", line );
    else {
      std::istringstream linestream( line );
      linestream >> x >> y >> z >> shapeIdx >> gauss >> mean >> k1 >> k2 >> d1 >> d2 >> normal >> pos;

      points.push_back( Point( std::stod( x ), std::stod( y ), std::stod( z ) ) );
      if (useABS){
        shapeIdxs.push_back( std::stod( shapeIdx ) );
        gausss.push_back( std::stod( gauss ) );
        means.push_back( std::stod( mean ) );
        k1s.push_back( std::stod( k1 ) );
        k2s.push_back( std::stod( k2 ) );
        d1s.push_back( std::stod( d1 ) );
        d2s.push_back( std::stod( d2 ) );
        normals.push_back( std::stod( normal ) );
        poss.push_back( std::stod( pos ) );
      }
      else {
        shapeIdxs.push_back( sqrt( std::stod( shapeIdx ) ) );
        gausss.push_back( sqrt( std::stod( gauss ) ) );
        means.push_back( sqrt( std::stod( mean ) ) );
        k1s.push_back( sqrt( std::stod( k1 ) ) );
        k2s.push_back( sqrt( std::stod( k2 ) ) );
        d1s.push_back( sqrt( std::stod( d1 ) ) );
        d2s.push_back( sqrt( std::stod( d2 ) ) );
        normals.push_back( sqrt( std::stod( normal ) ) );
        poss.push_back( sqrt( std::stod( pos ) ) );
      }
    }
  }
  spdlog::info( "Loaded {} points.", points.size() );
  file.close();

  std::pair<Point, Point> bbox = processBbox( points );
  if ( center )
    centering( points, bbox.first, bbox.second );

  auto ps = polyscope::registerPointCloud( filename, points );
  if ( pointRadius != 0.0 )
    ps->setPointRadius( pointRadius, false );

  addScalarQuantity( ps, "Shape Index error", "Shape Index", property, shapeIdxs );
  addScalarQuantity( ps, "Gauss error", "kGauss", property, gausss );
  addScalarQuantity( ps, "Mean error", "kMean", property, means );
  addScalarQuantity( ps, "kMin error", "kMin", property, k1s );
  addScalarQuantity( ps, "kMax error", "kMax", property, k2s );
  addScalarQuantity( ps, "kMin direction error", "kMin direction", property, d1s );
  addScalarQuantity( ps, "kMax direction error", "kMax direction", property, d2s );
  addScalarQuantity( ps, "Normal error", "Normal", property, normals );
  addScalarQuantity( ps, "Position error", "Position", property, poss );

}

std::pair<Point, Point> processFile( std::string & filename, double pointRadius, bool absolute_values, std::string &property, bool fast_render = false, bool vectorAsColor = false)
{
  PointCloudDiff<Scalar> pc;
  std::cout << "Process file called with " << filename << std::endl;
  pc.loadPointCloud( filename );
  std::cout << "Successfully loaded " << filename << std::endl;

  std::pair<Point, Point> bbox = processBbox ( pc.points );

  if ( center )
    centering( pc.points, bbox.first, bbox.second );

  std::cout << filename << " " << pc.points.size() << " points." << std::endl;

  if ( absolute_values ){
    std::cout << "Using absolute values." << std::endl;
    pc.switch_to_abs();
  }

  std::vector<Scalar> shapeIdx ( pc.points.size() );
  for ( auto i = 0u; i < pc.points.size(); ++i )
    shapeIdx[ i ] = (std::fabs( pc.k1 [ i ] - pc.k2[ i ] ) > 1e-6) ? ( 2.0 / M_PI ) * std::atan( ( pc.k1[ i ] + pc.k2[ i ] ) / ( pc.k1 [ i ] - pc.k2[ i ] ) ) : 0. ;

  auto ps = polyscope::registerPointCloud( filename, pc.points );

  if ( fast_render )
    ps->setPointRenderMode( polyscope::PointRenderMode::Quad );

  // Overriding point radius.
  if ( pointRadius != 0.0 )
    ps->setPointRadius( pointRadius, false );

  addScalarQuantity( ps, "kMin", "kMin", property, pc.k1 );
  addScalarQuantity( ps, "kMax", "kMax", property, pc.k2 );
  addScalarQuantity( ps, "kMean", "kMean", property, pc.mean );
  addScalarQuantity( ps, "kGauss", "kGauss", property, pc.gauss );
  addScalarQuantity( ps, "Shape Index", "Shape Index", property, shapeIdx );
  if ( !vectorAsColor )
  {
    addVectorQuantity( ps, "Normal", "Normal", property, pc.normals );
    addVectorQuantity( ps, "kMin Direction", "kMin Direction", property, pc.v1 );
    addVectorQuantity( ps, "kMax Direction", "kMax Direction", property, pc.v2 );
  } else
  {
    auto normals = pc.normals;
    auto kMinDirs = pc.v1;
    auto kMaxDirs = pc.v2;
    for ( auto i = 0u; i < pc.points.size(); ++i )
    {
      normals[ i ] = 0.5 * ( pc.normals[ i ] + Point( 1., 1., 1. ) );
      kMinDirs[ i ] = 0.5 * ( pc.v1[ i ] + Point( 1., 1., 1. ) );
      kMaxDirs[ i ] = 0.5 * ( pc.v2[ i ] + Point( 1., 1., 1. ) );
    }
    addColorQuantity( ps, "Normal", "Normal", property, normals);
    addColorQuantity( ps, "kMin Direction", "kMin Direction", property, kMinDirs);
    addColorQuantity( ps, "kMax Direction", "kMax Direction", property, kMaxDirs);
  }

  return bbox;
}

void createColorBarAsMesh () {
  // Create a mesh
  std::vector<Point> points = { { -1.5, 0, 0 }, {-1.5, 0.5, 0 }, {1.5, 0.5, 0 }, { 1.5, 0, 0 } };
  std::vector<std::vector<size_t>> faces = { { 0, 1, 2 }, { 0, 2, 3 } };
  auto mesh = polyscope::registerSurfaceMesh( "ColorBar", points, faces );
  mesh->setBackFacePolicy( polyscope::BackFacePolicy::Identical );
  mesh->setMaterial( "flat" );
  // min and max vertex scalar quantities
  double minValue = minBound + 0.02 * ( maxBound - minBound );
  double maxValue = maxBound - 0.02 * ( maxBound - minBound );
  auto qt = mesh->addVertexScalarQuantity( "ColorBar", std::vector<double>{ minValue, minValue, maxValue, maxValue }, polyscope::DataType::SYMMETRIC );
  qt->setEnabled( true );
  if ( ! errorFile )
    qt->setColorMap( "coolwarm" );
  else
    qt->setColorMap( "viridis" );

  qt->setMapRange( std::pair<double, double>( minBound, maxBound ) );
}

// ----------------- Callback GUI -----------------
void callback() {
    ImGui::PushItemWidth(100);
    ImGui::SetNextWindowSize(ImVec2(300, 600), ImGuiCond_FirstUseEver);

    // Slider global pour le rayon des points
    if (ImGui::SliderFloat("Point Size", &pointSize, 0.001f, 0.1f)) {
        updatePointSizes();
    }

    // Bouton pour sauvegarder caméra + pointSize + upDir
    if (ImGui::Button("Save camera settings")) {
        std::string base_name = "cameraSettings";
        std::string num = "0";
        std::string extension = ".txt";

        while (std::filesystem::exists(base_name + num + extension))
            num = std::to_string(std::stoi(num) + 1);

        std::ofstream file(base_name + num + extension);
        if (!file.is_open()) {
            spdlog::error("Could not open file to save {}", base_name + num + extension);
            return;
        }
        saveCamera( base_name + num + extension );
        spdlog::info("Saved camera settings to {}", base_name + num + extension);
    }
}

// ----------------- Main -----------------
int main(int argc, char** argv) {
    CLI::App app{ "displayPTS" };
    std::vector<std::string> filenames;
    app.add_option("-i,--input,1", filenames, "Input point clouds");

    bool colorbar = false;
    app.add_flag("--colorbar", colorbar, "Display a colorbar.");

    std::string camFilename = "";
    app.add_option("-c, --camera", camFilename, "Camera filename.");

    std::string screenshotFilename = "";
    app.add_option("--screenshot", screenshotFilename, "Screenshot filename.");

    double pointRadius = 0.0;
    app.add_option("--pointRadius", pointRadius, "Override the default point radius.");

    bool absolute_values = false;
    app.add_flag("--absolute_values", absolute_values, "Use absolute values.");

    errorFile = false;
    app.add_flag("--errorFile", errorFile, "Load a .pts file containing the errors.");

    std::string property = "kMean";
    app.add_option("--property", property, "Property to display.");

    minBound = -1.0;
    app.add_option("--minBound", minBound, "Minimum bound for the colormap.");

    maxBound = 1.0;
    app.add_option("--maxBound", maxBound, "Maximum bound for the colormap.");

    bool fast_render = false;
    app.add_flag("--fast", fast_render, "Fast render mode.");

    bool useABS = false;
    app.add_flag("--abs", useABS, "Use absolute errors, default RMSE for each points.");

    center = false;
    app.add_flag("--center", center, "Center the point cloud.");

    float shadow = 0.5;
    app.add_option("--shadowDark", shadow, "Drop shadow darkness in [0, 1]. Default 0.5.");

    int blur = 10;
    app.add_option("--shadowBlur", blur, "Drop shadow blurring iterations. Default 10.");

    bool diagCrop = false;
    app.add_flag("--diagCrop", diagCrop, "Crop screenshot diagonally. (Top diag for estimation, bottom diag for error)");

    bool convertToAbs = false;
    app.add_flag("--convertToAbs", convertToAbs, "Convert input point cloud to absolute values.");

    bool vectorAsColor = false;
    app.add_flag("--vectorAsColor", vectorAsColor, "Display vector quantities as colors");

    CLI11_PARSE(app, argc, argv);

    if ( convertToAbs )
    {
      PointCloudDiff<Scalar> pc;
        for ( auto & filename : filenames )
        {
            if ( filename.find( ".txt" ) != std::string::npos || std::filesystem::is_directory( filename ) )
                continue;
            if ( filename.find( ".pts" ) == std::string::npos )
            {
                std::cerr << "Error: " << filename << " is not a .pts file." << std::endl;
                continue;
            }
            // New filename
            std::string newFilename = filename.substr( 0, filename.find_last_of( "." ) ) + "_abs.pts";
            pc.loadPointCloud( filename );
            pc.switch_to_abs();
            pc.savePointCloud( newFilename );
            std::cout << "Saved absolute point cloud to " << newFilename << std::endl;
        }
      exit(0);
    }

    polyscope::options::allowHeadlessBackends = true;
    polyscope::options::autocenterStructures = false;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::ShadowOnly;
    polyscope::options::shadowDarkness = shadow;
    polyscope::options::shadowBlurIters = blur;
    polyscope::options::verbosity = 0;
    polyscope::init();

    if (!camFilename.empty()) {
        loadCamera(camFilename);
    }

    Point minBox = { std::numeric_limits<Scalar>::max(), std::numeric_limits<Scalar>::max(), std::numeric_limits<Scalar>::max() };
    Point maxBox = { std::numeric_limits<Scalar>::min(), std::numeric_limits<Scalar>::min(), std::numeric_limits<Scalar>::min() };

    if (colorbar) createColorBarAsMesh();

    // Process all files
    for (auto& filename : filenames) {
        if (filename.find(".txt") != std::string::npos || std::filesystem::is_directory(filename))
            continue;

        std::pair<Point, Point> bbox;
        if (errorFile)
            processErrorFile(filename, pointSize, property, useABS);
        else
            bbox = processFile(filename, pointSize, absolute_values, property, fast_render, vectorAsColor);

        // Ajouter le nuage de points à notre vecteur global
        auto psCloud = polyscope::getPointCloud(filename);
        if (psCloud) {
            g_pointClouds.push_back(psCloud);
            psCloud->setPointRadius(pointSize, false); // appliquer pointSize global
        }

        for (auto i = 0u; i < 3; ++i) {
            if (bbox.first[i] < minBox[i]) minBox[i] = bbox.first[i];
            if (bbox.second[i] > maxBox[i]) maxBox[i] = bbox.second[i];
        }
    }

    spdlog::info("Total BBOX : min = ({}, {}, {}), max = ({}, {}, {})",
        minBox[0], minBox[1], minBox[2],
        maxBox[0], maxBox[1], maxBox[2]);

    if (!screenshotFilename.empty()) {
        polyscope::options::ssaaFactor = 4;
        polyscope::screenshot(screenshotFilename, true);
        PNGHandler::Image img;
        if(!PNGHandler::readPNG(screenshotFilename, img)) {
          printf("image reading error\n");
        }
        const auto box = PNGHandler::findBoundingBox(img);
        PNGHandler::Image crop;
        if(!PNGHandler::cropToBox(img, crop, box)) {
          std::cout << "Error cropping image " << screenshotFilename << std::endl;
        }
        if ( diagCrop ) {
          PNGHandler::Image cropDiag;
          if ( errorFile ) {
            if(!PNGHandler::cropDiagonal(crop, false, cropDiag)) {
              std::cout << "Error cropping image " << screenshotFilename << std::endl;
            }
          } else {
            if(!PNGHandler::cropDiagonal(crop, true, cropDiag)) {
              std::cout << "Error cropping image " << screenshotFilename << std::endl;
            }
          }
          crop = cropDiag;
        }
        if(!PNGHandler::writePNG(screenshotFilename, crop)) {
          std::cout << "Error writing image " << screenshotFilename << std::endl;
        }
    } else {
        polyscope::state::userCallback = callback;
        polyscope::show();
    }

    return EXIT_SUCCESS;
}