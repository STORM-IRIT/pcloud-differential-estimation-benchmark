#include <iostream>
#include <sstream>
#include <random>
#include <vector>
#include <array>
#include <fstream>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/helpers/Shortcuts.h>
#include <DGtal/helpers/ShortcutsGeometry.h>
#include <Eigen/Dense>

// includes for std::clamp
#include <algorithm>
#include <cmath>


#include <SamplerLloyd.hpp>
#include <CLI11.hpp>

// DGtal
using namespace DGtal;
using namespace Z3i;
typedef Shortcuts<Z3i::KSpace> SH3;
typedef ShortcutsGeometry<Z3i::KSpace> SHG3;

// The implicit shape
CountedPtr<SH3::ImplicitShape3D> implicitShape;

// Main naive pcloud container
std::vector<SH3::RealPoint> pc;
// Quantities
std::vector<double> exactGauss;
std::vector<double> exactMean;
std::vector<RealPoint> projPoints;
std::vector<double> exactK1;
std::vector<double> exactK2;
std::vector<RealPoint> exactDir1;
std::vector<RealPoint> exactDir2;
std::vector<RealPoint> exactN;

// generate random seed
std::random_device random_seed;    
std::default_random_engine generator;
unsigned int seed = 0;


std::string formatDouble(double value, int precision = 6) {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(precision) << value;
    std::string s = stream.str();

    // Erase trailing zeros and the decimal point if necessary
    s.erase(s.find_last_not_of('0') + 1, std::string::npos);
    if (s.back() == '.') {
        s.pop_back();
    }
    
    // change the "." to ","
    std::replace( s.begin(), s.end(), '.', ',');
    return s;
}

std::vector<RealPoint> rescalePoints (){
  std::vector<RealPoint> rescaledPoints;
  Eigen::Vector3d baryCenter({0.0, 0.0, 0.0});
  // Compute barycenter
  for (const auto &p : pc){
    baryCenter[0] += p[0];
    baryCenter[1] += p[1];
    baryCenter[2] += p[2];
  }
  baryCenter /= double(pc.size());

  // Compute max distance
  Eigen::Vector3d maxDist({0.0, 0.0, 0.0});
  for (const auto &p : pc){
    Eigen::Vector3d current = Eigen::Vector3d(p[0], p[1], p[2]) - baryCenter;
    maxDist = maxDist.cwiseMax(current.cwiseAbs());
  }
  double maxDistNorm = maxDist.maxCoeff();

  // Rescale
  rescaledPoints.reserve(pc.size());
  for (int i = 0 ; i < pc.size() ; i++){
    RealPoint p = pc[i];
    // rescaledPoints.push_back((p - baryCenter) / maxDist);
    rescaledPoints.push_back((p - RealPoint(baryCenter[0], baryCenter[1], baryCenter[2])) / maxDistNorm);
    exactGauss [i] = exactGauss [i] * (maxDistNorm * maxDistNorm);
    exactMean [i] = exactMean [i] * maxDistNorm;
    exactK1 [i] = exactK1 [i] * maxDistNorm;
    exactK2 [i] = exactK2 [i] * maxDistNorm;
  }

  
  return rescaledPoints;
}

std::pair<std::vector<SH3::RealPoint>, std::vector<RealPoint>> addNoise( double noise_position, double noise_normal, double noise_flip, unsigned int& current_seed )
{

  std::vector <SH3::RealPoint> noisy_pc;
  std::vector<RealPoint> noisy_normals;

  current_seed = (seed == 0)? random_seed() : seed;
  
  generator.seed( current_seed );

  std::cout << " =======================>  seed = " << current_seed << std::endl;

  if (noise_position != 0){
    std::normal_distribution<double> dis_pos( 0, noise_position );
    
    for ( auto & v : pc )
    {
      SH3::RealPoint v_noisy = v;
      v_noisy[ 0 ] += dis_pos( generator );
      v_noisy[ 1 ] += dis_pos( generator );
      v_noisy[ 2 ] += dis_pos( generator );
      noisy_pc.push_back(v_noisy);
    }
  }
  else {
    noisy_pc = pc;
  }

  generator.seed(current_seed);

  if (noise_normal != 0){
    std::normal_distribution<double> dis_norm( 0, noise_normal );
    double max_degree = 0;
    for ( auto & n : exactN )
    {
      RealPoint n_noisy = n;
      n_noisy[ 0 ] += dis_norm( generator );
      n_noisy[ 1 ] += dis_norm( generator );
      n_noisy[ 2 ] += dis_norm( generator );
      double norm = sqrt(n_noisy[0]*n_noisy[0] + n_noisy[1]*n_noisy[1] + n_noisy[2]*n_noisy[2]);
      n_noisy[0] /= norm;
      n_noisy[1] /= norm;
      n_noisy[2] /= norm;

      double dot_prod = n_noisy[0] * n[0] + n_noisy[1] * n[1] + n_noisy[2] * n[2];
      dot_prod = std::clamp(dot_prod, -1.0, 1.0);

      double degree = acos(dot_prod) * 180 / M_PI;
      if (degree > max_degree){
        max_degree = degree;
      }

      noisy_normals.push_back(n_noisy);
    }
    std::cout << "For standard deviation on normal = " << noise_normal << " the max degree = " << max_degree << std::endl;
  }
  else {
    noisy_normals = exactN;
  }

  generator.seed(current_seed);

  if (noise_flip != 0){
    std::uniform_real_distribution<double> dis_flip( 0, 1 );
    for ( auto & n : noisy_normals )
    {
      if (dis_flip( generator ) < noise_flip){
        n[ 0 ] *= -1;
        n[ 1 ] *= -1;
        n[ 2 ] *= -1;
      }
    }
  }

  return std::make_pair(noisy_pc, noisy_normals);

}

std::vector<SH3::RealPoint> addOutliers( double ratio, double noise_position, unsigned int& current_seed )
{
  std::vector <SH3::RealPoint> noisy_pc;
  
  current_seed = (seed == 0) ? random_seed() : seed;
  
  generator.seed( current_seed );

  std::cout << " =======================>  seed = " << current_seed << std::endl;

  std::uniform_real_distribution<double> dis( 0, 1 );
  std::normal_distribution<double> dis_noise( 0, noise_position );
  for ( auto & p : pc )
  {
    SH3::RealPoint p_noisy = p;
    if (dis(generator) < ratio){
      p_noisy[ 0 ] += dis_noise( generator );
      p_noisy[ 1 ] += dis_noise( generator );
      p_noisy[ 2 ] += dis_noise( generator );
    }
    noisy_pc.push_back(p_noisy);
  }
  return noisy_pc;
}


// hard-coded bounds used for random sampling
void get_bounds( const std::string & surface, double & min, double & max )
{
  min = -15.5;
  max = +15.5;
  if ( surface == "sphere1" )
  {
    min = -2;
    max = +2;
    return;
  }
  if ( surface == "sphere9" )
  {
    min = -10;
    max = +10;
    return;
  }
  if ( surface == "ellipsoid" )
  {
    min = -10;
    max = +10;
    return;
  }
  if ( surface == "cylinder" )
  {
    min = -15.5;
    max = +15.5;
    return;
  }
  if ( surface == "torus" )
  {
    min = -10;
    max = +10;
    return;
  }
  if ( surface == "rcube" )
  {
    min = -15.5;
    max = +15.5;
    return;
  }
  if ( surface == "goursat" )
  {
    min = -15.5;
    max = +15.5;
    return;
  }
  if ( surface == "goursat-hole" )
  {
    min = -15.5;
    max = +15.5;
    return;
  }
  if ( surface == "distel" )
  {
    min = -15.5;
    max = +15.5;
    return;
  }
  if ( surface == "leopold" )
  {
    min = -10;
    max = +10;
    return;
  }
  if ( surface == "diabolo" )
  {
    min = -15.5;
    max = +15.5;
    return;
  }
  if ( surface == "heart" )
  {
    min = -2;
    max = +2;
    return;
  }
  if ( surface == "crixxi" )
  {
    min = -2;
    max = +2;
    return;
  }
  std::cout << "Warning: unkown surface type, bounds set to dedault (-15.5,+15.5)" << std::endl;
}

// ImGUI handler
void samplingAndGenerate( const std::string & surface, unsigned int nbPts, unsigned int nbSteps, double epsilonRejection )
{
  auto params = SH3::defaultParameters() | SHG3::defaultParameters() | SHG3::parametersGeometryEstimation();
  params( "polynomial", surface );
  implicitShape = SH3::makeImplicitShape3D( params );
  double mind, maxd;
  get_bounds( surface, mind, maxd );
  std::cout << "Extent = " << mind << " " << maxd << std::endl;

  // std::default_random_engine generator_base;
  generator.seed(1);
  std::uniform_real_distribution<double> distribution( mind, maxd );

  pc.clear();

  std::cout << "Generating samples by rejection..." << std::endl;
  // Rejection sampling close to the boundary
  unsigned int cpt = 0;
  while ( cpt < nbPts )
  {
    RealPoint p( { distribution( generator ), distribution( generator ), distribution( generator ) } );
    auto q = implicitShape->nearestPoint( p );
    if ( ( p - q ).norm() < epsilonRejection )
    {
      pc.push_back( p );
      ++cpt;
    }
  }
  std::cout << nbPts << " samples." << std::endl;
  // Optimization
  std::cout << "Optimizing..." << std::endl;
  SamplerLloyd optimizer;
  cpt = 0u;
  std::vector<RealPoint> bary( nbPts );
  while ( cpt < nbSteps )
  {
    std::cout << "+" << std::flush;
    optimizer.init( pc, mind, maxd );
    optimizer.Lloyd_step();
    optimizer.getSamples( bary );
    for ( auto i = 0u; i < pc.size(); ++i )
    {
      auto q = implicitShape->nearestPoint( pc[ i ] );
      pc[ i ] += 0.8 * ( q - pc[ i ] ) + 0.2 * ( bary[ i ] - pc[ i ] ); // weighting the gradients
    }
    ++cpt;
  }
  // Last step is a projection one.
  for ( auto i = 0u; i < pc.size(); ++i )
  {
    auto q  = implicitShape->nearestPoint( pc[ i ] );
    pc[ i ] = q;
  }
  std::cout << std::endl;
  std::cout << "done " << pc.size() << std::endl;

  // Reproject
  projPoints.clear();
  exactMean.clear();
  exactGauss.clear();
  exactK1.clear();
  exactK2.clear();
  exactN.clear();
  exactDir1.clear();
  exactDir2.clear();

  for ( auto p : pc )
  {
    exactMean.push_back( implicitShape->meanCurvature( p ) );
    exactGauss.push_back( implicitShape->gaussianCurvature( p ) );
    double k1, k2;
    implicitShape->principalCurvatures( p, k1, k2 );
    exactK1.push_back( k1 );
    exactK2.push_back( k2 );
    RealPoint d1, d2;
    implicitShape->principalDirections( p, d1, d2 );
    exactDir1.push_back( d1 );
    exactDir2.push_back( d2 );
    exactN.push_back( implicitShape->gradient( p ).getNormalized() );
  }
}

int main( int argc, char ** argv )
{
  CLI::App app{ "Sampling implicit polynomial surface" };

  unsigned int nbPts = 10;
  app.add_option( "-n,--nbPts", nbPts, "Number of samples" );
  double epsilon = .5;
  app.add_option( "-e,--epsilon", epsilon, "Epsilon for the rejection sampling" );
  double nbSteps = 40;
  app.add_option( "--nbSteps", nbSteps, "Number of steps of the Lloyd relaxation" );
  bool XYZOnly = false;
  app.add_flag( "--XYZOnly", XYZOnly, "Export point coordinates only (default: false)" );
  std::string outputFilename;
  app.add_option( "-o,--output", outputFilename, "Output filename" )->required();
  std::string surface = "goursat";
  app.add_option( "-s,--surface", surface, "Surface [goursat|...]" );
  bool rescale = false;
  app.add_flag("-r, --rescale", rescale, "Rescale the surface to [-1, 1] (default: false)");
  bool pos_norm = false;
  app.add_flag("-p, --pos_norm", pos_norm, "Generate the point cloud, containing only pos and normal informations (default: false)");
  std::vector<double> noise_normal;
  app.add_option("-N, --noise_normal", noise_normal, "Add noise to the normals as a list of noise (default: None)");
  std::vector<double> noise_pos;
  app.add_option("-P, --noise_pos", noise_pos, "Add noise to the positions as a list of noise (default: 0.0)");
  std::vector<double> noise_flip;
  app.add_option("-F, --noise_flip", noise_flip, "Add flip to the normals as a list of flip ratio (default: 0.0)");
  std::vector<double> outlier_rate;
  app.add_option("--outlier_rate", outlier_rate, "Add outliers to the point cloud as a list of outlier rate (default: 0.0)");
  std::vector<double> outlier_noise;
  app.add_option("--outlier-noise", outlier_noise, "Add noise to the outliers as a list of noise (default: None)");
  int nb_gen_noise = 0;
  app.add_option("--nb_gen_noise", nb_gen_noise, "Must be set to > 0 to apply noise nb_gen_noise times (default: 0)");
  seed = 0;
  app.add_option("--seed", seed, "Seed for the random generator (default: 0 (random seed) NO NEGATIVE VALUE" );

  CLI11_PARSE( app, argc, argv );

  noise_pos.push_back (0);
  noise_normal.push_back (0);
  noise_flip.push_back (0);
  outlier_rate.push_back (0);

  // Main loop
  samplingAndGenerate( surface, nbPts, nbSteps, epsilon );

  std::string file_extension = outputFilename.substr(outputFilename.find_last_of(".") + 1);

  if (rescale){
    pc = rescalePoints ();
  }

    std::ofstream ofs( outputFilename, std::ofstream::out );
    ofs << "# x y z Gauss_Curvature Mean_Curvature nx ny nz k1 k2 d1x d1y d1z "
          "d2x d2y d2z\n";
    for ( auto i = 0u; i < nbPts; ++i )
    {
      ofs << pc[ i ][ 0 ] << " " << pc[ i ][ 1 ] << " " << pc[ i ][ 2 ];
      if ( !XYZOnly )
      {
        ofs << " " << exactGauss[ i ] << " " << exactMean[ i ] << " ";
        ofs << exactN[ i ][ 0 ] << " " << exactN[ i ][ 1 ] << " " << exactN[ i ][ 2 ] << " ";
        ofs << exactK1[ i ] << " " << exactK2[ i ] << " ";
        ofs << exactDir1[ i ][ 0 ] << " " << exactDir1[ i ][ 1 ] << " " << exactDir1[ i ][ 2 ] << " ";
        ofs << exactDir2[ i ][ 0 ] << " " << exactDir2[ i ][ 1 ] << " " << exactDir2[ i ][ 2 ];
      }
      ofs << std::endl;
    }
    ofs.close();
  
  if (pos_norm){
    std::string name = outputFilename.substr(0, outputFilename.size()-4) + "_pos_norm.pts";
    std::ofstream ofs( name, std::ofstream::out );
    ofs << "# x y z nx ny nz\n";
    for ( auto i = 0u; i < nbPts; ++i )
    {
      ofs << pc[ i ][ 0 ] << " " << pc[ i ][ 1 ] << " " << pc[ i ][ 2 ];
        ofs << " " << exactN[ i ][ 0 ] << " " << exactN[ i ][ 1 ] << " " << exactN[ i ][ 2 ];
      ofs << std::endl;
    }
    ofs.close();
  }

  for (int i = 0; i < nb_gen_noise; i++){
    for (int idx_noise_pos = 0; idx_noise_pos < noise_pos.size() ; idx_noise_pos++)
      for (int idx_noise_normale = 0; idx_noise_normale < noise_normal.size() ; idx_noise_normale++)
        for (int idx_noise_flip = 0; idx_noise_flip < noise_flip.size() ; idx_noise_flip++)
          {
            if (noise_pos[idx_noise_pos] == 0 && noise_normal[idx_noise_normale] == 0 && noise_flip[idx_noise_flip] == 0){
              continue;
            }
            unsigned int current_seed = 0;
            std::pair <std::vector<SH3::RealPoint>, std::vector<RealPoint>> noisy_pc_norm = addNoise(noise_pos[idx_noise_pos], noise_normal[idx_noise_normale], noise_flip[idx_noise_flip], current_seed);
            std::vector<SH3::RealPoint> noisy_pc = noisy_pc_norm.first;
            std::vector<RealPoint> noisy_normals = noisy_pc_norm.second;

            std::string noise_name = formatDouble(noise_pos[idx_noise_pos]) + "_" + formatDouble(noise_normal[idx_noise_normale]) + "_" + formatDouble(noise_flip[idx_noise_flip]);
            // replace the "." in the string by "," (to avoid problems with the extension)
            std::replace( noise_name.begin(), noise_name.end(), '.', ',');

            std::string name = outputFilename.substr(0, outputFilename.size()-4) + "_" + noise_name + "." + file_extension;
            std::ofstream ofs( name, std::ofstream::out );
            
            if ( !XYZOnly )
              ofs << "# x y z Gauss_Curvature Mean_Curvature nx ny nz k1 k2 d1x d1y d1z "
                    "d2x d2y d2z\n";
            else 
              ofs << "# x y z\n";
            for ( auto i = 0u; i < nbPts; ++i )
            {
              ofs << noisy_pc[ i ][ 0 ] << " " << noisy_pc[ i ][ 1 ] << " " << noisy_pc[ i ][ 2 ];
                if ( !XYZOnly )
                {
                  ofs << " " << exactGauss[ i ] << " " << exactMean[ i ] << " ";
                  ofs << noisy_normals[ i ][ 0 ] << " " << noisy_normals[ i ][ 1 ] << " " << noisy_normals[ i ][ 2 ] << " ";
                  ofs << exactK1[ i ] << " " << exactK2[ i ] << " ";
                  ofs << exactDir1[ i ][ 0 ] << " " << exactDir1[ i ][ 1 ] << " " << exactDir1[ i ][ 2 ] << " ";
                  ofs << exactDir2[ i ][ 0 ] << " " << exactDir2[ i ][ 1 ] << " " << exactDir2[ i ][ 2 ];
                }
              ofs << std::endl;
            }
            ofs.close();
          }
    for (int idx_noise_outlier = 0; idx_noise_outlier < outlier_rate.size(); idx_noise_outlier++){
      for (int idx_strength_outlier = 0; idx_strength_outlier < outlier_noise.size(); idx_strength_outlier++){
        if (outlier_rate[idx_noise_outlier] == 0 || outlier_noise[idx_strength_outlier] == 0){
          continue;
        }
          
        unsigned int current_seed = 0;
        std::vector<SH3::RealPoint> noisy_pc = addOutliers(outlier_rate[idx_noise_outlier], outlier_noise[idx_strength_outlier], current_seed);
        std::string noise_name = formatDouble(outlier_rate[idx_noise_outlier]) + "_" + formatDouble(outlier_noise[idx_strength_outlier]);
        // replace the "." in the string by "," (to avoid problems with the extension)
        std::replace( noise_name.begin(), noise_name.end(), '.', ',');
        std::string name = outputFilename.substr(0, outputFilename.size()-4) + "_outliers_" + noise_name + ".pts";

        std::ofstream ofs( name, std::ofstream::out );
        if ( !XYZOnly )
          ofs << "# x y z Gauss_Curvature Mean_Curvature nx ny nz k1 k2 d1x d1y d1z "
                "d2x d2y d2z\n";
        else 
          ofs << "# x y z\n";
        for ( auto i = 0u; i < nbPts; ++i )
        {
          ofs << noisy_pc[ i ][ 0 ] << " " << noisy_pc[ i ][ 1 ] << " " << noisy_pc[ i ][ 2 ];
            if ( !XYZOnly )
            {
              ofs << " " << exactGauss[ i ] << " " << exactMean[ i ] << " ";
              ofs << exactN[ i ][ 0 ] << " " << exactN[ i ][ 1 ] << " " << exactN[ i ][ 2 ] << " ";
              ofs << exactK1[ i ] << " " << exactK2[ i ] << " ";
              ofs << exactDir1[ i ][ 0 ] << " " << exactDir1[ i ][ 1 ] << " " << exactDir1[ i ][ 2 ] << " ";
              ofs << exactDir2[ i ][ 0 ] << " " << exactDir2[ i ][ 1 ] << " " << exactDir2[ i ][ 2 ];
            }
          ofs << std::endl;
        }
        ofs.close();
      }
    }    

  } 

  return 0;
}
