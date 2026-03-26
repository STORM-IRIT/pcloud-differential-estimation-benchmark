#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <random>

#include <Eigen/Dense>
#include <spdlog/spdlog.h>
#include <DGtal/math/Statistic.h>

#include "Readers.h"
#include "Writters.h"


using Scalar = double;
// using Scalar = float;
unsigned int seed = 0;
std::random_device random_seed; 
std::default_random_engine generator;

template <typename _Scalar>
struct PointCloudDiff
{
  /// Types
  typedef Eigen::VectorX<_Scalar>      DenseVector;
  typedef Eigen::MatrixX<_Scalar>      DenseMatrix; 
  typedef Eigen::Vector<_Scalar, 3>    Point;
  typedef Eigen::Matrix<_Scalar, 3, 3> Tensor;

  // PointCloudDiff()  = default;
  ~PointCloudDiff() = default;

  PointCloudDiff () {
    errorNormal = DGtal::Statistic<_Scalar> (true);
    errorPos = DGtal::Statistic<_Scalar> (true);
    errorMean = DGtal::Statistic<_Scalar> (true);
    errorGauss = DGtal::Statistic<_Scalar> (true);
    errorK1 = DGtal::Statistic<_Scalar> (true);
    errorK2 = DGtal::Statistic<_Scalar> (true);
    errorD1 = DGtal::Statistic<_Scalar> (true);
    errorD2 = DGtal::Statistic<_Scalar> (true);
    errorShapeIndex = DGtal::Statistic<_Scalar> (true);
    statTimings = DGtal::Statistic<_Scalar> (true);
    statNeighbors = DGtal::Statistic<_Scalar> (true);
  }

  PointCloudDiff( const PointCloudDiff & apc ){
    points = apc.points;
    normals = apc.normals;
    k1 = apc.k1;
    k2 = apc.k2;
    mean = apc.mean;
    gauss = apc.gauss;
    v1 = apc.v1;
    v2 = apc.v2;
    errorNormal = apc.errorNormal;
    errorPos = apc.errorPos;
    errorMean = apc.errorMean;
    errorGauss = apc.errorGauss;
    errorK1 = apc.errorK1;
    errorK2 = apc.errorK2;
    errorD1 = apc.errorD1;
    errorD2 = apc.errorD2;
    errorShapeIndex = apc.errorShapeIndex;
    statTimings = apc.statTimings;
    statNeighbors = apc.statNeighbors;
    non_stable_idx = apc.non_stable_idx;
  }

  void loadPointCloud ( std::string & input )
  {
    IO::loadPointCloud<PointCloudDiff<_Scalar>, _Scalar>( *this, input );
  }

  void savePointCloud ( std::string & filename )
  {
    IO::savePointCloud<PointCloudDiff<_Scalar>>( *this, filename );
  }
  
  void savePointCloudAsPositions ( std::string & filename )
  {
    IO::savePointCloudAsPositions<PointCloudDiff<_Scalar>>( *this, filename );
  }

  void savePointCloudAsErrors ( std::string & filename )
  {
    IO::savePointCloudAsErrors<PointCloudDiff<_Scalar>, _Scalar>( *this, filename );
  }

  void init_zeros( int nVerts )
  {
    k1 = std::vector<_Scalar>( nVerts, 0 );
    k2 = std::vector<_Scalar>( nVerts, 0 );
    mean = std::vector<_Scalar>( nVerts, 0 );
    gauss = std::vector<_Scalar>( nVerts, 0 );
    v1 = std::vector<Point>( nVerts, Point( 0, 0, 0 ) );
    v2 = std::vector<Point>( nVerts, Point( 0, 0, 0 ) );
    normals = std::vector<Point>( nVerts, Point( 0, 0, 0 ) );
    points = std::vector<Point>( nVerts, Point( 0, 0, 0 ) );
  }

  std::string getStats ( std::string & name, _Scalar radius, _Scalar noise_position, _Scalar noise_normal, _Scalar flip_ratio, bool absolute_err = true )
  {
    return IO::getStats<PointCloudDiff<_Scalar>>( *this, name, radius, noise_position, noise_normal, flip_ratio, absolute_err );
  }

  void addOutliers( _Scalar ratio, _Scalar noise_position )
  {
    generator.seed( ( seed == 0 ) ? random_seed() : seed );
    std::uniform_real_distribution<_Scalar> dis( 0, 1 );
    std::normal_distribution<_Scalar> dis_noise( 0, noise_position );
    for ( auto & p : points )
    {
      if (dis(generator) < ratio){
        p[ 0 ] += dis_noise( generator );
        p[ 1 ] += dis_noise( generator );
        p[ 2 ] += dis_noise( generator );
      }
    }
  }

  void flipNormal(_Scalar ratio){
    generator.seed( ( seed == 0 ) ? random_seed() : seed );
    std::uniform_real_distribution<_Scalar> dis( 0, 1 );
    for ( auto & n : normals )
    {
      if (dis(generator) < ratio){
        n[ 0 ] *= -1;
        n[ 1 ] *= -1;
        n[ 2 ] *= -1;
      }
    }
  }

  void addNoisePosition( _Scalar noise_position )
  {
    generator.seed( ( seed == 0 ) ? random_seed() : seed );
    std::normal_distribution<_Scalar> dis( 0, noise_position );
    for ( auto & v : points )
    {
      v[ 0 ] += dis( generator );
      v[ 1 ] += dis( generator );
      v[ 2 ] += dis( generator );
    }
  }

  void addNoiseNormal( _Scalar noise_normal )
  {
    generator.seed( ( seed == 0 ) ? random_seed() : seed );
    std::normal_distribution<_Scalar> dis( 0, noise_normal );
    for ( auto & n : normals )
    {
      n[ 0 ] += dis( generator );
      n[ 1 ] += dis( generator );
      n[ 2 ] += dis( generator );
      _Scalar norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
      n[0] /= norm;
      n[1] /= norm;
      n[2] /= norm;
    }
  }

  void check_k1_k2 ()
  {
    #pragma omp parallel for
    for ( int i = 0; i < points.size(); i++ )
    {
      if ( k1[ i ]  >  k2[ i ] )
      {
        std::swap (k1[i], k2[i]);
        std::swap (v1[i], v2[i]);
      }
    }
  }

  void switch_to_abs ()
  {
    #pragma omp parallel for
    for ( int i = 0; i < points.size(); i++ )
    {
      // particular treatment for kmean and kgauss
      mean[i] = std::abs(mean[i]);
      gauss[i] = std::abs(gauss[i]);

      k1[i] = std::abs(k1[i]);
      k2[i] = std::abs(k2[i]);

      if ( k1[ i ]  >  k2[ i ] )
      {
        std::swap (k1[i], k2[i]);
        std::swap (v1[i], v2[i]);
      }

    }
  }

  _Scalar squared_error( _Scalar error ) {
    return std::pow( ( error ), 2 );
  }

  _Scalar abs_error( _Scalar error ) {
    return std::fabs( error );
  }

  template <typename Functor>
  _Scalar compute_scalar_error ( _Scalar val, _Scalar other, Functor metric ) {
    return metric( val - other );
  }

  _Scalar compute_radian_error ( Point vec, Point other, bool oriented = false ) {
    _Scalar cos_angle = vec.dot( other ) / ( vec.norm() * other.norm() );
    if ( !oriented && cos_angle < 0. ) {
      cos_angle = -cos_angle;
    }
    cos_angle = ( cos_angle > 1.0 ) ? 1.0 : cos_angle;
    cos_angle = ( cos_angle < -1.0 ) ? -1.0 : cos_angle;
    return std::acos( cos_angle );
  }

  template <typename Functor>
  _Scalar compute_vector_error( Point vec, Point other, Functor metric, bool oriented = false ) {
    _Scalar radian_err = compute_radian_error ( vec, other, oriented );
    return metric ( radian_err * 180.0 / M_PI ); // degree
  }

  template <typename Functor>
  _Scalar compute_point_error ( Point p, Point other, Functor metric ) {
    return metric( ( p - other ).norm() );
  }

  void correct_vectors ( int idx ) {
    if ( v1[ idx ].norm() < 1e-6 )
    {
      v1[ idx ] = Point( 1, 0, 0 );
      v2[ idx ] = Point( 0, 1, 0 );
    }
    if ( normals[ idx ].norm() < 1e-6 )
    {
      normals[ idx ] = Point( 1, 0, 0 );
    }
  }

  template <typename Functor>
  void compare_at ( const PointCloudDiff & other, int i, Functor metric ){

    errorPos.addValue( compute_point_error( points[ i ], other.points[ i ], metric ) );
    errorMean.addValue( compute_scalar_error( mean[ i ], other.mean[ i ], metric ) );
    errorGauss.addValue( compute_scalar_error( gauss[ i ], other.gauss[ i ], metric ) );

    errorK1.addValue( compute_scalar_error( k1[ i ], other.k1[ i ], metric ) );
    errorK2.addValue( compute_scalar_error( k2[ i ], other.k2[ i ], metric ) );

    errorNormal.addValue( compute_vector_error( normals[ i ], other.normals[ i ], metric, oriented ) );

    if ( std::fabs( other.k1[ i ] - other.k2[ i ] ) > 1e-6 )
    {
      errorD1.addValue( compute_vector_error( v1[ i ], other.v1[ i ], metric ) );
      errorD2.addValue( compute_vector_error( v2[ i ], other.v2[ i ], metric ) );

      _Scalar sIndex       = ( 2.0 / M_PI ) * std::atan( ( k1[ i ] + k2[ i ] ) / ( k1 [ i ] - k2[ i ] ) );
      _Scalar sIndex_other = ( 2.0 / M_PI ) * std::atan( ( other.k1[ i ] + other.k2[ i ] ) / ( other.k1 [ i ] - other.k2[ i ] ) );
      errorShapeIndex.addValue( compute_scalar_error( sIndex, sIndex_other, metric ) );
    }
    else
    {
      errorD1.addValue( 0. );
      errorD2.addValue( 0. );
      errorShapeIndex.addValue( 0. );
    }
  } 

  template <typename Functor>
  void estimate_at ( int i, Functor metric ){

    errorPos.addValue( 0 );
    errorMean.addValue( metric( metric(mean[i]) ));
    errorGauss.addValue( metric(gauss[ i ]) );

    errorK1.addValue( metric(k1[ i ]) );
    errorK2.addValue( metric(k2[ i ]) );

    errorNormal.addValue( 0 );

    errorD1.addValue( 0 );
    errorD2.addValue( 0 );
    _Scalar sIndex       = ( 2.0 / M_PI ) * std::atan( ( k1[ i ] + k2[ i ] ) / ( k1 [ i ] - k2[ i ] ) );
    if ( std::fabs( k1[ i ] - k2[ i ] ) > 1e-6 ){
      errorShapeIndex.addValue( metric( sIndex ) );
    }
  } 

  void compute_unique_error ( const PointCloudDiff & other, int unique_idx, bool absolute_err = true ){
      if ( non_stable_idx[ unique_idx ] == 1 )
      {
        statNeighbors.addValue( 1 );
        return;
      }
      correct_vectors( unique_idx );
      compare_at( other, unique_idx, [this, &absolute_err] ( _Scalar error ) { return ( absolute_err ) ? abs_error ( error ) : squared_error ( error ); } );
  }

  void compare( const PointCloudDiff & other, bool absolute_err = true, bool statsAsEstimations = false )
  {
    for ( int i = 0; i < points.size(); i++ )
    {
      // if ( non_stable_idx.size() > 0 && non_stable_idx[ i ] == 1 )
      // {
      //   statNeighbors.addValue( 1 );
      //   continue;
      // }
      correct_vectors( i );
      if ( statsAsEstimations )
      {
        estimate_at( i, 
          [this, &absolute_err] ( _Scalar estim ) 
          { 
            return ( absolute_err ) ? abs_error ( estim ) : estim; 
          } 
        );
      }
      else {
        compare_at( other, i, 
          [this, &absolute_err] ( _Scalar error ) 
          { 
            return ( absolute_err ) ? abs_error ( error ) : squared_error ( error ); 
          } 
        );
      }
      
    }
  }

  bool oriented = true;

  std::vector<Point> points;
  std::vector<Point> normals;
  std::vector<_Scalar> gauss, mean, k1, k2;
  std::vector<Point> v1;
  std::vector<Point> v2;
  std::vector<unsigned int> non_stable_idx;
  
  DGtal::Statistic<_Scalar> errorShapeIndex;

  DGtal::Statistic<_Scalar> errorPos;
  DGtal::Statistic<_Scalar> errorNormal;
  DGtal::Statistic<_Scalar> errorMean;
  DGtal::Statistic<_Scalar> errorGauss;
  DGtal::Statistic<_Scalar> errorK1;
  DGtal::Statistic<_Scalar> errorK2;
  DGtal::Statistic<_Scalar> errorD1;
  DGtal::Statistic<_Scalar> errorD2;
  DGtal::Statistic<_Scalar> statTimings;
  DGtal::Statistic<_Scalar> statNeighbors;
};
