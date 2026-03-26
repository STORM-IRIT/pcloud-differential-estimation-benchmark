#pragma once

#include <utility>
#include <nanoflann.hpp>
#include <Eigen/Dense>
#include "KDTreeVectorOfVectorsAdaptor.h"

template<typename _Scalar>
struct KDTree
{
  using Scalar = _Scalar;
  using num_t = _Scalar;
  using Point = Eigen::Vector<_Scalar, 3>;
  using Index = size_t;

  using MyTree = KDTreeVectorOfVectorsAdaptor<std::vector<Point>, Scalar>;

  MyTree * mat_index;

  KDTree( std::vector<Point> & samples )
  {
    mat_index = new MyTree( 3, samples, 10 );
  }

  std::vector<std::pair<Index, Scalar>> knnSearch( Point & query, size_t k )
  {
    std::vector<size_t> ret_indices( k );
    std::vector<Scalar> out_dists_sqr( k );
    Scalar query_pt[ 3 ] = { query( 0 ), query( 1 ), query( 2 ) };

    nanoflann::KNNResultSet<Scalar> resultSet( k );

    resultSet.init( &ret_indices[ 0 ], &out_dists_sqr[ 0 ] );
    mat_index->index->findNeighbors( resultSet, &query_pt[ 0 ], nanoflann::SearchParams( 10 ) );

    std::vector<std::pair<Index, Scalar>> ret_indices_out;
    for ( size_t i = 0; i < k; i++ )
    {
      ret_indices_out.push_back( std::make_pair( ret_indices[ i ], out_dists_sqr[ i ] ) );
    }

    return ret_indices_out;
  }

  std::vector<std::pair<Index, Scalar>> radiusSearch( Point & query, Scalar radius )
  {
    std::vector<std::pair<Index, Scalar>> indices_dists;
    nanoflann::RadiusResultSet<num_t, size_t> resultSet( radius, indices_dists );
    Scalar query_pt[ 3 ] = { query( 0 ), query( 1 ), query( 2 ) };

    mat_index->index->findNeighbors( resultSet, &query_pt[ 0 ], nanoflann::SearchParams( 10, 0.0f, true ) );

    return indices_dists;
  }
};
