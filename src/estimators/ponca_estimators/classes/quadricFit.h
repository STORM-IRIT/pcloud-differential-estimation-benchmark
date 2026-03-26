#pragma once


#include <Eigen/Core>
#include <Ponca/Fitting>
#include "quadric.h"
#include "localFrame.h"

namespace Ponca
{

template < class DataPoint, class _NFilter, typename T >
class QuadricFitImpl : public T
{
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    PONCA_FITTING_DECLARE_MATRIX_TYPE

protected:
    enum
    {
        check = Base::PROVIDES_ALGEBRAIC_QUADRIC
    };

public:
    PONCA_EXPLICIT_CAST_OPERATORS(QuadricFitImpl, quadricFitImpl)

public:

    using VectorA = typename Base::VectorA;
    using MatrixA = typename Base::MatrixA;

    using Matrix32 = Eigen::Matrix<Scalar, 3, 2>;

    using Vector2 = Eigen::Matrix<Scalar, 2, 1>;
    using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;

// results
public:

    MatrixA m_matA;
    int m_minId;
    // EigenSolver of m_matA
    Eigen::EigenSolver<MatrixA> m_solver;

public:
    PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE

    PONCA_MULTIARCH inline FIT_RESULT constructTensor ();

protected:

    PONCA_MULTIARCH Matrix32 tangentPlane(const VectorType& n);

}; //class Quadric

template < class DataPoint, class _NFilter, typename T>
using QuadricFit = 
    QuadricFitImpl<DataPoint, _NFilter, 
        Quadric<DataPoint, _NFilter,
            Ponca::LocalFrame<DataPoint, _NFilter,
                Ponca::PrimitiveBase<DataPoint,_NFilter,T>>>>;


#include "quadricFit.hpp"
} //namespace Ponca
