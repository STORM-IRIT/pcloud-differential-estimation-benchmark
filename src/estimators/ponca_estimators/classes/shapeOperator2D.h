#pragma once

#include <Ponca/Fitting>

namespace Ponca
{

//!
//! 2D Shape operator fitting
//!
//! Fit the shape operator as a 2x2 matrix W defined by 
//!         W dp_i = dn_i
//! where 
//!         dp_i = P^T (p_i - p)
//!         dn_i = P^T (n_i - n) 
//! with 
//!     - p the 3D point of evaluation
//!     - n the 3D normal vector at p
//!     - p_i a 3D neighbor point
//!     - n_i the 3D normal vector at p_i
//!     - P the 3x2 matrix that contains 
//!       two 3D orthonormal tangent vectors as columns
//!
//! Least-squares fit of 
//!     E_W = sum_i w_i |W dp_i - dn_i|^2
//! 
//! Solution
//!     A W = B
//! where
//!     A = sum_i w_i dp_i dp_i^T
//!     B = sum_i w_i dn_i dp_i^T
//! 
//! See Equation 2 of
//! Robust statistical estimation of curvature on discretized surfaces
//! Kalogerakis et al. 2007
//! 
template < class DataPoint, class _NFilter, typename T >
class ShapeOperator2DFitImpl : public T
{
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    PONCA_FITTING_DECLARE_MATRIX_TYPE
    using Matrix32 = Eigen::Matrix<Scalar,3,2>;
    using Matrix22 = Eigen::Matrix<Scalar,2,2>;
    using Vector2 = Eigen::Matrix<Scalar,2,1>;

    enum { check = Base::PROVIDES_LOCAL_FRAME };
    static_assert(DataPoint::Dim == 3); // only valid in 3D

public:
    enum Status : int {
        Success           = 0,
        UnderDetermined, // 1
        NonInvertible,   // 2
        EigenSolverError // 3
    };

public:
    // input
    Matrix32 m_P;

    // computation data
    Matrix22 m_A;
    Matrix22 m_B;

    // results
    Matrix22 m_W;
    Scalar m_kmin, m_kmax;
    VectorType m_dmin, m_dmax;

    bool m_planeIsReady {false};

public:
    PONCA_EXPLICIT_CAST_OPERATORS(ShapeOperator2DFitImpl,shapeOperator2DFit)

    PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE
    
    PONCA_MULTIARCH inline Matrix22& dNormal() const {return m_W;}
    PONCA_MULTIARCH inline Matrix32 tangentPlane() const;

    //! \brief Returns an estimate of the mean curvature
    PONCA_MULTIARCH inline Scalar kMean() const { return (m_kmin + m_kmax) / Scalar(2); }

    //! \brief Returns an estimate of the Gaussian curvature
    PONCA_MULTIARCH inline Scalar GaussianCurvature() const { return m_kmin * m_kmax; };

    //! \brief Returns an estimate of the minimum curvature
    PONCA_MULTIARCH inline Scalar kmin() const { return m_kmin; };

    //! \brief Returns an estimate of the maximum curvature
    PONCA_MULTIARCH inline Scalar kmax() const { return m_kmax; }

    //! \brief Returns an estimate of the minimum curvature direction $
    PONCA_MULTIARCH inline VectorType kminDirection() const { return m_dmin; };

    //! \brief Returns an estimate of the maximum curvature direction
    PONCA_MULTIARCH inline VectorType kmaxDirection() const { return m_dmax; }

    PONCA_MULTIARCH inline VectorType primitiveGradient() const {
        VectorType covNormal = Base::primitiveGradient();
        VectorType meanNormal = Base::meanNormalVector();
        if (covNormal.dot(meanNormal) < Scalar(0))
            covNormal *= Scalar(-1);
        return covNormal;
    }

protected:
    // Solve X A = B 
    // by inverting A
    static bool solve(
        const Matrix22& A,
        const Matrix22& B,
        Matrix22& X);

}; //class ShapeOperatorFitImpl

template < class DataPoint, class _NFilter, typename T>
using ShapeOperator2DFit =
    ShapeOperator2DFitImpl<DataPoint, _NFilter, 
        Ponca::CovariancePlaneFitImpl<DataPoint, _NFilter,
            Ponca::CovarianceFitBase<DataPoint, _NFilter,
        // Ponca::MeanPlaneFitImpl<DataPoint, _NFilter,
                Ponca::MeanNormal<DataPoint, _NFilter,
                    Ponca::MeanPosition<DataPoint, _NFilter,
                        Ponca::LocalFrame<DataPoint, _NFilter,
                            Ponca::Plane<DataPoint, _NFilter,
                                Ponca::PrimitiveBase<DataPoint,_NFilter,T>>>>>>>>;


#include "shapeOperator2D.hpp"
} //namespace Ponca
