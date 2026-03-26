#pragma once

#include <Ponca/src/Fitting/defines.h>
#include "Ponca/src/Fitting/covariancePlaneFit.h"

#include <Eigen/Dense>

namespace Ponca
{

template < class DataPoint, class _NFilter, typename T>
class NormalWeingarten : public T
{
PONCA_FITTING_DECLARE_DEFAULT_TYPES

protected:
    enum { Check = Base::PROVIDES_LOCAL_FRAME };

public:
    using Matrix2       = Eigen::Matrix<Scalar, 2, 2>;
    using Matrix3       = Eigen::Matrix<Scalar, 3, 3>;

    using Vector2       = Eigen::Matrix<Scalar, 2, 1>;
    using Vector3       = Eigen::Matrix<Scalar, 3, 1>;
    using Vector4       = Eigen::Matrix<Scalar, 4, 1>;

    using Matrix32      = Eigen::Matrix<Scalar, 3, 2>;
    using Solver = Eigen::SelfAdjointEigenSolver<Matrix2>;
protected:

    Matrix3   m_XtX;     
    Vector3   m_XtY1;        
    Vector3   m_XtY2;

    Scalar m_mu1;
    Scalar m_mu2;
    Vector4 m_w_tilde;

    Scalar    m_kmin; 
    Scalar    m_kmax;

    VectorType m_dmin;
    VectorType m_dmax;
    
    Matrix32  m_P;

    VectorType m_normal;
    VectorType m_normal_no_proj;

    bool m_planeIsReady {false};

private:
    
    PONCA_MULTIARCH inline void computeCurvature();

    PONCA_MULTIARCH inline Vector4 projectOnto4DPlane(const Vector4& point, const Vector4& normal, const Scalar& offset) const;

public:
    PONCA_EXPLICIT_CAST_OPERATORS(NormalWeingarten,normalWeingarten)
    PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE

    //! \brief Approximation of the scalar field gradient at \f$ \mathbf{q} \f$
    PONCA_MULTIARCH inline VectorType primitiveGradient (const VectorType& _q) const { return primitiveGradient(); };

    /*! \brief Approximation of the scalar field gradient at the evaluation point */
    PONCA_MULTIARCH inline VectorType primitiveGradient () const {
        return m_normal;
    }
    
    //! \brief Returns an estimate of the mean curvature
    PONCA_MULTIARCH inline Scalar kMean() const;

    //! \brief Returns an estimate of the Gaussian curvature
    PONCA_MULTIARCH inline Scalar GaussianCurvature() const;

    //! \brief Returns an estimate of the minimum curvature
    PONCA_MULTIARCH inline Scalar kmin() const;

    //! \brief Returns an estimate of the maximum curvature
    PONCA_MULTIARCH inline Scalar kmax() const;

    //! \brief Returns an estimate of the minimum curvature direction $
    PONCA_MULTIARCH inline VectorType kminDirection() const;

    //! \brief Returns an estimate of the maximum curvature direction
    PONCA_MULTIARCH inline VectorType kmaxDirection() const;


    //! \brief Orthogonal projecting on the patch, such that h = f(u,v)
    PONCA_MULTIARCH inline VectorType project (const VectorType& _q) const
    {
        VectorType x = Base::worldToLocalFrame(_q);
        *(x.data()) = Scalar(0);
        return Base::localFrameToWorld(x);
    }

    PONCA_MULTIARCH inline void setTangentPlane()
    {
        m_P.col(0) = Base::getFrameU();
        m_P.col(1) = Base::getFrameV();
    }

    PONCA_MULTIARCH inline Scalar f_u () const
    {
        return - m_normal(0) * sqrt( Scalar(1.0) - m_normal(0) * m_normal(0) - m_normal(1) * m_normal(1) );
    }

    PONCA_MULTIARCH inline Scalar f_v () const
    {
        return - m_normal(1) * sqrt( Scalar(1.0) - m_normal(0) * m_normal(0) - m_normal(1) * m_normal(1) );
    }

    PONCA_MULTIARCH inline Scalar f_uu () const
    {
        Scalar X = m_normal_no_proj(0);
        Scalar Y = m_normal_no_proj(1);
        Scalar Xu = m_w_tilde(0);
        Scalar Xv = m_w_tilde(1);
        Scalar Yu = m_w_tilde(2);
        Scalar Yv = m_w_tilde(3);

        return - ( X * Y * Yu + ( Scalar(1.0) - Y * Y ) * Xu ) * pow( 1.0 - X * X - Y * Y, Scalar(3.0)/Scalar(2.0) );
    }
    
    PONCA_MULTIARCH inline Scalar f_vu () const
    {
        Scalar X = m_normal_no_proj(0);
        Scalar Y = m_normal_no_proj(1);
        Scalar Xu = m_w_tilde(0);
        Scalar Xv = m_w_tilde(1);
        Scalar Yu = m_w_tilde(2);
        Scalar Yv = m_w_tilde(3);

        return - ( X * Y * Xu + ( Scalar(1.0) - X * X ) * Yu ) * pow( 1.0 - X * X - Y * Y, Scalar(3.0)/Scalar(2.0) );
    }

    PONCA_MULTIARCH inline Scalar f_uv () const
    {
        Scalar X = m_normal_no_proj(0);
        Scalar Y = m_normal_no_proj(1);
        Scalar Xu = m_w_tilde(0);
        Scalar Xv = m_w_tilde(1);
        Scalar Yu = m_w_tilde(2);
        Scalar Yv = m_w_tilde(3);

        return - ( X * Y * Yv + ( Scalar(1.0) - Y * Y ) * Xv ) * pow( 1.0 - X * X - Y * Y, Scalar(3.0)/Scalar(2.0) );
    }

    PONCA_MULTIARCH inline Scalar f_vv () const
    {
        Scalar X = m_normal_no_proj(0);
        Scalar Y = m_normal_no_proj(1);
        Scalar Xu = m_w_tilde(0);
        Scalar Xv = m_w_tilde(1);
        Scalar Yu = m_w_tilde(2);
        Scalar Yv = m_w_tilde(3);

        return - ( X * Y * Xv + ( Scalar(1.0) - X * X ) * Yv ) * pow( 1.0 - X * X - Y * Y, Scalar(3.0)/Scalar(2.0) );
    }

};

/// \brief Helper alias for Covariance2DFit on points
//! [Covariance2DFit Definition]
template < class DataPoint, class _NFilter, typename T>
    using NormalWeingartenFit =
        Ponca::NormalWeingarten<DataPoint, _NFilter,
            Ponca::MeanPlaneFitImpl<DataPoint, _NFilter,
                    Ponca::MeanNormal<DataPoint, _NFilter,
                        Ponca::MeanPosition<DataPoint, _NFilter,
                            Ponca::LocalFrame<DataPoint, _NFilter,
                                Ponca::Plane<DataPoint, _NFilter,
                                    Ponca::PrimitiveBase<DataPoint,_NFilter,T>>>>>>>;


#include "normalW.hpp"

} //namespace Ponca
