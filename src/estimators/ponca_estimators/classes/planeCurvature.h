#pragma once

#include <Ponca/src/Fitting/defines.h>

#include <Eigen/Dense>

namespace Ponca
{

template < class DataPoint, class _NFilter, typename T>
class CovPlaneCurvature : public T
{
PONCA_FITTING_DECLARE_DEFAULT_TYPES
PONCA_FITTING_DECLARE_MATRIX_TYPE

protected:
    enum { Check = Base::PROVIDES_LOCAL_FRAME };

public:
    using Solver = Eigen::SelfAdjointEigenSolver<MatrixType>;
protected:

    using Matrix32 = Eigen::Matrix<Scalar, 3, 2>;
    using Vector2 = Eigen::Matrix<Scalar, 2, 1>;

    VectorType   m_dmin;
    VectorType   m_dmax;

    Scalar       m_kmin;
    Scalar       m_kmax;


public:
    PONCA_EXPLICIT_CAST_OPERATORS(CovPlaneCurvature,covPlaneCurvature)
    // PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE

    PONCA_MULTIARCH void init()
    {
        Base::init();
        m_dmin = VectorType::Zero();
        m_dmax = VectorType::Zero();
        m_kmin = Scalar(0);
        m_kmax = Scalar(0);
    }

    PONCA_MULTIARCH bool addLocalNeighbor( Scalar w, const VectorType &localQ, const DataPoint &attributes ){
        return Base::addLocalNeighbor(w, localQ, attributes);
    }

    PONCA_MULTIARCH FIT_RESULT finalize()
    {
        auto res = Base::finalize();
        if ( res != FIT_RESULT::STABLE )
            return res;
        m_dmin = Base::solver().eigenvectors().col(1);
        m_dmax = Base::solver().eigenvectors().col(2);
        m_kmin = ( Base::solver().eigenvalues()(1) / Base::solver().eigenvalues()(0) ) * Base::getNeighborFilter().evalScale();
        m_kmax = ( Base::solver().eigenvalues()(2) / Base::solver().eigenvalues()(0) ) * Base::getNeighborFilter().evalScale();
        return res;
    }

    //! \brief Returns an estimate of the mean curvature
    PONCA_MULTIARCH inline Scalar kMean() const { return ( m_kmin + m_kmax ) / Scalar(2); }

    //! \brief Returns an estimate of the Gaussian curvature
    PONCA_MULTIARCH inline Scalar GaussianCurvature() const { return m_kmin * m_kmax; }

    //! \brief Returns an estimate of the minimum curvature
    PONCA_MULTIARCH inline Scalar kmin() const { return m_kmin; }

    //! \brief Returns an estimate of the maximum curvature
    PONCA_MULTIARCH inline Scalar kmax() const { return m_kmax; }

    //! \brief Returns an estimate of the minimum curvature direction $
    PONCA_MULTIARCH inline VectorType kminDirection() const { return m_dmin; }

    //! \brief Returns an estimate of the maximum curvature direction
    PONCA_MULTIARCH inline VectorType kmaxDirection() const { return m_dmax; }

    //! \brief Orthogonal projecting on the patch, such that h = f(u,v)
    PONCA_MULTIARCH inline VectorType project (const VectorType& _q) const { return Base::plane().project(_q); }

    PONCA_MULTIARCH inline VectorType primitiveGradient() const { return Base::plane().primitiveGradient(); }

    PONCA_MULTIARCH inline Matrix32 tangentPlane (  ) const { 
        Matrix32 B;
        B.col(0) = Base::getFrameU();
        B.col(1) = Base::getFrameV();
        return B;
     }
};

  template < class DataPoint, class _NFilter, typename T>
class MeanPlaneCurvature : public T
{
PONCA_FITTING_DECLARE_DEFAULT_TYPES
PONCA_FITTING_DECLARE_MATRIX_TYPE

protected:
    enum { Check = Base::PROVIDES_LOCAL_FRAME };

protected:

    using Matrix32 = Eigen::Matrix<Scalar, 3, 2>;
    using Vector2 = Eigen::Matrix<Scalar, 2, 1>;

    VectorType   m_dmin;
    VectorType   m_dmax;

    Scalar       m_kmin;
    Scalar       m_kmax;


public:
    PONCA_EXPLICIT_CAST_OPERATORS(MeanPlaneCurvature,meanPlaneCurvature)

    PONCA_MULTIARCH void init()
    {
        Base::init();
        m_dmin = VectorType::Zero();
        m_dmax = VectorType::Zero();
        m_kmin = Scalar(0);
        m_kmax = Scalar(0);
    }

    PONCA_MULTIARCH bool addLocalNeighbor( Scalar w, const VectorType &localQ, const DataPoint &attributes ){
        return Base::addLocalNeighbor(w, localQ, attributes);
    }

    PONCA_MULTIARCH FIT_RESULT finalize()
    {
        auto res = Base::finalize();
        if ( res != FIT_RESULT::STABLE )
            return res;

        Base::computeFrameFromNormalVector(Base::plane().primitiveGradient());
      return res;
    }

    PONCA_MULTIARCH inline Scalar kMean() const { return ( m_kmin + m_kmax ) / Scalar(2); }
    PONCA_MULTIARCH inline Scalar GaussianCurvature() const { return m_kmin * m_kmax; }
    PONCA_MULTIARCH inline Scalar kmin() const { return m_kmin; }
    PONCA_MULTIARCH inline Scalar kmax() const { return m_kmax; }
    PONCA_MULTIARCH inline VectorType kminDirection() const { return Base::getFrameU(); }
    PONCA_MULTIARCH inline VectorType kmaxDirection() const { return Base::getFrameV(); }
    PONCA_MULTIARCH inline VectorType project (const VectorType& _q) const { return Base::plane().project(_q); }
    PONCA_MULTIARCH inline VectorType primitiveGradient() const { return Base::plane().primitiveGradient(); }
    PONCA_MULTIARCH inline Matrix32 tangentPlane (  ) const {
        Matrix32 B;
        B.col(0) = Base::getFrameU();
        B.col(1) = Base::getFrameV();
        return B;
     }
};

template < class DataPoint, class _NFilter, typename T>
    using CovPlaneCurvatureFit =
        Ponca::CovPlaneCurvature<DataPoint, _NFilter,
            Ponca::CovariancePlaneFitImpl<DataPoint, _NFilter,
                Ponca::CovarianceFitBase<DataPoint, _NFilter,
                    Ponca::MeanPosition<DataPoint, _NFilter,
                        Ponca::LocalFrame<DataPoint, _NFilter,
                            Ponca::Plane<DataPoint, _NFilter,
                                Ponca::PrimitiveBase<DataPoint,_NFilter,T>>>>>>>;

  template < class DataPoint, class _NFilter, typename T>
    using MeanPlaneCurvatureFit =
        Ponca::MeanPlaneCurvature<DataPoint, _NFilter,
          Ponca::MeanPlaneFitImpl<DataPoint, _NFilter,
            Ponca::MeanNormal<DataPoint, _NFilter,
                Ponca::MeanPosition<DataPoint, _NFilter,
                      Ponca::LocalFrame<DataPoint, _NFilter,
                            Ponca::Plane<DataPoint, _NFilter,
                                Ponca::PrimitiveBase<DataPoint,_NFilter,T>>>>>>>;
} //namespace Ponca
