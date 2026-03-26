#pragma once

#include <Ponca/src/Fitting/defines.h>
#include <Eigen/Dense>

namespace Ponca
{

template < class DataPoint, class _NFilter, typename T>
class NormalCovariance3D : public T
{
PONCA_FITTING_DECLARE_DEFAULT_TYPES
PONCA_FITTING_DECLARE_MATRIX_TYPE

protected:
    enum { Check = Base::PROVIDES_LOCAL_FRAME };

public:
    using Matrix32 = Eigen::Matrix<Scalar, 3, 2>;
    using Matrix22 = Eigen::Matrix<Scalar, 2, 2>;
    using Solver = Eigen::SelfAdjointEigenSolver<Matrix22>;
protected:

    MatrixType   m_cov;
    Solver       m_solver;
    VectorType   m_normal_centroid = VectorType::Zero();
    bool         m_planeIsReady = false;
    int          m_minDirIndex = 0;
    int          m_maxDirIndex = 0;
    int          m_normalIndex = 0;

    Matrix32     m_P;
    Matrix22     m_W;

public:
    PONCA_EXPLICIT_CAST_OPERATORS(NormalCovariance3D,normalCovariance3D)
    PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE

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
        return _q;
    }

    PONCA_MULTIARCH inline VectorType primitiveGradient() const
    {
        VectorType normal = Base::normal();
        normal.normalize();
        return normal;
        // return Base::normal();
    }


private:

    PONCA_MULTIARCH inline void setTangentPlane()
    {
        m_P.col(0) = Base::getFrameU();
        m_P.col(1) = Base::getFrameV();
    }

};

/// \brief Helper alias for NormalCovariance3DFit on points
template < class DataPoint, class _NFilter, typename T>
    using NormalCovariance3DFit =
        Ponca::NormalCovariance3D<DataPoint, _NFilter,
            Ponca::CovariancePlaneFitImpl<DataPoint, _NFilter,
                Ponca::CovarianceFitBase<DataPoint, _NFilter,
            // Ponca::MeanPlaneFitImpl<DataPoint, _NFilter,
                    Ponca::MeanNormal<DataPoint, _NFilter,
                        Ponca::MeanPosition<DataPoint, _NFilter,
                            Ponca::LocalFrame<DataPoint, _NFilter,
                                Ponca::Plane<DataPoint, _NFilter,
                                    Ponca::PrimitiveBase<DataPoint,_NFilter,T>>>>>>>>;


#include "normalCovariance3DFit.hpp"

} //namespace Ponca
