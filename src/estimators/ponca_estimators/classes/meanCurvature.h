#pragma once

#include <Ponca/src/Fitting/defines.h>
#include <Eigen/Dense>

namespace Ponca
{

template < class DataPoint, class _NFilter, typename T>
class MeanCurvature : public T
{
PONCA_FITTING_DECLARE_DEFAULT_TYPES
PONCA_FITTING_DECLARE_MATRIX_TYPE

protected:
    enum { Check = Base::PROVIDES_LOCAL_FRAME };

public:
    using Solver = Eigen::SelfAdjointEigenSolver<MatrixType>;
protected:

    VectorType   m_pos;

public:
    PONCA_EXPLICIT_CAST_OPERATORS(MeanCurvature,meanCurvature)
    // PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE

    PONCA_MULTIARCH void init()
    {
        Base::init();
        m_pos = Base::getNeighborFilter().evalPos();
    }

    PONCA_MULTIARCH bool addLocalNeighbor( Scalar w, const VectorType &localQ, const DataPoint &attributes ){
        auto res = Base::addLocalNeighbor(w, localQ, attributes);
        m_pos += w * attributes.pos();
        return res;

    }

    PONCA_MULTIARCH FIT_RESULT finalize()
    {
        auto res = Base::finalize();
        
        m_pos /= Base::getWeightSum();

        return res;
    }

    //! \brief Returns an estimate of the mean curvature
    PONCA_MULTIARCH inline Scalar kMean() const { 
        PONCA_MULTIARCH_STD_MATH(pow);
        PONCA_MULTIARCH_STD_MATH(abs);
        Scalar four = Scalar(4);
        Scalar two = Scalar(2);
        VectorType proj = m_pos - Base::getNeighborFilter().evalPos();
        Scalar dist = proj.norm();
        return four * dist / ( pow( Base::getNeighborFilter().evalScale(), two ) );
}

    //! \brief Returns an estimate of the Gaussian curvature
    PONCA_MULTIARCH inline Scalar GaussianCurvature() const { return Scalar(0); }

    //! \brief Returns an estimate of the minimum curvature
    PONCA_MULTIARCH inline Scalar kmin() const { return Scalar(0); }

    //! \brief Returns an estimate of the maximum curvature
    PONCA_MULTIARCH inline Scalar kmax() const { return Scalar(0); }

    //! \brief Returns an estimate of the minimum curvature direction $
    PONCA_MULTIARCH inline VectorType kminDirection() const { return VectorType::Zero(); }

    //! \brief Returns an estimate of the maximum curvature direction
    PONCA_MULTIARCH inline VectorType kmaxDirection() const { return VectorType::Zero(); }

    //! \brief Orthogonal projecting on the patch, such that h = f(u,v)
    PONCA_MULTIARCH inline VectorType project (const VectorType& _q) const
    {
        return Base::project(_q);
    }

    PONCA_MULTIARCH inline VectorType primitiveGradient() const
    {
        return Base::normal();
    }
};

/// \brief Helper alias for MeanCurvatureFit on points
template < class DataPoint, class _NFilter, typename T>
    using MeanCurvatureFit =
        Ponca::MeanCurvature<DataPoint, _NFilter,
            Ponca::MeanPlaneFitImpl<DataPoint, _NFilter,
                    Ponca::MeanNormal<DataPoint, _NFilter,
                        Ponca::MeanPosition<DataPoint, _NFilter,
                            Ponca::LocalFrame<DataPoint, _NFilter,
                                Ponca::Plane<DataPoint, _NFilter,
                                    Ponca::PrimitiveBase<DataPoint,_NFilter,T>>>>>>>;

} //namespace Ponca
