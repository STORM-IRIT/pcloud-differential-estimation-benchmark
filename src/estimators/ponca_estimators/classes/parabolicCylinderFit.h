#pragma once

#include "parabolicCylinder.h"
#include "localFrame.h"
#include <Ponca/Fitting>
#include <Eigen/Dense>

namespace Ponca
{

template < class DataPoint, class _NFilter, typename T >
class ParabolicCylinderFitImpl: public T
{
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    // PONCA_FITTING_DECLARE_MATRIX_TYPE

protected:
    enum
    {
        Check = Base::PROVIDES_PRIMITIVE_BASE &&
                Base::PROVIDES_LOCAL_FRAME &&
                Base::PROVIDES_PARABOLIC_CYLINDER, /*!< \brief Requires PrimitiveBase and plane*/
    };

public:
    using SampleMatrix  = Eigen::Matrix<Scalar, 7, 7>;
    using SampleVector  = Eigen::Matrix<Scalar, 7, 1>;
    using Vector7       = Eigen::Matrix<Scalar, 7, 1>;

    // Todo : remove when differential properties are implemented in other class
    using Matrix2       = Eigen::Matrix<Scalar, 2, 2>;
    using Vector2       = Eigen::Matrix<Scalar, 2, 1>;
// results
protected:

    SampleMatrix m_A_cov;
    SampleVector m_F_cov;

    Vector2      m_vector_uq;

    bool m_planeIsReady {false};
public:
    
    PONCA_EXPLICIT_CAST_OPERATORS(ParabolicCylinderFitImpl,parabolicCylinderFitImpl)
    PONCA_FITTING_DECLARE_INIT_ADD_FINALIZE

private:

    PONCA_MULTIARCH inline void    m_fitting_process         ();

    PONCA_MULTIARCH inline void    m_ellipsoid_fitting       ();
    PONCA_MULTIARCH inline void    m_parabolic_fitting       ();
    
}; //class ParabolicCylinderFitImpl

/// \brief Helper alias for ParabolicCylinder fitting on points
//! [ParabolicCylinderFit Definition]
template < class DataPoint, class _NFilter, typename T>
    using BaseParabolicCylinderFit =
    Ponca::ParabolicCylinderFitImpl<DataPoint, _NFilter,
        Ponca::ParabolicCylinder<DataPoint, _NFilter,
            Ponca::CovariancePlaneFitImpl<DataPoint, _NFilter,
                Ponca::CovarianceFitBase<DataPoint, _NFilter,
                    Ponca::MeanNormal<DataPoint, _NFilter,
                        Ponca::MeanPosition<DataPoint, _NFilter,
                            Ponca::LocalFrame<DataPoint, _NFilter,
                                Ponca::Plane<DataPoint, _NFilter,
                                    Ponca::PrimitiveBase<DataPoint,_NFilter,T>>>>>>>>>;
                                            
// Base Method PC-MLS cov + oriented
template < class DataPoint, class _NFilter, typename T>
    using BaseOrientedParabolicCylinderFit =
        Ponca::ParabolicCylinderFitImpl<DataPoint, _NFilter,
            Ponca::ParabolicCylinder<DataPoint, _NFilter,
                Ponca::MeanPlaneFitImpl<DataPoint, _NFilter,
                    Ponca::MeanNormal<DataPoint, _NFilter,
                        Ponca::MeanPosition<DataPoint, _NFilter,
                            Ponca::LocalFrame<DataPoint, _NFilter,
                                Ponca::Plane<DataPoint, _NFilter,
                                    Ponca::PrimitiveBase<DataPoint,_NFilter,T>>>>>>>>;
//! [ParabolicCylinderFit Definition]

#include "parabolicCylinderFit.hpp"
}
