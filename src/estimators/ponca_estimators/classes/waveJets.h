#pragma once


#include <Eigen/Core>
#include <Ponca/Fitting>
#include "localFrame.h"

namespace Ponca
{

template < class DataPoint, class _NFilter, typename T >
class WaveJets : public T
{
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    PONCA_FITTING_DECLARE_MATRIX_TYPE

protected:
    enum
    {
        check = Base::PROVIDES_PRIMITIVE_BASE
        && Base::PROVIDES_LOCAL_FRAME,
        PROVIDES_WAVEJETS
    };

public:
    PONCA_EXPLICIT_CAST_OPERATORS(WaveJets, waveJets)
    PONCA_FITTING_DECLARE_FINALIZE
    PONCA_FITTING_DECLARE_ADDNEIGHBOR

protected:

    using Vector2 = Eigen::Matrix<Scalar, 2, 1>;
    using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;

    Scalar     m_order; // Order of the polynomial, default is 2
    int        m_ncolM; // Number of columns of the matrix M

    Eigen::Matrix<std::complex<Scalar>, Eigen::Dynamic, Eigen::Dynamic> m_M; // Matrix of the system
    Eigen::Matrix<std::complex<Scalar>, Eigen::Dynamic, 1>              m_b; // Vector of the system
    int m_idx_j = 0; // Index of the current neighbor  when computing addNeighbors

    // Curvature info
    Scalar     m_k1;
    Scalar     m_k2;
    VectorType m_v1;
    VectorType m_v2;

    // TODO : Check if it's necessary, or change the LocalFrame already fitted.
    MatrixType m_P; // Corrected Local Frame 

    bool m_planeIsReady {false};
    bool m_normalIsCorrect {false};

public:
    PONCA_MULTIARCH inline void init()
    {
        Base::init();

        m_order = Scalar(2);
        m_planeIsReady = false;
        m_normalIsCorrect = false;
        m_ncolM = m_order*m_order/2 + 3*m_order/2+1;

        m_M.setZero();
        m_b.setZero();
        m_P.setZero();
        m_idx_j = 0;

        m_v1 = VectorType::Zero();
        m_v2 = VectorType::Zero();
        m_k1 = Scalar(0);
        m_k2 = Scalar(0);
    }

    PONCA_MULTIARCH inline VectorType project(const VectorType& _point) const
    {
        return _point; // This is done to avoid classical MLS steps, and avoid issue when using MLS methods. 
    }

    PONCA_MULTIARCH inline VectorType primitiveGradient () const
    {
        // Uniform gradient defined only by the orientation of the plane
        VectorType n = m_P.col(0);
        return n;
    }

    PONCA_MULTIARCH inline VectorType primitiveGradient (const VectorType&) const
    {
        return primitiveGradient();
    }

    PONCA_MULTIARCH inline Scalar kmin () const { return m_k1; }

    PONCA_MULTIARCH inline Scalar kmax () const { return m_k2; }

    PONCA_MULTIARCH inline Scalar kMean () const { return (m_k1 + m_k2) / Scalar(2); }

    PONCA_MULTIARCH inline Scalar GaussianCurvature () const { return m_k1 * m_k2; }

    PONCA_MULTIARCH inline VectorType kminDirection() const { return m_v1; }
    
    PONCA_MULTIARCH inline VectorType kmaxDirection() const { return m_v2; }

private:

    
    PONCA_MULTIARCH inline bool       m_local_jet_pass         (const VectorType& localPos);

    PONCA_MULTIARCH inline FIT_RESULT m_correct_normal_process ();

    PONCA_MULTIARCH inline FIT_RESULT m_jet_process            ();

}; //class WaveJets


/// \brief Helper alias for WaveJets fitting on points
//! [WaveJets Definition]
//! // wavejets using mean
template < class DataPoint, class _NFilter, typename T>
    using OrientedWaveJetsFit =
        Ponca::WaveJets<DataPoint, _NFilter,
                Ponca::MeanPlaneFitImpl<DataPoint, _NFilter,
                    Ponca::MeanNormal<DataPoint, _NFilter,
                        Ponca::MeanPosition<DataPoint, _NFilter,
                                Ponca::LocalFrame<DataPoint, _NFilter,
                                    Ponca::Plane<DataPoint, _NFilter,
                                        Ponca::PrimitiveBase<DataPoint,_NFilter,T>>>>>>>;
//  wavejets using COV
template < class DataPoint, class _NFilter, typename T>
    using WaveJetsFit =
    Ponca::WaveJets<DataPoint, _NFilter,
            Ponca::CovariancePlaneFitImpl<DataPoint, _NFilter,
                Ponca::CovarianceFitBase<DataPoint, _NFilter,
                        Ponca::MeanNormal<DataPoint, _NFilter,
                            Ponca::MeanPosition<DataPoint, _NFilter,
                                Ponca::LocalFrame<DataPoint, _NFilter,
                                    Ponca::Plane<DataPoint, _NFilter,
                                        Ponca::PrimitiveBase<DataPoint,_NFilter,T>>>>>>>>;


#include "waveJets.hpp"

}
