/**
Copyright (c) 2022
 Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr)
 Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France,

All rights reserved.

*/

#pragma once

#include <Ponca/src/Common/defines.h>
#include <Ponca/src/Fitting/cnc.h>
#include <Ponca/src/Fitting/cncFormulaEigen.h>
#include <Ponca/src/Fitting/weightFunc.h>
#include <Ponca/src/Fitting/weightKernel.h>

namespace Ponca
{

template < class P, TriangleGenerationMethod _method = UniformGeneration>
class CNC_separate : public CNC<P, _method> {
protected:
    enum
    {
        PROVIDES_PRINCIPAL_CURVATURES
    };
public:
    using DataPoint = P;
    using MatrixType = typename DataPoint::MatrixType;
    using Scalar     = typename DataPoint::Scalar;
    using VectorType = typename DataPoint::VectorType;
    typedef Eigen::VectorXd  DenseVector;
    typedef Eigen::MatrixXd  DenseMatrix;
    using NeighborFilter = NeighborFilterStoreNormal<DataPoint, NoWeightFunc<DataPoint>>;
    using BaseCNC = CNC<P, _method>;

public:

    template <typename PointContainer>
    PONCA_MULTIARCH inline void constructTriangles(const PointContainer& points)
    {
      BaseCNC::init();
      std::vector<unsigned int> indicesSample(points.size());
      std::iota(indicesSample.begin(), indicesSample.end(), 0);

      this->m_eCurrentState = internal::TriangleGenerator<_method, P>::generate( indicesSample, points, this->m_nFilter, this->m_triangles);
      if (this->m_eCurrentState != STABLE) return;
      this->m_nb_vt = this->m_triangles.size();
    }

    PONCA_MULTIARCH inline FIT_RESULT finalize()
    {
      if (this->m_eCurrentState != STABLE) return this->m_eCurrentState;
      return BaseCNC::finalize();
    }

    //! \brief Comparison operator
    PONCA_MULTIARCH [[nodiscard]] bool operator==(const CNC_separate& other) const {
        // We use the matrix to compare the fitting results
        return (this->m_eCurrentState == other.m_eCurrentState)
            && (BaseCNC::kMean() == other.kMean())
            && (BaseCNC::kmin() == other.kmin())
            && (BaseCNC::kmax() == other.kmax())
            && (BaseCNC::kminDirection() == other.kminDirection())
            && (BaseCNC::kmaxDirection() == other.kmaxDirection())
            && (BaseCNC::GaussianCurvature() == other.GaussianCurvature())
            && (BaseCNC::m_T11 == other.m_T11) && (BaseCNC::m_T12 == other.m_T12) && (BaseCNC::m_T13 == other.m_T13)
            && (BaseCNC::m_T22 == other.m_T22) && (BaseCNC::m_T23 == other.m_T23)
            && (BaseCNC::m_T33 == other.m_T33);
    }

    //! \brief Comparison operator, convenience function
    PONCA_MULTIARCH [[nodiscard]] bool operator!=(const CNC_separate& other) const {
        // We use the matrix to compare the fitting results
        return !(this == &other);
    }

}; //class CNC_separate

} // namespace Ponca
