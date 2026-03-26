#pragma once

#include "types.h"
#include "../classes/cnc_separate.h"
#include "../classes/sphereCurvature.h"
#include "../classes/covariance2DFit.h"
#include "../classes/planeCurvature.h"
#include "../classes/meanCurvature.h"
#include "../classes/normalCovariance2DFit.h"
#include "../classes/normalCovariance3DFit.h"
#include "../classes/normalW.h"
#include "../classes/shapeOperator2D.h"
#include "../classes/parabolicCylinderFit.h"
#include "../classes/varifolds.h"
#include "../classes/waveJets.h"
#include "../classes/quadricFit.h"

template <typename WeightFunc>
using fit_ASO = Ponca::BasketDiff<
            Ponca::Basket<PPAdapter, WeightFunc, Ponca::OrientedSphereFit>,
            Ponca::FitSpaceDer,
            Ponca::OrientedSphereDer, Ponca::MlsSphereFitDer,
            Ponca::CurvatureEstimatorDer, Ponca::NormalDerivativeWeingartenEstimator,
            Ponca::WeingartenCurvatureEstimatorDer>;

template <typename WeightFunc>
using fit_PCA = Ponca::Basket<PPAdapter, WeightFunc, Ponca::CovPlaneCurvatureFit>;

template <typename WeightFunc>
using fit_MeanPLANE = Ponca::Basket<PPAdapter, WeightFunc, Ponca::MeanPlaneCurvatureFit>;

template <typename WeightFunc>
using fit_APSS = Ponca::Basket<PPAdapter, WeightFunc, Ponca::SimpleOrientedSphereFit>;

template <typename WeightFunc>
using fit_Sphere = Ponca::Basket<PPAdapter, WeightFunc, Ponca::SimpleSphereFit>;

template <typename WeightFunc>
using fit_UnorientedSphere = Ponca::Basket<PPAdapter, WeightFunc, Ponca::SimpleUnorientedSphereFit>;

template <typename WeightFunc>
using fit_MeanCurvature = Ponca::Basket<PPAdapter, WeightFunc, Ponca::MeanCurvatureFit>;

template <typename WeightFunc>
using fit_Cov2D = Ponca::Basket<PPAdapter, WeightFunc, Ponca::Covariance2DFit>;

template <typename WeightFunc>
using fit_NormCov2D = Ponca::Basket<PPAdapter, WeightFunc, Ponca::NormalCovariance2DFit>;

template <typename WeightFunc>
using fit_NormCov3D = Ponca::Basket<PPAdapter, WeightFunc, Ponca::NormalCovariance3DFit>;

template<typename WeightFunc>
using fit_ShapeOperator = Ponca::Basket<PPAdapter, WeightFunc, Ponca::ShapeOperator2DFit>;

// using fit_CNC_Uniform = Ponca::CNC<PPAdapter, Ponca::UniformGeneration>;
// using fit_CNC_Independent = Ponca::CNC<PPAdapter, Ponca::IndependentGeneration>;
// using fit_CNC_AvgHexagram = Ponca::CNC<PPAdapter, Ponca::AvgHexagramGeneration>;
// using fit_CNC_Hexagram = Ponca::CNC<PPAdapter, Ponca::HexagramGeneration>;

using fit_CNC_Uniform = Ponca::CNC_separate<PPAdapter, Ponca::UniformGeneration>;
using fit_CNC_Independent = Ponca::CNC_separate<PPAdapter, Ponca::IndependentGeneration>;
using fit_CNC_AvgHexagram = Ponca::CNC_separate<PPAdapter, Ponca::AvgHexagramGeneration>;
using fit_CNC_Hexagram = Ponca::CNC_separate<PPAdapter, Ponca::HexagramGeneration>;

template <typename WeightFunc>
using fit_WaveJets = Ponca::Basket<PPAdapter, WeightFunc, Ponca::WaveJetsFit>;

template <typename WeightFunc>
using fit_OrientedWaveJets = Ponca::Basket<PPAdapter, WeightFunc, Ponca::OrientedWaveJetsFit>;

using fit_Varifolds = Ponca::Basket<PPAdapter, VarifoldWeightFunc, Ponca::VarifoldsCovPlane>;

using fit_VarifoldsMeanPlane = Ponca::Basket<PPAdapter, VarifoldWeightFunc, Ponca::VarifoldsMeanPlane>;

template <typename WeightFunc>
using fit_Monge = Ponca::Basket<PPAdapter, WeightFunc, Ponca::MongePatchQuadraticFit>;

template <typename WeightFunc>
using fit_PCMLS = Ponca::Basket<PPAdapter, WeightFunc, Ponca::BaseParabolicCylinderFit>;

template <typename WeightFunc>
using fit_quadric = Ponca::Basket<PPAdapter, WeightFunc, Ponca::QuadricFit>;
