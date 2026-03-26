#pragma once

#include <Eigen/Dense>
#include <Ponca/Common>
#include <Ponca/Fitting>

#include "poncaAdapters.hpp"

#include "../classes/varifolds.h"

// Types definition
// using Scalar             = double; => Defined in the PointCloudDiff.h
using VectorType         = Eigen::Matrix<Scalar, 3,1>;
using PPAdapter          = DataPoint<Scalar>;
using KnnGraph           = Ponca::KnnGraph<PPAdapter>;
using KdTree             = Ponca::KdTreeDense<PPAdapter>;

// Weighting functions
using SmoothWeightFunc      = Ponca::DistWeightFunc<PPAdapter, Ponca::SmoothWeightKernel<Scalar> >;
using ConstWeightFunc       = Ponca::DistWeightFunc<PPAdapter, Ponca::ConstantWeightKernel<Scalar> >;
using WendlandWeightFunc    = Ponca::DistWeightFunc<PPAdapter, Ponca::WendlandWeightKernel<Scalar> >;
using SingularWeightFunc    = Ponca::DistWeightFunc<PPAdapter, Ponca::SingularWeightKernel<Scalar> >;

using VarifoldWeightFunc    = Ponca::DistWeightFunc<PPAdapter, Ponca::VarifoldWeightKernel<Scalar> >;
