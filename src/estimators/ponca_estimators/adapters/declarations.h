#pragma once

#include "types.h"
#include "DifferentialQuantities.hpp"
#include "fitDefinitions.h"
#include "poncaAdapters.hpp"


std::unique_ptr<KdTree> ponca_kdtree;
KnnGraph* ponca_knnGraph = nullptr;
bool use_kNNGraph = false;
Scalar radius = 0.1;
int mls_iter = 3;
int kNN = -1;
