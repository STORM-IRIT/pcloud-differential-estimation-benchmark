#pragma once

#include <Eigen/Dense>
#include <Ponca/Common>
#include <Ponca/Fitting>
#include "../../IO/PointCloudDiff.h"

#include "adapters/EstimatorFactory.h"
#include "adapters/declarations.h"

#include <omp.h>

/// Convenience function measuring and printing the processing time of F
template <typename Functor>
void measureTime( const std::string &actionName, Functor F ){
    using namespace std::literals; // enables the usage of 24h instead of e.g. std::chrono::hours(24)

    const std::chrono::time_point<std::chrono::steady_clock> start =
            std::chrono::steady_clock::now();
    F(); // run process
    const auto end = std::chrono::steady_clock::now();
}

/// Generic processing function: traverse point cloud, compute fitting, and use functor to process fitting output
/// \note Functor is called only if fit is stable
template<typename Estimator, typename Functor>
void processPointCloud( Estimator& estimator, int uniqueIdx, Functor f ){
    int nvert = ponca_kdtree->points().size();

    if ( uniqueIdx > -1 ){
        Quantity<Scalar> q;
        estimator( uniqueIdx, q );
        f ( uniqueIdx, q );
        return;
    }

#ifdef NDEBUG
#pragma omp parallel for
#endif
    for (int i = 0; i < nvert; ++i) {
#ifdef DEBUG
        DGtal::trace.progressBar( i, nvert );
#endif
        Quantity<Scalar> q;
        estimator(i, q);
        f (i, q);
    }
}

/// Generic processing function: traverse point cloud and compute mean, first and second curvatures + their direction
/// \tparam FitT Defines the type of estimator used for computation
template<typename Estimator>
DifferentialQuantities<Scalar> estimateDifferentialQuantities_impl( Estimator &estimator, int uniqueIdx ) {

    int nvert = ponca_kdtree->points().size();

    VectorType defaultVectorType = VectorType(1, 0, 0); // Keep the norm to 1 to avoid numerical issues for non-stable points

    DGtal::Statistic<Scalar> statTime, statNei;
    std::vector<Scalar> mean ( nvert ), kmin ( nvert ), kmax ( nvert ), gauss ( nvert );
    std::vector<Eigen::Vector<Scalar, 3>> proj ( nvert ), normal( nvert, defaultVectorType ), dmin( nvert, defaultVectorType ), dmax( nvert,defaultVectorType );
    std::vector <unsigned int> non_stable_idx(nvert, 0);

    measureTime( "[Ponca] Compute differential quantities using " + estimator.getName(),
                 [&uniqueIdx, &estimator, &proj, &mean, &kmin, &kmax, &gauss, &normal, &dmin, &dmax, &statTime, &statNei, &non_stable_idx]() {
                            processPointCloud(estimator, uniqueIdx, 
                                [&proj, &mean, &kmin, &kmax, &gauss, &normal, &dmin, &dmax, &statTime, &statNei, &non_stable_idx]
                                (
                                    int i, Quantity<Scalar> &quantity
                                ) {
                                    statNei.addValue( quantity.neighbor_count );
                                    // statNei.addValue( fit.getNumNeighbors() );
                                    statTime.addValue( quantity.computing_time.count() );

                                    if (quantity.non_stable >= 1) {
                                        non_stable_idx[i] = 1;
                                        return;
                                    }
                                    mean  [ i ]   = quantity.mean;
                                    kmax  [ i ]   = quantity.k2;
                                    kmin  [ i ]   = quantity.k1;
                                    dmax  [ i ]   = quantity.d2;
                                    dmin  [ i ]   = quantity.d1;
                                    gauss [ i ]   = quantity.gauss;
                                    normal[ i ]   = quantity.normal;
                                    proj  [ i ]   = quantity.projection;
                                });
                            });

    // if there is unstable points, we don't return the results

    DifferentialQuantities<Scalar> dq = { proj, kmin, kmax, mean, gauss, dmin, dmax, normal, statNei, statTime };
    dq.setNonStableVector(non_stable_idx);
    return dq;
}

template<typename WeightFunc>
DifferentialQuantities<Scalar> estimateDifferentialQuantities( std::string name, int uniqueIdx = -1 ) {

    std::shared_ptr<EstimatorFactory< WeightFunc >> factory = getEstimatorFactory<WeightFunc>();

    auto estimator = factory->getEstimator(name);

    DifferentialQuantities<Scalar> diff = estimateDifferentialQuantities_impl( *estimator, uniqueIdx ); 
    diff.setOriented(estimator->isOriented());
    return diff;
}

template <typename WeightFunc>
std::vector<Scalar> neiRequest(int index) {
    int nvert = ponca_kdtree->points().size();
    std::vector <Scalar> nei_value(nvert, 0);

    VectorType pos = ponca_kdtree->points()[index].pos();
    WeightFunc w(pos, radius);

    processNeighbors(pos, index, [&w, &pos, &nei_value, &index](int j) {
        // fit.addNeighbor(ponca_kdtree.point_data()[j]);
        if (kNN != -1){
            nei_value[j] = ( j != index ) ? 1 : 2;
        }
        else {
            const auto& q = ponca_kdtree->points()[j];
            nei_value[j] = w(q).first;
        }
    });
    
    return nei_value;
}

Scalar estimateRadiusFromKNN(int k) {
    Scalar radius = 0;
    int nvert = ponca_kdtree->points().size();
    for (int i = 0; i < nvert; ++i) {
        VectorType pos = ponca_kdtree->points()[i].pos();
        Scalar max_dist = 0;
        
        for (int j : ponca_kdtree->kNearestNeighbors(i, k)){
            max_dist = std::max(max_dist, (pos - ponca_kdtree->points()[j].pos()).norm());
        }

        radius += max_dist;
    }
    radius /= nvert;
    return radius;
}
