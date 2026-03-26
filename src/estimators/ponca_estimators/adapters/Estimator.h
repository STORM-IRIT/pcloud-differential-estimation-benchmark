#pragma once

#include "fitDefinitions.h"
#include "types.h"
#include <chrono>

#include "declarations.h"

/// Function to compute the neighbors of a point
template <typename Functor>
void processNeighbors(const int &idx, const Functor f){
    // VectorType pos = ponca_kdtree.point_data()[idx].pos();
    if(use_kNNGraph){
        if (kNN != -1){
            for (int j : ponca_knnGraph->kNearestNeighbors(idx)){
                f(j);
            }
            f(idx);
        }
        else {
            for (int j : ponca_knnGraph->rangeNeighbors(idx, radius)){
                f(j);
            }
            f(idx);
        }
    } else {
        if (kNN != -1){
            for (int j : ponca_kdtree->kNearestNeighbors(idx, kNN)){
                f(j);
            }
            f(idx);
        }
        else {
            for (int j : ponca_kdtree->rangeNeighbors(idx, radius)){
                f(j);
            }
            f(idx);
        }
    }
}

/// Function to compute the neighbors of a point
template <typename Functor>
void processNeighbors(const VectorType &pos, const int &init_idx, const Functor f){
    if(use_kNNGraph){
        if (kNN != -1){
            for (int j : ponca_knnGraph->kNearestNeighbors(init_idx)){
                f(j);
            }
            f(init_idx);
        }
        else {
            for (int j : ponca_knnGraph->rangeNeighbors(init_idx, radius)){
                f(j);
            }
            f(init_idx);
        }
    } else {
        if (kNN != -1){
            for (int j : ponca_kdtree->kNearestNeighbors(pos, kNN)){
                f(j);
            }
        }
        else {
            for (int j : ponca_kdtree->rangeNeighbors(pos, radius)){
                f(j);
            }
                
        }
    }
}

template < typename _Fit > 
class InterfaceProcessHandler {

    protected : 
        using FitT = _Fit;
        using WeightFunc = typename FitT::NeighborFilter;
    
    public : 

        ~InterfaceProcessHandler() = default;
        InterfaceProcessHandler() {
            init_idx = 0;
            // current_pos = VectorType::Zero();
            current_fitResult = Ponca::UNDEFINED;
            current_NeighborsPos.clear();
            current_NeighborsNormal.clear();
            current_NeighborsCount = 0;
            current_NeighborsTiming = std::chrono::nanoseconds::zero();
            current_ComputeTiming = std::chrono::nanoseconds::zero();
        };

        void initHandler(int query_idx) {
          init_idx = query_idx;
          current_data = ponca_kdtree->points()[ query_idx ];
        }

        void initEstimator( FitT &fit ) {
            fit.setNeighborFilter({current_data, radius});
            fit.init();
            current_NeighborsCount = 0;
            current_NeighborsPos.clear();
            current_NeighborsNormal.clear();
        }

        virtual void applyNeighbors(FitT &fit) = 0;

        Ponca::FIT_RESULT finalize ( FitT &fit ) {
            
            auto tCompute_start = std::chrono::high_resolution_clock::now();

            current_fitResult = fit.finalize();
            if (current_fitResult == Ponca::STABLE) {
              current_data.m_pos = fit.project(current_data.m_pos);
              current_data.m_nor = fit.primitiveGradient();
            }

            auto tCompute_end = std::chrono::high_resolution_clock::now();
            current_ComputeTiming += std::chrono::duration_cast<std::chrono::nanoseconds>(tCompute_end - tCompute_start);
            return current_fitResult;
        }

        void applyFunctor( FitT &fit, Quantity<Scalar> &diffQuantity ) {
            diffQuantity.neighbor_count = current_NeighborsCount;
            diffQuantity.computing_time = current_ComputeTiming;
            diffQuantity.non_stable = (current_fitResult != Ponca::STABLE || current_NeighborsCount < 3);

            if (diffQuantity.non_stable) return;

            diffQuantity.k1 = fit.kmin();
            diffQuantity.k2 = fit.kmax();
            diffQuantity.mean = fit.kMean();
            diffQuantity.gauss = fit.GaussianCurvature();
            diffQuantity.d1 = fit.kminDirection();
            diffQuantity.d2 = fit.kmaxDirection();
            diffQuantity.normal = fit.primitiveGradient();
            diffQuantity.projection = current_data.m_pos;
        }

        PPAdapter                current_data = PPAdapter(VectorType::Zero(), VectorType::Zero());

    protected :

        int                      init_idx = 0;
        Ponca::FIT_RESULT        current_fitResult = Ponca::UNDEFINED;
        std::vector<VectorType>  current_NeighborsPos;
        std::vector<VectorType>  current_NeighborsNormal;
        int                      current_NeighborsCount = 0;

        std::chrono::nanoseconds current_NeighborsTiming = std::chrono::nanoseconds::zero();
        std::chrono::nanoseconds current_ComputeTiming = std::chrono::nanoseconds::zero();

}; // InterfaceProcessHandler

// Generic ProcessHandler
template < typename _Fit > 
class ProcessHandler : public InterfaceProcessHandler< _Fit > {

    using Base = InterfaceProcessHandler< _Fit >;
    using FitT = typename Base::FitT;

    public :
        std::string name;

    public : 

        ProcessHandler (const std::string &name) : name(name), Base() {}
        ~ProcessHandler() = default;

        void applyNeighbors( FitT & fit ) override {
            Base::current_NeighborsCount = 0;
            fit.startNewPass();
            auto tNei_start = std::chrono::high_resolution_clock::now();

            processNeighbors(Base::current_data.m_pos, Base::init_idx, [&]( int idx ){
                fit.addNeighbor( ponca_kdtree->points()[idx] );
                Base::current_NeighborsCount ++;
            });

            auto tNei_end = std::chrono::high_resolution_clock::now();
            Base::current_NeighborsTiming += std::chrono::duration_cast<std::chrono::nanoseconds>( tNei_end - tNei_start );
        }

}; // ProcessHandler

template <typename T>
concept IsCNCType = std::is_same_v<T, fit_CNC_Uniform> ||
                    std::is_same_v<T, fit_CNC_AvgHexagram> ||
                    std::is_same_v<T, fit_CNC_Hexagram> ||
                    std::is_same_v<T, fit_CNC_Independent>;

// Specific ProcessHandler for CNC
template <IsCNCType FitT>
class ProcessHandler< FitT > : public InterfaceProcessHandler< FitT > {

    using Base = InterfaceProcessHandler< FitT >;

    public :
        std::string name;
    
    public : 

        ProcessHandler (const std::string &name) : name( name ), Base() {}
        ~ProcessHandler() = default;

        void initEstimator( FitT &fit ) {
          Base::initEstimator(fit);
        }

        void applyNeighbors(FitT &fit) override {
          Base::current_NeighborsCount = 0;
          // fit.startNewPass();

          std::vector<PPAdapter> neighborData;
          processNeighbors(Base::current_data.m_pos, Base::init_idx, [&]( int idx ){
                neighborData.push_back(ponca_kdtree->points()[idx]);
                Base::current_NeighborsCount ++;
            });

          const auto t_start = std::chrono::high_resolution_clock::now();

          fit.constructTriangles(neighborData);

          const auto t_end = std::chrono::high_resolution_clock::now();
          Base::current_NeighborsTiming += std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start);
        }

        Ponca::FIT_RESULT finalize(FitT &fit) {
          auto tCompute_start = std::chrono::high_resolution_clock::now();

          Base::current_fitResult = fit.finalize();

          const auto tCompute_end = std::chrono::high_resolution_clock::now();
          Base::current_ComputeTiming += std::chrono::duration_cast<std::chrono::nanoseconds>(tCompute_end - tCompute_start);

          return Base::current_fitResult;
        }

  void applyFunctor( FitT &fit, Quantity<Scalar> &diffQuantity ) {
          diffQuantity.neighbor_count = Base::current_NeighborsCount;
          diffQuantity.computing_time = Base::current_ComputeTiming;
          diffQuantity.non_stable = (Base::current_fitResult != Ponca::STABLE || Base::current_NeighborsCount < 3);

          if (diffQuantity.non_stable) return;

          diffQuantity.k1 = fit.kmin();
          diffQuantity.k2 = fit.kmax();
          diffQuantity.mean = fit.kMean();
          diffQuantity.gauss = fit.GaussianCurvature();
          diffQuantity.d1 = fit.kminDirection();
          diffQuantity.d2 = fit.kmaxDirection();
          diffQuantity.normal = Base::current_data.m_nor; // <- It is not a normal estimator
          diffQuantity.projection = Base::current_data.m_pos;
        }

        
}; // ProcessHandler for CNC


class BaseEstimator {
    
    public :

        virtual ~BaseEstimator() = default;
        virtual std::string getName () const = 0;
        virtual std::string toString() const = 0;
        virtual bool isOriented() const = 0;
        
        virtual void operator() ( const int &query_idx, Quantity<Scalar> &quantity ) {
            std::cout << " Only here to check. " << std::endl;
        }
}; // BaseEstimator


// The Handler store the functors, called in the main loop
template < typename _FitT >
class Estimator : public BaseEstimator {
    using Self = Estimator<_FitT>;

    protected :
        bool oriented;
        std::string name = "None";
        int mls_max = 1;

    public : 
        // Overload the FitT
        using FitT = _FitT;

        using WeightFunc = typename FitT::NeighborFilter;

        Estimator() = default;
        Estimator(std::string name, bool oriented, int mls_max) : name(name), oriented(oriented), mls_max(mls_max) {}

        std::string getName () const override {
            return name;
        }

        std::string toString() const override {
            return name;
        }

        bool isOriented() const override {
            return oriented;
        }

        void operator() ( const int &query_idx, Quantity<Scalar> &quantity ) override {
            
            ProcessHandler<FitT> pHandler( name );

            int mls_current = 0; 

            pHandler.initHandler ( query_idx );

            for ( mls_current = 0 ; mls_current < mls_max ; mls_current ++ ){
                FitT fit;
                pHandler.initEstimator( fit );
                Ponca::FIT_RESULT res;
                do {
                    pHandler.applyNeighbors( fit );                
                    res = pHandler.finalize( fit );
                } while ( res == Ponca::NEED_OTHER_PASS );

                if ( mls_current == mls_max - 1 ){
                    pHandler.applyFunctor( fit, quantity );
                }
            }
        }
}; // Estimator
