#pragma once

#include <vector>
#include <Eigen/Dense>
#include <DGtal/math/Statistic.h>


template <typename _Scalar>
struct Quantity {
    using Scalar = _Scalar;
    using VectorType = Eigen::Matrix<Scalar, 3, 1>;

public :
    Scalar k1 = 0; 
    Scalar k2 = 0; 
    Scalar mean = 0; 
    Scalar gauss = 0;
    Scalar neighbor_count = 0;
    VectorType d1 = VectorType(1, 0, 0);
    VectorType d2 = VectorType(1, 0, 0);
    VectorType normal = VectorType(1, 0, 0);
    VectorType projection = VectorType::Zero();

    int non_stable = 0;
    std::chrono::nanoseconds computing_time = std::chrono::nanoseconds(0);

    Quantity() = default; 

    Quantity(const Quantity& q) {
        this->k1 = q.k1;
        this->k2 = q.k2;
        this->mean = q.mean;
        this->gauss = q.gauss;
        this->d1 = q.d1;
        this->d2 = q.d2;

        this->neighbor_count = q.neighbor_count;

        this->projection = q.projection;
        this->normal = q.normal;

        this->non_stable = q.non_stable;
        this->computing_time = q.computing_time;
    }

    // To string
    std::string toString() const {
        std::stringstream ss;
        ss << "k1: " << k1 << ", k2: " << k2 << ", mean: " << mean << ", gauss: " << gauss
           << ", neighbor_count: " << neighbor_count
           << ", d1: [" << d1.transpose() << "], d2: [" << d2.transpose() << "]"
           << ", normal: [" << normal.transpose() << "]"
           << ", projection: [" << projection.transpose() << "]"
           << ", non_stable: " << non_stable
           << ", computing_time: " << computing_time.count() << " ns";
        return ss.str();
    }

};


/// Store the differential quantities
template <typename _Scalar>
struct DifferentialQuantities {
public:
    enum {Dim = 3};
    using Scalar     = _Scalar;
    using VectorSample = std::vector<Scalar>;
    using MatrixSample = std::vector<Eigen::Vector<Scalar, 3>>;

    DifferentialQuantities()  = default;
    ~DifferentialQuantities() = default;

    inline DifferentialQuantities(MatrixSample position, VectorSample k1, VectorSample k2, VectorSample mean, VectorSample gauss, MatrixSample d1, MatrixSample d2, MatrixSample normal, DGtal::Statistic<Scalar> statNei, DGtal::Statistic<Scalar> statTime) {
        m_k1 = k1;
        m_k2 = k2;
        m_mean = mean;
        m_gauss = gauss;
        m_d1 = d1;
        m_d2 = d2;
        m_normal = normal;
        m_statNeighbors = statNei;
        m_statTimings = statTime;
        m_position = position;
        m_nonstable = 0;
    }

    inline void setNonStableVector(std::vector<unsigned int> non_stables ) { 
        m_non_stable_ratio = non_stables; 
        m_nonstable = std::count_if(non_stables.begin(), non_stables.end(), [](Scalar x) { return x == 1; });
    }

    inline int getNumNonstable()  const { return m_nonstable; }

    inline Scalar getNonStableRatio () const { return Scalar( m_nonstable ) / Scalar( m_position.size() ); }

    inline std::vector<unsigned int> const getNonStableVector() { return m_non_stable_ratio;} 

    inline VectorSample k1()     const { return m_k1; }
    inline VectorSample k2()     const { return m_k2; }
    inline VectorSample mean()   const { return m_mean; }
    inline VectorSample gauss()  const { return m_gauss; }

    inline MatrixSample d1()     const { return m_d1; }
    inline MatrixSample d2()     const { return m_d2; }
    inline MatrixSample normal() const { return m_normal; }
    inline MatrixSample position() const { return m_position; }

    inline void setOriented(bool oriented) { m_oriented = oriented; }
    inline bool isOriented() const { return m_oriented; }

    inline DGtal::Statistic<Scalar> statNeighbors() const { return m_statNeighbors; }
    inline DGtal::Statistic<Scalar> statTimings()   const { return m_statTimings; }

private:
    VectorSample m_k1, m_k2, m_mean, m_gauss;
    MatrixSample m_d1, m_d2, m_normal, m_position;

    DGtal::Statistic<Scalar> m_statNeighbors;
    DGtal::Statistic<Scalar> m_statTimings;

    int m_nonstable = 0;
    std::vector<unsigned int> m_non_stable_ratio;

    bool m_oriented = true;
};
