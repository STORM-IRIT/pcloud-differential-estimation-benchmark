#pragma once

#include <Ponca/SpatialPartitioning>
#include <Eigen/Core>

/// Map a block to a Ponca point
template <typename _Scalar>
struct DataPoint {

  enum {Dim = 3};
  using Scalar     = _Scalar;
  using VectorType = Eigen::Matrix<Scalar, Dim, 1>;
  using MatrixType = Eigen::Matrix<Scalar, Dim, Dim>;

  PONCA_MULTIARCH inline DataPoint() : m_pos(VectorType::Zero()), m_nor(VectorType::Zero()) {}

  /// \brief Map a vector as ponca Point
  PONCA_MULTIARCH inline DataPoint(VectorType v, VectorType n) :
      m_pos(v), m_nor(n) {
  }

  PONCA_MULTIARCH inline VectorType pos() const { return m_pos; }
  PONCA_MULTIARCH inline VectorType normal() const { return m_nor; }

  // operator =
  PONCA_MULTIARCH inline DataPoint &operator=(const DataPoint &other) = default;

  // private:
  VectorType m_pos;
  VectorType m_nor;
};

template<typename Data, typename VectorType>
void
buildKdTree(const std::vector<VectorType>& cloudV, const std::vector<VectorType>& cloudN, std::unique_ptr<Ponca::KdTreeDense<Data>>& tree){
    std::vector<Data> bufs (cloudV.size());

#pragma omp parallel for
    for (int i = 0 ; i < cloudV.size(); ++i) {
        bufs[i] = Data(cloudV[i], cloudN[i]);
    }
  tree = std::make_unique<Ponca::KdTreeDense<Data>>(bufs);
}
