#include <Eigen/SVD>
#include <Eigen/Geometry>

template < class DataPoint, class _NFilter, typename T>
void
NormalCovariance2D<DataPoint, _NFilter, T>::init()
{
    Base::init();

    m_planeIsReady = false;
    m_cov2D.setZero();
    m_normal_centroid2D.setZero();
    m_P.setZero();
}

template < class DataPoint, class _NFilter, typename T>
bool
NormalCovariance2D<DataPoint, _NFilter, T>::addLocalNeighbor(Scalar w,
                                                      const VectorType &localQ,
                                                      const DataPoint &attributes)
{
    auto res = Base::addLocalNeighbor(w, localQ, attributes);
    if(! m_planeIsReady)
    {
        return res;
    }
    else // base plane is ready, we can now fit the patch
    {
        const Vector2 n_proj = m_P.transpose() * (attributes.normal());
        m_cov2D += w * n_proj * n_proj.transpose();
        m_normal_centroid2D += w * n_proj;
        return true;
    }
    return false;
}

template < class DataPoint, class _NFilter, typename T>
FIT_RESULT
NormalCovariance2D<DataPoint, _NFilter, T>::finalize ()
{
    if (! m_planeIsReady) {
        FIT_RESULT res = Base::finalize();

        if(res == STABLE) {  // plane is ready
            Base::computeFrameFromNormalVector(Base::plane().primitiveGradient());
            setTangentPlane();
            m_planeIsReady = true;
            return Base::m_eCurrentState = NEED_OTHER_PASS;
        }
        return res;
    }
    else {
        m_normal_centroid2D /= Base::getWeightSum();
        m_cov2D = ( m_cov2D / Base::getWeightSum() ) - ( m_normal_centroid2D * m_normal_centroid2D.transpose() );
        m_solver.compute(m_cov2D);

        if (m_solver.info() != Eigen::Success) {
            return Base::m_eCurrentState = UNDEFINED;
        }
        return Base::m_eCurrentState = STABLE;
    }
}

template < class DataPoint, class _NFilter, typename T>
typename NormalCovariance2D<DataPoint, _NFilter, T>::Scalar
NormalCovariance2D<DataPoint, _NFilter, T>::kMean() const {
  return ( kmin() + kmax() ) / Scalar(2);
}

template < class DataPoint, class _NFilter, typename T>
typename NormalCovariance2D<DataPoint, _NFilter, T>::Scalar
NormalCovariance2D<DataPoint, _NFilter, T>::GaussianCurvature() const {
    return kmin() * kmax();
}

template < class DataPoint, class _NFilter, typename T>
typename NormalCovariance2D<DataPoint, _NFilter, T>::Scalar
NormalCovariance2D<DataPoint, _NFilter, T>::kmin() const {
    PONCA_MULTIARCH_STD_MATH(pow);
    PONCA_MULTIARCH_STD_MATH(sqrt);
    constexpr Scalar four = Scalar(4);
    constexpr Scalar PI = Scalar(M_PI);
    Scalar factor = pow ( Base::getNeighborFilter().evalScale(), four ) * PI / four;
    return sqrt( m_solver.eigenvalues()(0) / factor );
    // return m_solver.eigenvalues()(0);
}

template < class DataPoint, class _NFilter, typename T>
typename NormalCovariance2D<DataPoint, _NFilter, T>::Scalar
NormalCovariance2D<DataPoint, _NFilter, T>::kmax() const {
    PONCA_MULTIARCH_STD_MATH(pow);
    PONCA_MULTIARCH_STD_MATH(sqrt);
    constexpr Scalar four = Scalar(4);
    constexpr Scalar PI = Scalar(M_PI);
    Scalar factor = pow ( Base::getNeighborFilter().evalScale(), four ) * PI / four;
    return sqrt( m_solver.eigenvalues()(1) / factor );
    // return m_solver.eigenvalues()(1);
}

template < class DataPoint, class _NFilter, typename T>
typename NormalCovariance2D<DataPoint, _NFilter, T>::VectorType
NormalCovariance2D<DataPoint, _NFilter, T>::kminDirection() const {
    Vector2 dir = m_solver.eigenvectors().col(0);
    VectorType vmin = VectorType(0, dir(0), dir(1));
    return Base::template localFrameToWorld<true>(vmin);
}

template < class DataPoint, class _NFilter, typename T>
typename NormalCovariance2D<DataPoint, _NFilter, T>::VectorType
NormalCovariance2D<DataPoint, _NFilter, T>::kmaxDirection() const {
    Vector2 dir = m_solver.eigenvectors().col(1);
    VectorType vmax = VectorType(0, dir(0), dir(1));
    return Base::template localFrameToWorld<true>(vmax);
}
