#include <Eigen/Geometry>

template < class DataPoint, class _NFilter, typename T>
void
ParabolicCylinderFitImpl<DataPoint, _NFilter, T>::init()
{
    Base::init();

    m_A_cov.setZero();
    m_F_cov.setZero();

    m_planeIsReady = false;
}

template < class DataPoint, class _NFilter, typename T>
bool
ParabolicCylinderFitImpl<DataPoint, _NFilter, T>::addLocalNeighbor(Scalar w,
                                                      const VectorType &localQ,
                                                      const DataPoint &attributes)
{
  auto res = Base::addLocalNeighbor(w, localQ, attributes);
  if(! m_planeIsReady)
  {
    return res;
  }
  else
  {
    VectorType local = Base::worldToLocalFrame(attributes.pos());

    const Scalar& f = *(local.data());
    const Scalar& x = *(local.data()+1);
    const Scalar& y = *(local.data()+2);

    Scalar xy = x*y;
    Scalar xx = x*x;
    Scalar yy = y*y;

    Eigen::Vector<Scalar, 7> v {1, x, y , xx, xy, xy, yy};

    m_A_cov += w * v * v.transpose();
    m_F_cov +=  w * v * f;

    return true;
  }
}



template < class DataPoint, class _NFilter, typename T>
FIT_RESULT
ParabolicCylinderFitImpl<DataPoint, _NFilter, T>::finalize () {
  PONCA_MULTIARCH_STD_MATH(abs);
  if (! m_planeIsReady) {
    FIT_RESULT res = Base::finalize();
    if (res == STABLE) {
      Base::computeFrameFromNormalVector(Base::plane().primitiveGradient());
      m_planeIsReady = true;
      return Base::m_eCurrentState = NEED_OTHER_PASS;
    }
    return res;
  }
  else {
    m_fitting_process();
    return Base::m_eCurrentState = STABLE;
  }
}

template < class DataPoint, class _NFilter, typename T>
void
ParabolicCylinderFitImpl<DataPoint, _NFilter, T>::m_fitting_process () {

  m_ellipsoid_fitting();
  m_parabolic_fitting();

}

template < class DataPoint, class _NFilter, typename T>
void
ParabolicCylinderFitImpl<DataPoint, _NFilter, T>::m_ellipsoid_fitting () {


  const Vector7 x = m_A_cov.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(m_F_cov);

  Base::m_uc      = x(0,0);
  Base::m_ul(0)   = x(1,0);
  Base::m_ul(1)   = x(2,0);
  Base::m_uq(0,0) = x(3,0);
  Base::m_uq(0,1) = x(4,0);
  Base::m_uq(1,0) = x(5,0);
  Base::m_uq(1,1) = x(6,0);
}

template < class DataPoint, class _NFilter, typename T>
void
ParabolicCylinderFitImpl<DataPoint, _NFilter, T>::m_parabolic_fitting () {
    PONCA_MULTIARCH_STD_MATH(abs);
    constexpr Scalar epsilon = Eigen::NumTraits<Scalar>::dummy_precision();

    Eigen::SelfAdjointEigenSolver<Matrix2> eig(Base::m_uq);
    Vector2 values = eig.eigenvalues();

    int higher = abs(values(0)) > abs(values(1)) ? 0 : 1;
    Scalar lambda0 = values(higher);
    Scalar lambda1 = values((higher + 1) % 2);

    Vector2 e_u = eig.eigenvectors().col(higher);
    Vector2 e_v = eig.eigenvectors().col((higher + 1) % 2);

    Scalar t = Base::getNeighborFilter().evalScale();
    Scalar alpha = Scalar(1);
    if (abs(lambda0 + Scalar(1)/t) > epsilon) {
        alpha = Scalar(2) * (abs(lambda0) - abs(lambda1)) / (abs(lambda0) + Scalar(1) / t);
    }
    alpha = (alpha < Scalar(1)) ? alpha : Scalar(1);
    alpha = (alpha > Scalar(0)) ? alpha : Scalar(0);

    Eigen::Matrix<Scalar, 7, 4> T_mat;
    T_mat.col(0) << 0, 0, 0, e_u(0)*e_u(0), e_u(0)*e_u(1), e_u(0)*e_u(1), e_u(1)*e_u(1);
    T_mat.col(1) << 0, e_u(0), e_u(1), 0, 0, 0, 0;
    T_mat.col(2) << 0, e_v(0), e_v(1), 0, 0, 0, 0;
    T_mat.col(3) << 1, 0, 0, 0, 0, 0, 0;

    Eigen::Matrix<Scalar, 4, 4> A_sub4 = T_mat.transpose() * m_A_cov * T_mat;
    Eigen::Matrix<Scalar, 4, 1> F_sub4 = T_mat.transpose() * m_F_cov;

    Eigen::Matrix<Scalar, 4, 1> P4 = A_sub4.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(F_sub4);

    Scalar a_raw = P4(0);

    Base::m_a = alpha * a_raw;

    Eigen::Matrix<Scalar, 7, 3> T_lin;
    T_lin << T_mat.col(1), T_mat.col(2), T_mat.col(3);

    Eigen::Matrix<Scalar, 7, 1> F_mod = m_F_cov - m_A_cov * (Base::m_a * T_mat.col(0));

    Eigen::Matrix<Scalar, 3, 3> A_sub3 = T_lin.transpose() * m_A_cov * T_lin;
    Eigen::Matrix<Scalar, 3, 1> F_sub3 = T_lin.transpose() * F_mod;

    Eigen::Matrix<Scalar, 3, 1> P3 = A_sub3.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(F_sub3);

    Base::m_uq = e_u * e_u.transpose();
    Base::m_ul = P3(0) * e_u + P3(1) * e_v;
    Base::m_uc = P3(2);

    VectorType v1 = VectorType(0, e_u(0), e_u(1));
    VectorType v2 = VectorType(0, e_v(0), e_v(1));
    Base::m_v1 = Base::template localFrameToWorld<true>(v1);
    Base::m_v2 = Base::template localFrameToWorld<true>(v2);
}
