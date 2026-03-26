template < class DataPoint, class _NFilter, typename T>
typename ParabolicCylinder<DataPoint, _NFilter, T>::Scalar
ParabolicCylinder<DataPoint, _NFilter, T>::potential( const VectorType &_q ) const
{
    VectorType x = Base::worldToLocalFrame(_q);
    return eval_quadratic_function(*(x.data() +1 ), *(x.data() + 2)) - ( *(x.data() ) );
}


template < class DataPoint, class _NFilter, typename T>
typename ParabolicCylinder<DataPoint, _NFilter, T>::VectorType
ParabolicCylinder<DataPoint, _NFilter, T>::project( const VectorType& _q ) const
{
    VectorType x = Base::worldToLocalFrame(_q);
    *(x.data()) = eval_quadratic_function(*(x.data()+1),*(x.data()+2));
    return Base::localFrameToWorld(x);
}

template < class DataPoint, class _NFilter, typename T>
typename ParabolicCylinder<DataPoint, _NFilter, T>::VectorType
ParabolicCylinder<DataPoint, _NFilter, T>::primitiveGradient( const VectorType& _q ) const
{
    VectorType proj = Base::worldToLocalFrame(_q);
    Vector2 df = evaldUV(proj(1), proj(2));
    VectorType local_gradient { 1, df(0) , df(1) };
    VectorType world_gradient = Base::template localFrameToWorld<true>(local_gradient);
    return world_gradient;
}

template < class DataPoint, class _NFilter, typename T>
typename ParabolicCylinder<DataPoint, _NFilter, T>::Scalar
ParabolicCylinder<DataPoint, _NFilter, T>::kMean() const {
  PONCA_MULTIARCH_STD_MATH(pow);
  static const Scalar one (1);
  static const Scalar two (2);
  static const Scalar threeOverTwo (Scalar(3)/Scalar(2));
  return ( dh_uu() * ( one + pow( dh_v(), two ) ) + dh_vv() * ( one + pow( dh_u(), two ) ) - two * dh_uv() * dh_u() * dh_v() ) /
    ( two * pow( one + pow( dh_u(), two ) + pow( dh_v(), two ), threeOverTwo ) );
}

template < class DataPoint, class _NFilter, typename T>
typename ParabolicCylinder<DataPoint, _NFilter, T>::Scalar
ParabolicCylinder<DataPoint, _NFilter, T>::GaussianCurvature() const {
    PONCA_MULTIARCH_STD_MATH(pow);
    static const Scalar one (1);
    static const Scalar two (2);
    return ( dh_uu()* dh_vv() - pow(dh_uv(),two)) /
        pow(( pow(dh_u(),two) + pow(dh_v(),two) + one ), two);
}

template < class DataPoint, class _NFilter, typename T>
typename ParabolicCylinder<DataPoint, _NFilter, T>::Scalar
ParabolicCylinder<DataPoint, _NFilter, T>::kmin() const {
    Scalar mean = kMean();
    Scalar gauss = GaussianCurvature();
    return mean - sqrt(mean*mean - gauss);
}

template < class DataPoint, class _NFilter, typename T>
typename ParabolicCylinder<DataPoint, _NFilter, T>::Scalar
ParabolicCylinder<DataPoint, _NFilter, T>::kmax() const {
    Scalar mean = kMean();
    Scalar gauss = GaussianCurvature();
    return mean + sqrt(mean*mean - gauss);
}

template < class DataPoint, class _NFilter, typename T>
typename ParabolicCylinder<DataPoint, _NFilter, T>::VectorType
ParabolicCylinder<DataPoint, _NFilter, T>::kminDirection() const {
    Eigen::Matrix<Scalar, 2, 2> I, II, W;  
    I << dE(), dF(), dF(), dG();
    II << dL(), dM(), dM(), dN();
    W = I.inverse() * II;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, 2, 2>> eigW(W);
    Eigen::Vector<Scalar, 2> dir = eigW.eigenvectors().col(0);

    VectorType v1 = VectorType(0, dir(0), dir(1));
    return Base::template localFrameToWorld<true>(v1);
}

template < class DataPoint, class _NFilter, typename T>
typename ParabolicCylinder<DataPoint, _NFilter, T>::VectorType
ParabolicCylinder<DataPoint, _NFilter, T>::kmaxDirection() const {

    Eigen::Matrix<Scalar, 2, 2> I, II, W;  
    I << dE(), dF(), dF(), dG();
    II << dL(), dM(), dM(), dN();
    W = I.inverse() * II;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, 2, 2>> eigW(W);
    Eigen::Vector<Scalar, 2> dir = eigW.eigenvectors().col(1);

    VectorType v2 = VectorType(0, dir(0), dir(1));
    return Base::template localFrameToWorld<true>(v2);
}
