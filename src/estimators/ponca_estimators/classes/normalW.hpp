#include <Eigen/SVD>
#include <Eigen/Geometry>

template < class DataPoint, class _NFilter, typename T>
void
NormalWeingarten<DataPoint, _NFilter, T>::init()
{
    Base::init();

    m_XtX = Matrix3::Zero();
    m_XtY1 = Vector3::Zero();
    m_XtY2 = Vector3::Zero();

    m_mu1 = Scalar(0);
    m_mu2 = Scalar(0);
    m_w_tilde = Vector4::Zero();
    m_normal = VectorType::Zero();
    m_normal_no_proj = VectorType::Zero();

    m_dmin = VectorType::Zero();
    m_dmax = VectorType::Zero();
    m_kmin = Scalar(0);
    m_kmax = Scalar(0);

    m_planeIsReady = false;
}

template < class DataPoint, class _NFilter, typename T>
bool
NormalWeingarten<DataPoint, _NFilter, T>::addLocalNeighbor(Scalar w,
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
        const Vector2 dP_proj = m_P.transpose() * localQ;
        const Vector2 dN_proj = m_P.transpose() * attributes.normal();

        Vector3 Xi = Vector3(Scalar(1), dP_proj(0), dP_proj(1));

        m_XtX += w * Xi * Xi.transpose();
        m_XtY1 += w * Xi * dN_proj(0);
        m_XtY2 += w * Xi * dN_proj(1);
    
        return true;
    }
    return false;
}

template < class DataPoint, class _NFilter, typename T>
FIT_RESULT
NormalWeingarten<DataPoint, _NFilter, T>::finalize ()
{
    if (! m_planeIsReady) {
        FIT_RESULT res = Base::finalize();

        if(res == STABLE) {  // plane is ready
            setTangentPlane();
            m_planeIsReady = true;
            return Base::m_eCurrentState = NEED_OTHER_PASS;
        }
        return res;
    }
    else {

        if( Base::getNumNeighbors() < 9 )
            return Base::m_eCurrentState = UNDEFINED;
        
        Eigen::JacobiSVD<Matrix3> svd(m_XtX, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Vector3 beta1 = svd.solve(m_XtY1);
        Vector3 beta2 = svd.solve(m_XtY2);

        m_mu1          = beta1(0);
        m_mu2          = beta2(0);


        if (m_normal == VectorType::Zero()) {
            m_normal = VectorType( sqrt(1 - pow( m_mu1, Scalar(2) ) - pow( m_mu2, Scalar(2) ) ), m_mu1, m_mu2 );
            m_normal_no_proj = VectorType( m_mu1, m_mu2, sqrt(1 - pow( m_mu1, Scalar(2) ) - pow( m_mu2, Scalar(2) ) ) );
            m_normal = Base::template localFrameToWorld<true>(m_normal);
            m_normal.normalize();
            Base::computeFrameFromNormalVector(m_normal);
            return Base::m_eCurrentState = NEED_OTHER_PASS;
        }

        Vector4 w_bar       = Vector4( beta1(1), beta1(2), beta2(1), beta2(2) );
        Vector4 planeNormal = Vector4( m_mu1*m_mu2, - ( 1 - m_mu2*m_mu2 ), (1 - m_mu1*m_mu1 ), - m_mu1*m_mu2  );
        Scalar  planeOffset = Scalar(0);

        m_w_tilde             = projectOnto4DPlane(w_bar, planeNormal, planeOffset);

        std::cout << "[NEED TO BE 0] Value to check : " << m_w_tilde.dot(planeNormal) << std::endl;

        computeCurvature();

        return Base::m_eCurrentState = STABLE;
    }
}

template < class DataPoint, class _NFilter, typename T>
void
NormalWeingarten<DataPoint, _NFilter, T>::computeCurvature() {

    Vector2 evalPos2D = m_P.transpose() * Base::m_w.basisCenter();

    Scalar Xu = - ( ( 1.0 + f_v() * f_v() ) * f_uu() - f_u() * f_v() * f_vu() ) * pow( m_normal_no_proj(2) , 3);
    Scalar Xv = - ( ( 1.0 + f_v() * f_v() ) * f_uv() - f_u() * f_v() * f_vv() ) * pow( m_normal_no_proj(2) , 3);
    Scalar Yu = - ( ( 1.0 + f_u() * f_u() ) * f_vu() - f_u() * f_v() * f_uu() ) * pow( m_normal_no_proj(2) , 3);
    Scalar Yv = - ( ( 1.0 + f_u() * f_u() ) * f_vv() - f_u() * f_v() * f_uv() ) * pow( m_normal_no_proj(2) , 3);
    
    Matrix2 W;
    // W << m_w_tilde(0), m_w_tilde(1), 
    //      m_w_tilde(2), m_w_tilde(3);
    W << Xu, Xv, 
         Yu, Yv;
    

    std::cout << "Xu : " << Xu << std::endl;
    std::cout << "Xv : " << Xv << std::endl;
    std::cout << "Yu : " << Yu << std::endl;
    std::cout << "Yv : " << Yv << std::endl;
    std::cout << "m_w_tilde : " << m_w_tilde << std::endl;

    // Decomposition of the Weingarten matrix
    Eigen::SelfAdjointEigenSolver<Matrix2> eigensolver(W);
    if (eigensolver.info() != Eigen::Success) {
        std::cerr << "EigenSolver failed" << std::endl;
    }

    Vector2 eivals = eigensolver.eigenvalues().real();
    Matrix2 eivecs = eigensolver.eigenvectors().real();

    if ( eivals(0) < eivals(1) ) {
        m_kmin = eivals(0);
        m_kmax = eivals(1);

        VectorType vmin = VectorType(0, eivecs.col(0)(0), eivecs.col(0)(1));
        VectorType vmax = VectorType(0, eivecs.col(1)(0), eivecs.col(1)(1));

        m_dmin = Base::template localFrameToWorld<true>(vmin);
        m_dmax = Base::template localFrameToWorld<true>(vmax);
    }
    else {
        m_kmin = eivals(1);
        m_kmax = eivals(0);

        VectorType vmin = VectorType(0, eivecs.col(1)(0), eivecs.col(1)(1));
        VectorType vmax = VectorType(0, eivecs.col(0)(0), eivecs.col(0)(1));

        m_dmin = Base::template localFrameToWorld<true>(vmin);
        m_dmax = Base::template localFrameToWorld<true>(vmax);
    }
        

}   


template < class DataPoint, class _NFilter, typename T>
typename NormalWeingarten<DataPoint, _NFilter, T>::Vector4
NormalWeingarten<DataPoint, _NFilter, T>::projectOnto4DPlane(const Vector4& point, const Vector4& normal, const Scalar& offset) const {
    Scalar t = ( offset - point.dot(normal) ) / normal.squaredNorm();
    return point + t * normal;
}


template < class DataPoint, class _NFilter, typename T>
typename NormalWeingarten<DataPoint, _NFilter, T>::Scalar
NormalWeingarten<DataPoint, _NFilter, T>::kMean() const {
  return ( kmin() + kmax() ) / Scalar(2);
}

template < class DataPoint, class _NFilter, typename T>
typename NormalWeingarten<DataPoint, _NFilter, T>::Scalar
NormalWeingarten<DataPoint, _NFilter, T>::GaussianCurvature() const {
    return kmin() * kmax();
}

template < class DataPoint, class _NFilter, typename T>
typename NormalWeingarten<DataPoint, _NFilter, T>::Scalar
NormalWeingarten<DataPoint, _NFilter, T>::kmin() const {
    return m_kmin;
}

template < class DataPoint, class _NFilter, typename T>
typename NormalWeingarten<DataPoint, _NFilter, T>::Scalar
NormalWeingarten<DataPoint, _NFilter, T>::kmax() const {
    return m_kmax;
}

template < class DataPoint, class _NFilter, typename T>
typename NormalWeingarten<DataPoint, _NFilter, T>::VectorType
NormalWeingarten<DataPoint, _NFilter, T>::kminDirection() const {
    return m_dmin;
}

template < class DataPoint, class _NFilter, typename T>
typename NormalWeingarten<DataPoint, _NFilter, T>::VectorType
NormalWeingarten<DataPoint, _NFilter, T>::kmaxDirection() const {
    return m_dmax;
}
