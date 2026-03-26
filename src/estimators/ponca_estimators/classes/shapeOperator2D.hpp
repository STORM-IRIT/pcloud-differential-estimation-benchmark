
template < class DataPoint, class _NFilter, typename T >
void ShapeOperator2DFitImpl<DataPoint, _NFilter, T>::init()
{
    Base::init();

    m_planeIsReady = false;

    m_P.setZero();
    m_A.setZero();
    m_B.setZero();
    m_W.setZero();

    m_kmin = 0;
    m_kmax = 0;
    m_dmin.setZero();
    m_dmax.setZero();
}

template < class DataPoint, class _NFilter, typename T >
bool ShapeOperator2DFitImpl<DataPoint, _NFilter, T>::addLocalNeighbor(
    Scalar w, const VectorType& localQ, const DataPoint& attributes)
{
    auto res = Base::addLocalNeighbor(w, localQ, attributes);
    if( ! m_planeIsReady ){
        return res;
    } 
    else {
        const VectorType dN = attributes.normal() - Base::normal();
        const Vector2 dN_proj = m_P.transpose() * dN;
        const Vector2 dP_proj = m_P.transpose() * localQ;
        m_A += w * dP_proj * dP_proj.transpose();
        m_B += w * dN_proj * dP_proj.transpose();
        return true;
    }
    return false;
}

// static
template < class DataPoint, class _NFilter, typename T >
bool ShapeOperator2DFitImpl<DataPoint, _NFilter, T>::solve(
    const Matrix22& A,
    const Matrix22& B,
    Matrix22& X)
{
    {
        Matrix22 invA = Matrix22::Zero();
        Scalar det = 0;
        bool invertible = false;
        A.computeInverseAndDetWithCheck(invA, det, invertible);
        if(invertible)
            X = B * invA;
        return invertible;
    }
}


template < class DataPoint, class _NFilter, typename T >
FIT_RESULT ShapeOperator2DFitImpl<DataPoint, _NFilter, T>::finalize()
{

    if ( !m_planeIsReady ){
        FIT_RESULT res = Base::finalize();
        if(res == STABLE) {
            Base::computeFrameFromNormalVector(Base::plane().primitiveGradient());
            m_P = tangentPlane();
            m_planeIsReady = true;
            return Base::m_eCurrentState = NEED_OTHER_PASS;
        }
        return res;
    }
    else {
        if(Base::getNumNeighbors() < 3) {
            std::cout << "ShapeOperator2DFitImpl<DataPoint, _NFilter, T>::finalize() : Base::finalize() != STABLE  || Base::getNumNeighbors() < 3" << std::endl;
            return Base::m_eCurrentState = UNSTABLE;
        }
        const bool ok = solve(m_A, m_B, m_W);
        if(not ok) {
            std::cout << "ShapeOperator2DFitImpl<DataPoint, _NFilter, T>::finalize() : not ok" << std::endl;
            return Base::m_eCurrentState = UNSTABLE;
        }
        // symmetrize
        m_W(0,1) = m_W(1,0) = (m_W(0,1) + m_W(1,0))/Scalar(2);
        Eigen::SelfAdjointEigenSolver<Matrix22> solver;
        solver.computeDirect(m_W);
        if(solver.info() != Eigen::Success) {
            std::cout << "ShapeOperator2DFitImpl<DataPoint, _NFilter, T>::finalize() : solver.info() != Eigen::Success" << std::endl;
            return Base::m_eCurrentState = UNSTABLE;
        }
        m_kmin = solver.eigenvalues()[0];
        m_dmin = m_P * solver.eigenvectors().col(0);
        m_kmax = solver.eigenvalues()[1];
        m_dmax = m_P * solver.eigenvectors().col(1);
        return Base::m_eCurrentState = STABLE;
    }
}


template<class DataPoint, class _NFilter, typename T>
typename ShapeOperator2DFitImpl<DataPoint, _NFilter, T>::Matrix32 
ShapeOperator2DFitImpl<DataPoint, _NFilter, T>::tangentPlane() const
{
    Matrix32 B;
    B.col(0) = Base::getFrameU();
    B.col(1) = Base::getFrameV();
    return B;
}
