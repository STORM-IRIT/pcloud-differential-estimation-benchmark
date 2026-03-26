template < class DataPoint, class _NFilter, typename T>
void
QuadricFitImpl<DataPoint, _NFilter, T>::init()
{
    Base::init();
    m_matA = MatrixA::Zero();
    m_minId = -1;
}

template < class DataPoint, class _NFilter, typename T>
bool
QuadricFitImpl<DataPoint, _NFilter, T>::addLocalNeighbor(Scalar w,
                                                     const VectorType &localQ,
                                                     const DataPoint &attributes)
{
    if( Base::addLocalNeighbor(w, localQ, attributes) ) {
        VectorA a = Base::convertToParameters(localQ);
        m_matA     += w * a * a.transpose();
        return true;
    }

    return false;
}


template < class DataPoint, class _NFilter, typename T>
FIT_RESULT
QuadricFitImpl<DataPoint, _NFilter, T>::finalize ()
{
    FIT_RESULT res = Base::finalize();

    if( res != STABLE )
        return Base::m_eCurrentState;

    if( Base::getNumNeighbors() < 9 )
        return Base::m_eCurrentState = UNDEFINED;
    
    m_solver.compute( m_matA );

    VectorA eivals = m_solver.eigenvalues().real();

    m_minId = 0;

    for(int i=0 ; i < eivals.size() ; ++i)
    {
        Scalar ev = eivals(i);
        if( ev < eivals(m_minId) )
            m_minId = i;
    }

    Base::m_coefficients = m_solver.eigenvectors().col(m_minId).real();
    Base::computeFrameFromNormalVector( Base::primitiveGradient() );

    FIT_RESULT solverRes = constructTensor();
    if ( solverRes != STABLE )
        return res = solverRes;

    return res = Base::m_eCurrentState;
}

template < class DataPoint, class _NFilter, typename T>
FIT_RESULT
QuadricFitImpl<DataPoint, _NFilter, T>::constructTensor()
{
  VectorType x = Base::project( Base::getNeighborFilter().evalPos() );

  VectorType g = Base::primitiveGradient(x);
  Scalar g_norm = g.norm();

  if (g_norm < 1e-6) return UNDEFINED;

  VectorType n = g / g_norm;

  MatrixType H;
  H << Base::f_xx(), Base::f_xy(), Base::f_xz(),
       Base::f_xy(), Base::f_yy(), Base::f_yz(),
       Base::f_xz(), Base::f_yz(), Base::f_zz();

  MatrixType I3 = MatrixType::Identity();
  MatrixType P = I3 - n * n.transpose();

  MatrixType W = - (P * H * P) / g_norm;

  Eigen::SelfAdjointEigenSolver<MatrixType> eig(W);
  VectorType evals = eig.eigenvalues().real();
  MatrixType evecs = eig.eigenvectors().real();

  int zero_idx = 0;
  if (std::abs(evals(1)) < std::abs(evals(zero_idx))) zero_idx = 1;
  if (std::abs(evals(2)) < std::abs(evals(zero_idx))) zero_idx = 2;

  int id1 = (zero_idx + 1) % 3;
  int id2 = (zero_idx + 2) % 3;

  if (evals(id1) < evals(id2)) std::swap(id1, id2);

  Base::m_kmax = evals(id1);
  Base::m_kmin = evals(id2);

  Base::m_dmax = evecs.col(id1).normalized();
  Base::m_dmin = evecs.col(id2).normalized();

  return STABLE;
}

template<class DataPoint, class _NFilter, typename T>
typename QuadricFitImpl<DataPoint, _NFilter, T>::Matrix32 
QuadricFitImpl<DataPoint, _NFilter, T>::tangentPlane(const VectorType& n)
{
    Matrix32 B;
    B.col(0) = Base::getFrameU();
    B.col(1) = Base::getFrameV();
    return B;
}
