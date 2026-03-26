template < class DataPoint, class _NFilter, typename T>
template <bool ignoreTranslation>
typename LocalFrame<DataPoint, _NFilter, T>::VectorType
LocalFrame<DataPoint, _NFilter, T>::worldToLocalFrame (const VectorType& _q) const
{
  MatrixType B;
  B << Base::primitiveGradient(), m_u, m_v;
  if (ignoreTranslation)
    return B.transpose() * _q;
  else {
    // apply rotation and translation to get uv coordinates
    return B.transpose() * (Base::getNeighborFilter().convertToLocalBasis(_q));
  }
}

template < class DataPoint, class _NFilter, typename T>
template <bool ignoreTranslation>
typename LocalFrame<DataPoint, _NFilter, T>::VectorType
LocalFrame<DataPoint, _NFilter, T>::localFrameToWorld (const VectorType& _lq) const
{
  MatrixType B;
  B << Base::primitiveGradient(), m_u, m_v;
  if (ignoreTranslation)
    return B * _lq;
  else {
    return B * _lq + Base::getNeighborFilter().evalPos();
  }
}

template < class DataPoint, class _NFilter, typename T>
void
LocalFrame<DataPoint, _NFilter, T>::computeFrameFromNormalVector (const VectorType& _norm)
{
    // Cas 2D
    if constexpr (VectorType::RowsAtCompileTime == 2) {
        VectorType m_u(-_norm.y(), _norm.x()); // Vecteur perpendiculaire en 2D
        m_u.normalize();
        setFrameUV(m_u, -m_u);
    } 
    
    // Cas 3D
    else if constexpr (VectorType::RowsAtCompileTime == 3) {
      // Creation of the vector 'a' non-collinear to the normal vector
      VectorType a;
      if (std::abs(_norm.x()) > std::abs(_norm.z())) {
          a = VectorType(-_norm.y(), _norm.x(), 0);
      } else {
          a = VectorType(0, -_norm.z(), _norm.y());
      }
      a.normalize();
      // Creation of the two vectors of the local frame (m_u and m_v) thanks to the cross product
      VectorType m_u = _norm.cross(a);
      VectorType m_v = _norm.cross(m_u);
      m_u.normalize();
      m_v.normalize();
      setFrameUV (m_u, m_v);
    }
}
