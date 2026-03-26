#pragma once


#include <Eigen/Core>
#include <Ponca/Fitting>

namespace Ponca
{

/*!
    \brief Algebraic quadric primitive

    \see Quadric
*/
template < class DataPoint, class _NFilter, typename T >
class Quadric : public T
{
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    PONCA_FITTING_DECLARE_MATRIX_TYPE

protected:
    enum
    {
        check = Base::PROVIDES_PRIMITIVE_BASE &&
                Base::PROVIDES_LOCAL_FRAME,
        PROVIDES_ALGEBRAIC_QUADRIC,
        PROVIDES_NORMAL_DERIVATIVE,
    };

public:
    PONCA_EXPLICIT_CAST_OPERATORS(Quadric, quadric)

public:

    using VectorA = Eigen::Matrix<Scalar, 10, 1>;
    using MatrixA = Eigen::Matrix<Scalar, 10, 10>;
// results
public:

    VectorA m_coefficients;
    
    // Curvature info
    VectorType m_dmin;
    VectorType m_dmax;

    Scalar     m_kmin;
    Scalar     m_kmax;

public:

    PONCA_MULTIARCH inline void init()
    {
        Base::init();

        m_coefficients = VectorA::Zero();

        m_dmin = VectorType::Zero();
        m_dmax = VectorType::Zero();
        
        m_kmin = Scalar(0);
        m_kmax = Scalar(0);
    }

    PONCA_MULTIARCH inline VectorA convertToParameters(const VectorType& _q) const{
        return VectorA( _q(0)*_q(0), _q(1)*_q(1), _q(2)*_q(2), _q(0)*_q(1), _q(1)*_q(2), _q(0)*_q(2), _q(0), _q(1), _q(2), Scalar(1) );
    }

    //! \brief Value of the scalar field at the location \f$ \mathbf{q} \f$
    PONCA_MULTIARCH inline Scalar potential (const VectorType& _q) const;

    /*! \brief Value of the scalar field at the evaluation point */
    PONCA_MULTIARCH inline Scalar potential() const { return potential( Base::m_w.basisCenter() ); }

    //! \brief Project a point on the ellipsoid
    PONCA_MULTIARCH inline VectorType project (const VectorType& _q) const;

    //! \brief Approximation of the scalar field gradient at \f$ \mathbf{q} \f$
    PONCA_MULTIARCH inline VectorType primitiveGradient (const VectorType& _q) const;

    /*! \brief Approximation of the scalar field gradient at the evaluation point */
    PONCA_MULTIARCH inline VectorType primitiveGradient () const {
        VectorType n = primitiveGradient( project( Base::getNeighborFilter().evalPos() ) );
        return n / n.norm();
    }
    
    PONCA_MULTIARCH inline const Scalar f_a ( ) const { return m_coefficients(0); }
    PONCA_MULTIARCH inline const Scalar f_b ( ) const { return m_coefficients(1); }
    PONCA_MULTIARCH inline const Scalar f_c ( ) const { return m_coefficients(2); }
    PONCA_MULTIARCH inline const Scalar f_e ( ) const { return m_coefficients(3); }
    PONCA_MULTIARCH inline const Scalar f_f ( ) const { return m_coefficients(4); }
    PONCA_MULTIARCH inline const Scalar f_g ( ) const { return m_coefficients(5); }
    PONCA_MULTIARCH inline const Scalar f_l ( ) const { return m_coefficients(6); }
    PONCA_MULTIARCH inline const Scalar f_m ( ) const { return m_coefficients(7); }
    PONCA_MULTIARCH inline const Scalar f_n ( ) const { return m_coefficients(8); }
    PONCA_MULTIARCH inline const Scalar f_d ( ) const { return m_coefficients(9); }


    PONCA_MULTIARCH inline const Scalar f_xx () const { return Scalar(2) * f_a(); }
    PONCA_MULTIARCH inline const Scalar f_yy () const { return Scalar(2) * f_b(); }
    PONCA_MULTIARCH inline const Scalar f_zz () const { return Scalar(2) * f_c(); }
    PONCA_MULTIARCH inline const Scalar f_xy () const { return f_e(); }
    PONCA_MULTIARCH inline const Scalar f_yz () const { return f_f(); }
    PONCA_MULTIARCH inline const Scalar f_xz () const { return f_g(); }
    
    PONCA_MULTIARCH inline const Scalar f_x  ( VectorType _q = VectorType::Zero() ) const { 
        Scalar two = Scalar(2);
        Scalar x = _q(0);
        Scalar y = _q(1);
        Scalar z = _q(2);
        return two * f_a() * x + f_e() * y + f_g() * z + f_l();
    }
    PONCA_MULTIARCH inline const Scalar f_y  ( VectorType _q = VectorType::Zero() ) const { 
        Scalar two = Scalar(2);
        Scalar x = _q(0);
        Scalar y = _q(1);
        Scalar z = _q(2);
        return two * f_b() * y + f_e() * x + f_f() * z + f_m();
    }
    PONCA_MULTIARCH inline const Scalar f_z  ( VectorType _q = VectorType::Zero() ) const { 
        Scalar two = Scalar(2);
        Scalar x = _q(0);
        Scalar y = _q(1);
        Scalar z = _q(2);
        return two * f_c() * z + f_f() * y + f_g() * x + f_n();
    }

    PONCA_MULTIARCH inline const MatrixType dNormal() const {

        VectorType n = primitiveGradient();
        MatrixType mat;
        mat << f_xx(), f_xy(), f_xz(), 
               f_xy(), f_yy(), f_yz(), 
               f_xz(), f_yz(), f_zz();
        mat /= n.norm();
        return mat;
    }

    PONCA_MULTIARCH inline const Scalar magnitude ( VectorType _q = VectorType::Zero() ) const { 
        PONCA_MULTIARCH_STD_MATH(pow);
        PONCA_MULTIARCH_STD_MATH(sqrt);
        Scalar two = Scalar(2);
        Scalar fx2 = pow ( f_x(_q), two );
        Scalar fy2 = pow ( f_y(_q), two );
        Scalar fz2 = pow ( f_z(_q), two );
        return sqrt ( fx2 + fy2 + fz2 );
    }

    PONCA_MULTIARCH inline const Scalar dE ( VectorType _q = VectorType::Zero() ) const { 
        PONCA_MULTIARCH_STD_MATH(pow);
        Scalar two = Scalar(2);
        Scalar fx2 = pow ( f_x(_q), two );
        Scalar fz2 = pow ( f_z(_q), two ); 
        if ( fz2 <= 1e-6 ) fz2 = Scalar(1);
        return Scalar(1) + ( fx2 / fz2 );
    }

    PONCA_MULTIARCH inline const Scalar dF ( VectorType _q = VectorType::Zero() ) const {
      PONCA_MULTIARCH_STD_MATH(pow);
      Scalar two = Scalar(2);
      Scalar fz2 = pow ( f_z(_q), two );
      if ( fz2 <= 1e-6 ) fz2 = Scalar(1);
      // return Scalar(1) + ( f_x(_q) * f_y(_q) ) / fz2;
      return ( f_x(_q) * f_y(_q) ) / fz2;
    }

    PONCA_MULTIARCH inline const Scalar dG ( VectorType _q = VectorType::Zero() ) const { 
        PONCA_MULTIARCH_STD_MATH(pow);
        Scalar two = Scalar(2);
        Scalar fy2 = pow ( f_y(_q), two );
        Scalar fz2 = pow ( f_z(_q), two ); 
        if ( fz2 <= 1e-6 ) fz2 = Scalar(1);
        return Scalar(1) + ( fy2 / fz2 );
    }

    PONCA_MULTIARCH inline const Scalar dL ( VectorType _q = VectorType::Zero() ) const { 
        PONCA_MULTIARCH_STD_MATH(pow);
        static const Scalar one (1);
        static const Scalar two (2);
        Scalar fz2 = pow ( f_z(_q), two );
        if ( fz2 <= 1e-6 ) fz2 = Scalar(1);

        MatrixType mat; 
        mat << f_xx(), f_xz(), f_x(_q), f_xz(), f_zz(), f_z(_q), f_x(_q), f_z(_q), Scalar(0);  
        
        return ( Scalar(1) / ( fz2 * magnitude(_q) ) ) * mat.determinant();
    }

    PONCA_MULTIARCH inline const Scalar dM ( VectorType _q = VectorType::Zero() ) const { 
        PONCA_MULTIARCH_STD_MATH(pow);
        static const Scalar one (1);
        static const Scalar two (2);
        Scalar fz2 = pow ( f_z(_q), two );
        if ( fz2 <= 1e-6 ) fz2 = Scalar(1);

        MatrixType mat; 
        mat << f_xy(), f_yz(), f_y(_q), f_xz(), f_zz(), f_z(_q), f_x(_q), f_z(_q), Scalar(0);  

        return ( Scalar(1) / ( fz2 * magnitude(_q) ) ) * mat.determinant();
    }

    PONCA_MULTIARCH inline const Scalar dN ( VectorType _q = VectorType::Zero() ) const { 
        PONCA_MULTIARCH_STD_MATH(pow);
        static const Scalar one (1);
        static const Scalar two (2);
        Scalar fz2 = pow ( f_z(_q), two );
        if ( fz2 <= 1e-6 ) fz2 = Scalar(1);

        MatrixType mat; 
        mat << f_yy(), f_yz(), f_y(_q), f_yz(), f_zz(), f_z(_q), f_y(_q), f_z(_q), Scalar(0);  

        return ( Scalar(1) / ( fz2 * magnitude(_q) ) ) * mat.determinant();
    }

    PONCA_MULTIARCH inline Scalar kmin () const { return m_kmin; }
    PONCA_MULTIARCH inline Scalar kmax () const { return m_kmax; }
    PONCA_MULTIARCH inline Scalar kMean () const { return (m_kmin + m_kmax) / Scalar(2); }
    PONCA_MULTIARCH inline Scalar GaussianCurvature () const { return m_kmin * m_kmax; }
    PONCA_MULTIARCH inline VectorType kminDirection() const { return m_dmin; }
    PONCA_MULTIARCH inline VectorType kmaxDirection() const { return m_dmax; }


}; //class Quadric

#include "quadric.hpp"
}
