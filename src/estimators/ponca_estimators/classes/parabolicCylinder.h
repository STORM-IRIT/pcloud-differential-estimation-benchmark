#pragma once

#include <Eigen/Core>
#include <Ponca/Fitting>

namespace Ponca
{

/*!
    \brief Algebraic parabolical-cylinder primitive

    This primitive provides:
    \verbatim PROVIDES_PARABOLIC_CYLINDER \endverbatim

    \see ParabolicCylinder
*/
template < class DataPoint, class _NFilter, typename T >
class ParabolicCylinder : public T
{
    PONCA_FITTING_DECLARE_DEFAULT_TYPES
    PONCA_FITTING_DECLARE_MATRIX_TYPE

protected:
    enum
    {
        check = Base::PROVIDES_PRIMITIVE_BASE 
        && Base::PROVIDES_PLANE,                                        /*!< \brief Requires PrimitiveBase      */
        PROVIDES_PARABOLIC_CYLINDER,                                    /*!< \brief Provides Parabolic cylinder */
    };

public:
    PONCA_EXPLICIT_CAST_OPERATORS(ParabolicCylinder,parabolicCylinder)

public:

// results
public:

    using Matrix2       = Eigen::Matrix<Scalar, 2, 2>;
    using Vector2       = Eigen::Matrix<Scalar, 2, 1>;



    // f(q) = m_uc + m_ul^T * q + m_a * q^T * (m_uq^T * m_uq) * q;
    Scalar     m_uc;       /*!< \brief Constant parameter of the Algebraic hyper-ellipsoid   */
    Vector2    m_ul;       /*!< \brief Linear parameter of the Algebraic hyper-ellipsoid     */
    Scalar     m_a;        /*!< \brief Alpha quadratic parameter, used for ambiguous point   */
    Matrix2    m_uq;       /*!< \brief Quadratic parameter of the Algebraic hyper-ellipsoid  */

    // Fix to correct the orientation of the primitive
    Scalar m_correctOrientation = Scalar(1);

    // Curvature info
    Scalar     m_k1;
    Scalar     m_k2;
    VectorType m_v1;
    VectorType m_v2;

public:
    PONCA_MULTIARCH inline void init()
    {
        Base::init();

        m_uc = Scalar(0);
        m_ul = Vector2::Zero();
        m_a  = Scalar(1);
        m_uq = Matrix2::Zero();

        m_correctOrientation = Scalar(1);

        m_v1 = VectorType::Zero();
        m_v2 = VectorType::Zero();
        m_k1 = Scalar(0);
        m_k2 = Scalar(0);
    }

    PONCA_MULTIARCH inline Scalar eval_quadratic_function(Scalar q1, Scalar q2) const {
        Vector2 q {q1, q2};
        return h_c() + h_u()*q1 + h_v()*q2 + m_a * ( h_uu() * q1 * q1 + h_vv() * q2 * q2 + Scalar(2) * h_uv() * q1 * q2 );
    }

    PONCA_MULTIARCH inline Vector2 evaldUV(Scalar q1, Scalar q2) const {
      return Vector2( h_u() + Scalar(2) * m_a * ( h_uu() * q1 +  h_uv() * q2 ),
                      h_v() + Scalar(2) * m_a * ( h_vv() * q2  +  h_uv() * q1 ) );
    }

    //! \brief Value of the scalar field at the location \f$ \mathbf{q} \f$
    PONCA_MULTIARCH inline Scalar potential (const VectorType& _q) const;

    /*! \brief Value of the scalar field at the evaluation point */
    PONCA_MULTIARCH inline Scalar potential() const { return potential(Base::getNeighborFilter().evalPos()); }

    //! \brief Project a point on the ellipsoid
    PONCA_MULTIARCH inline VectorType project (const VectorType& _q) const;

    //! \brief Approximation of the scalar field gradient at \f$ \mathbf{q} \f$
    PONCA_MULTIARCH inline VectorType primitiveGradient (const VectorType& _q) const;

    /*! \brief Approximation of the scalar field gradient at the evaluation point */
    PONCA_MULTIARCH inline VectorType primitiveGradient () const {
        VectorType n = primitiveGradient(Base::getNeighborFilter().evalPos());
        return n / n.norm();
    }

    PONCA_MULTIARCH inline const Scalar h_uu () const { return m_uq(0,0);   }
    PONCA_MULTIARCH inline const Scalar h_vv () const { return m_uq(1,1); }
    PONCA_MULTIARCH inline const Scalar h_uv () const { return m_uq(1,0); }
    PONCA_MULTIARCH inline const Scalar h_u  () const { return m_ul(0); }
    PONCA_MULTIARCH inline const Scalar h_v  () const { return m_ul(1); }
    PONCA_MULTIARCH inline const Scalar h_c  () const { return m_uc; }

    // Partial Derivatives at (u,v) = (0,0)
    PONCA_MULTIARCH inline const Scalar dh_uu () const { return Scalar(-2) * m_a * h_uu(); }
    PONCA_MULTIARCH inline const Scalar dh_vv () const { return Scalar(-2) * m_a * h_vv(); }
    PONCA_MULTIARCH inline const Scalar dh_uv () const { return Scalar(-2) * m_a * h_uv(); }
    PONCA_MULTIARCH inline const Scalar dh_u  () const { return h_u(); }
    PONCA_MULTIARCH inline const Scalar dh_v  () const { return h_v(); }
    PONCA_MULTIARCH inline const Scalar dh_c  () const { return Scalar(0); }

    PONCA_MULTIARCH inline const Scalar dE () const { return Scalar (1) + ( dh_u() * dh_u() ); }
    PONCA_MULTIARCH inline const Scalar dF () const { return dh_u() * dh_v(); }
    PONCA_MULTIARCH inline const Scalar dG () const { return Scalar (1) + ( dh_v() * dh_v() ); }
    PONCA_MULTIARCH inline const Scalar dL () const { 
        PONCA_MULTIARCH_STD_MATH(pow);
        static const Scalar one (1);
        static const Scalar two (2);
        static const Scalar oneOverTwo (Scalar(1)/Scalar(2));
        return ( dh_uu() ) / ( pow( pow( dh_u(), two) + pow( dh_v(), two) + one , oneOverTwo ) );
    }
    PONCA_MULTIARCH inline const Scalar dM () const { 
        PONCA_MULTIARCH_STD_MATH(pow);
        static const Scalar one (1);
        static const Scalar two (2);
        static const Scalar oneOverTwo (Scalar(1)/Scalar(2));
        return ( dh_uv() ) / ( pow( pow( dh_u(), two) + pow( dh_v(), two) + one , oneOverTwo ) );
    }
    PONCA_MULTIARCH inline const Scalar dN () const { 
        PONCA_MULTIARCH_STD_MATH(pow);
        static const Scalar one (1);
        static const Scalar two (2);
        static const Scalar oneOverTwo (Scalar(1)/Scalar(2));
        return ( dh_vv() ) / ( pow( pow( dh_u(), two) + pow( dh_v(), two) + one , oneOverTwo ) );
    }

    PONCA_MULTIARCH inline Scalar alpha_curvature () const { return m_a; }
    PONCA_MULTIARCH inline Scalar kmin () const;
    PONCA_MULTIARCH inline Scalar kmax () const;
    PONCA_MULTIARCH inline Scalar kMean () const;
    PONCA_MULTIARCH inline Scalar GaussianCurvature () const;
    PONCA_MULTIARCH inline VectorType kminDirection() const;
    PONCA_MULTIARCH inline VectorType kmaxDirection() const;


}; //class ParabolicCylinder

#include "parabolicCylinder.hpp"
}
