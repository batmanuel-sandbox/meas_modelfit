#ifndef NDARRAY_vectorize_hpp_INCLUDED
#define NDARRAY_vectorize_hpp_INCLUDED

/** 
 *  \file ndarray/vectorize.hpp \brief Code to apply arbitrary scalar functors to arrays.
 */

#include "ndarray_fwd.hpp"
#include "ndarray/detail/UnaryOp.hpp"

#include <boost/mpl/and.hpp>

namespace ndarray {
namespace result_of {

template <typename T1, typename T2, typename T3=void>
struct vectorize {
    typedef T1 BinaryFunction;
    typedef T2 Argument1;
    typedef T3 Argument2;

    typename boost::mpl::if_<
        boost::mpl::and_<
            typename detail::ExpressionTraits<Argument1>::IsScalar,
            typename detail::ExpressionTraits<Argument2>::IsScalar
        >,
        typename BinaryFunction::result_type,
        detail::BinaryOpExpression<Argument1,Argument2,BinaryFunction>
    >::type type;
    
};

template <typename T1, typename T2>
struct vectorize<T1,T2,void> {
    typedef T1 UnaryFunction;
    typedef T2 Argument;

    typedef typename boost::mpl::if_<
        typename detail::ExpressionTraits<Argument>::IsScalar,
        typename UnaryFunction::result_type,
        detail::UnaryOpExpression<Argument,UnaryFunction>
    >::type type;
};

} // namespace ndarray::result_of

/// \addtogroup MainGroup
/// @{

/** 
 *  \brief Apply a non-mutating unary function object to a scalar.
 *
 *  This overload exists to allow recursive usage of the Array-argument vectorize functions.
 */
template <typename Scalar, typename UnaryFunction>
#ifndef DOXYGEN
typename boost::enable_if<typename detail::ExpressionTraits<Scalar>::IsScalar,
                          typename UnaryFunction::result_type>::type
#else
typename UnaryFunction::result_type
#endif
vectorize(
    UnaryFunction const & functor,
    Scalar const & scalar
) {
    return functor(scalar);
}

/** 
 *  \brief Apply a non-mutating unary function object to each element of a multidimensional Expression.
 *
 *  Evaluation is lazy.
 */
template <typename Derived, typename UnaryFunction>
detail::UnaryOpExpression<Derived,UnaryFunction>
vectorize(
    UnaryFunction const & functor,
    Expression<Derived> const & operand
) {
    return detail::UnaryOpExpression<Derived,UnaryFunction>(
        static_cast<Derived const &>(operand),
        functor
    );
}

/** 
 *  \brief Apply a non-mutating binary function object to a pair of scalars.
 *
 *  This overload exists to allow recursive usage of the Array-argument vectorize functions.
 */
template <typename Scalar1, typename Scalar2, typename BinaryFunction>
#ifndef DOXYGEN
typename boost::enable_if_c<
    (detail::ExpressionTraits<Scalar1>::IsScalar::value
     && detail::ExpressionTraits<Scalar2>::IsScalar::value),
    typename BinaryFunction::result_type
    >::type
#else
typename BinaryFunction::result_type
#endif
vectorize(
    BinaryFunction const & functor,
    Scalar1 const & scalar1,
    Scalar2 const & scalar2
) {
    return functor(scalar1,scalar2);
}

/** 
 *  \brief Apply a non-mutating binary function object pairwise to 
 *  the elements of two multidimensional Expressions.
 *
 *  Evaluation is lazy.
 */
template <typename Derived1, typename Derived2, typename BinaryFunction>
detail::BinaryOpExpression<Derived1,Derived2,BinaryFunction>
vectorize(
    BinaryFunction const & functor,
    Expression<Derived1> const & operand1,
    Expression<Derived2> const & operand2
) {
    return detail::BinaryOpExpression<Derived1,Derived2,BinaryFunction>(
        static_cast<Derived1 const &>(operand1),
        static_cast<Derived2 const &>(operand2),
        functor
    );
}

/// @}

} // namespace ndarray

#endif // !NDARRAY_vectorize_hpp_INCLUDED