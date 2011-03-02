
#ifndef NDARRAY_operators_hpp_INCLUDED
#define NDARRAY_operators_hpp_INCLUDED

/** 
 *  \file ndarray/operators.hpp \brief Arithmetic and logical operators for Array.
 */

#include "ndarray/Array.hpp"
#include <boost/call_traits.hpp>
#include <boost/functional.hpp>

#include "ndarray/detail/UnaryOp.hpp"
#include "ndarray/detail/BinaryOp.hpp"

namespace ndarray {

/// \cond INTERNAL
namespace detail {

/** 
 *  \internal @class PromotingBinaryFunction
 *  \brief A typedef-providing base class for binary functors with numeric type promotion.
 *
 *  \ingroup InternalGroup
 */
template <typename A, typename B>
struct PromotingBinaryFunction {
    typedef A ElementA;
    typedef B ElementB;
    typedef A first_argument_type;
    typedef B second_argument_type;
    typedef typename boost::call_traits<A>::param_type ParamA;
    typedef typename boost::call_traits<B>::param_type ParamB;
    typedef typename Promote<A,B>::Type result_type;
};

/** 
 *  \internal @class BinaryPredicate
 *  \brief A typedef-providing base class for binary predicates.
 *
 *  \ingroup InternalGroup
 */
template <typename A, typename B>
struct BinaryPredicate {
    typedef A ElementA;
    typedef B ElementB;
    typedef A first_argument_type;
    typedef B second_argument_type;
    typedef typename boost::call_traits<A>::param_type ParamA;
    typedef typename boost::call_traits<B>::param_type ParamB;
    typedef bool result_type;
};

/** 
 *  \internal @class AdaptableFunctionTag
 *  \brief A CRTP base class for non-template classes that contain a templated functor.
 *
 *  \ingroup InternalGroup
 */
template <typename Derived>
struct AdaptableFunctionTag {

    template <typename OperandB, typename A>
    struct ScalarExpr {
        typedef typename Derived::template ScalarFunction<
            A, typename detail::ExpressionTraits<OperandB>::Element
            > BinaryFunction;
        typedef boost::binder1st<BinaryFunction> Bound;
        static Bound bind(A const & scalar) {
            return Bound(BinaryFunction(),scalar);
        }
    };

    template <typename OperandA, typename B>
    struct ExprScalar {
        typedef typename Derived::template ScalarFunction<
            typename detail::ExpressionTraits<OperandA>::Element, B
            > BinaryFunction;
        typedef boost::binder2nd<BinaryFunction> Bound;
        static Bound bind(B const & scalar) {
            return Bound(BinaryFunction(),scalar);
        }
    };

    template <typename OperandA, typename OperandB>
    struct ExprExpr {
        typedef typename Derived::template ScalarFunction<
            typename detail::ExpressionTraits<OperandA>::Element,
            typename detail::ExpressionTraits<OperandB>::Element
            > BinaryFunction;
    };

};

/**
 *  \internal @class BitwiseNot
 *  \ingroup InternalGroup
 *  \brief An STL Unary Function class for bitwise NOT (unary ~).
 */
template <typename T>
struct BitwiseNot {
    typedef T argument_type;
    typedef T result_type;

    result_type operator()(argument_type arg) const { return ~arg; }
};


    struct PlusTag : public AdaptableFunctionTag<PlusTag> {
        template <typename A, typename B>
        struct ScalarFunction : public PromotingBinaryFunction<A,B> {
            typename PromotingBinaryFunction<A,B>::result_type operator()(
                typename PromotingBinaryFunction<A,B>::ParamA a,
                typename PromotingBinaryFunction<A,B>::ParamB b
            ) const {
                return a + b;
            }
        };
    };

    struct MinusTag : public AdaptableFunctionTag<MinusTag> {
        template <typename A, typename B>
        struct ScalarFunction : public PromotingBinaryFunction<A,B> {
            typename PromotingBinaryFunction<A,B>::result_type operator()(
                typename PromotingBinaryFunction<A,B>::ParamA a,
                typename PromotingBinaryFunction<A,B>::ParamB b
            ) const {
                return a - b;
            }
        };
    };

    struct MultipliesTag : public AdaptableFunctionTag<MultipliesTag> {
        template <typename A, typename B>
        struct ScalarFunction : public PromotingBinaryFunction<A,B> {
            typename PromotingBinaryFunction<A,B>::result_type operator()(
                typename PromotingBinaryFunction<A,B>::ParamA a,
                typename PromotingBinaryFunction<A,B>::ParamB b
            ) const {
                return a * b;
            }
        };
    };

    struct DividesTag : public AdaptableFunctionTag<DividesTag> {
        template <typename A, typename B>
        struct ScalarFunction : public PromotingBinaryFunction<A,B> {
            typename PromotingBinaryFunction<A,B>::result_type operator()(
                typename PromotingBinaryFunction<A,B>::ParamA a,
                typename PromotingBinaryFunction<A,B>::ParamB b
            ) const {
                return a / b;
            }
        };
    };

    struct ModulusTag : public AdaptableFunctionTag<ModulusTag> {
        template <typename A, typename B>
        struct ScalarFunction : public PromotingBinaryFunction<A,B> {
            typename PromotingBinaryFunction<A,B>::result_type operator()(
                typename PromotingBinaryFunction<A,B>::ParamA a,
                typename PromotingBinaryFunction<A,B>::ParamB b
            ) const {
                return a % b;
            }
        };
    };

    struct BitwiseXorTag : public AdaptableFunctionTag<BitwiseXorTag> {
        template <typename A, typename B>
        struct ScalarFunction : public PromotingBinaryFunction<A,B> {
            typename PromotingBinaryFunction<A,B>::result_type operator()(
                typename PromotingBinaryFunction<A,B>::ParamA a,
                typename PromotingBinaryFunction<A,B>::ParamB b
            ) const {
                return a ^ b;
            }
        };
    };

    struct BitwiseOrTag : public AdaptableFunctionTag<BitwiseOrTag> {
        template <typename A, typename B>
        struct ScalarFunction : public PromotingBinaryFunction<A,B> {
            typename PromotingBinaryFunction<A,B>::result_type operator()(
                typename PromotingBinaryFunction<A,B>::ParamA a,
                typename PromotingBinaryFunction<A,B>::ParamB b
            ) const {
                return a | b;
            }
        };
    };

    struct BitwiseAndTag : public AdaptableFunctionTag<BitwiseAndTag> {
        template <typename A, typename B>
        struct ScalarFunction : public PromotingBinaryFunction<A,B> {
            typename PromotingBinaryFunction<A,B>::result_type operator()(
                typename PromotingBinaryFunction<A,B>::ParamA a,
                typename PromotingBinaryFunction<A,B>::ParamB b
            ) const {
                return a | b;
            }
        };
    };

    struct BitwiseLeftShiftTag : public AdaptableFunctionTag<BitwiseLeftShiftTag> {
        template <typename A, typename B>
        struct ScalarFunction : public PromotingBinaryFunction<A,B> {
            typename PromotingBinaryFunction<A,B>::result_type operator()(
                typename PromotingBinaryFunction<A,B>::ParamA a,
                typename PromotingBinaryFunction<A,B>::ParamB b
            ) const {
                return a << b;
            }
        };
    };

    struct BitwiseRightShiftTag : public AdaptableFunctionTag<BitwiseRightShiftTag> {
        template <typename A, typename B>
        struct ScalarFunction : public PromotingBinaryFunction<A,B> {
            typename PromotingBinaryFunction<A,B>::result_type operator()(
                typename PromotingBinaryFunction<A,B>::ParamA a,
                typename PromotingBinaryFunction<A,B>::ParamB b
            ) const {
                return a >> b;
            }
        };
    };


    struct EqualToTag : public AdaptableFunctionTag<EqualToTag> {
        template <typename A, typename B>
        struct ScalarFunction : public BinaryPredicate<A,B> {
            typename BinaryPredicate<A,B>::result_type operator()(
                typename BinaryPredicate<A,B>::ParamA a,
                typename BinaryPredicate<A,B>::ParamB b
            ) const {
                return a == b;
            }
        };
    };

    struct NotEqualToTag : public AdaptableFunctionTag<NotEqualToTag> {
        template <typename A, typename B>
        struct ScalarFunction : public BinaryPredicate<A,B> {
            typename BinaryPredicate<A,B>::result_type operator()(
                typename BinaryPredicate<A,B>::ParamA a,
                typename BinaryPredicate<A,B>::ParamB b
            ) const {
                return a != b;
            }
        };
    };

    struct LessTag : public AdaptableFunctionTag<LessTag> {
        template <typename A, typename B>
        struct ScalarFunction : public BinaryPredicate<A,B> {
            typename BinaryPredicate<A,B>::result_type operator()(
                typename BinaryPredicate<A,B>::ParamA a,
                typename BinaryPredicate<A,B>::ParamB b
            ) const {
                return a < b;
            }
        };
    };

    struct GreaterTag : public AdaptableFunctionTag<GreaterTag> {
        template <typename A, typename B>
        struct ScalarFunction : public BinaryPredicate<A,B> {
            typename BinaryPredicate<A,B>::result_type operator()(
                typename BinaryPredicate<A,B>::ParamA a,
                typename BinaryPredicate<A,B>::ParamB b
            ) const {
                return a > b;
            }
        };
    };

    struct LessEqualTag : public AdaptableFunctionTag<LessEqualTag> {
        template <typename A, typename B>
        struct ScalarFunction : public BinaryPredicate<A,B> {
            typename BinaryPredicate<A,B>::result_type operator()(
                typename BinaryPredicate<A,B>::ParamA a,
                typename BinaryPredicate<A,B>::ParamB b
            ) const {
                return a <= b;
            }
        };
    };

    struct GreaterEqualTag : public AdaptableFunctionTag<GreaterEqualTag> {
        template <typename A, typename B>
        struct ScalarFunction : public BinaryPredicate<A,B> {
            typename BinaryPredicate<A,B>::result_type operator()(
                typename BinaryPredicate<A,B>::ParamA a,
                typename BinaryPredicate<A,B>::ParamB b
            ) const {
                return a >= b;
            }
        };
    };

    struct LogicalAnd : public AdaptableFunctionTag<LogicalAnd> {
        template <typename A, typename B>
        struct ScalarFunction : public BinaryPredicate<A,B> {
            typename BinaryPredicate<A,B>::result_type operator()(
                typename BinaryPredicate<A,B>::ParamA a,
                typename BinaryPredicate<A,B>::ParamB b
            ) const {
                return a && b;
            }
        };
    };

    struct LogicalOr : public AdaptableFunctionTag<LogicalOr> {
        template <typename A, typename B>
        struct ScalarFunction : public BinaryPredicate<A,B> {
            typename BinaryPredicate<A,B>::result_type operator()(
                typename BinaryPredicate<A,B>::ParamA a,
                typename BinaryPredicate<A,B>::ParamB b
            ) const {
                return a || b;
            }
        };
    };

} // namespace ndarray::detail
/// \endcond

/// \addtogroup OpGroup
/// @{

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::PlusTag::template ExprScalar<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator +(Expression<Operand> const & operand, Scalar const & scalar) {
        return vectorize(detail::PlusTag::template ExprScalar<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::PlusTag::template ScalarExpr<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator +(Scalar const & scalar, Expression<Operand> const & operand) {
        return vectorize(detail::PlusTag::template ScalarExpr<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand1, typename Operand2>
#ifndef DOXYGEN
    detail::BinaryOpExpression< 
         Operand1, Operand2,
         typename detail::PlusTag::template ExprExpr<Operand1,Operand2>::BinaryFunction
    >
#else
    <unspecified-expression-type>
#endif
    operator +(Expression<Operand1> const & operand1, Expression<Operand2> const & operand2) {
        return vectorize(
            typename detail::PlusTag::template ExprExpr<Operand1,Operand2>::BinaryFunction(),
            operand1,
            operand2
        );
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::MinusTag::template ExprScalar<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator -(Expression<Operand> const & operand, Scalar const & scalar) {
        return vectorize(detail::MinusTag::template ExprScalar<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::MinusTag::template ScalarExpr<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator -(Scalar const & scalar, Expression<Operand> const & operand) {
        return vectorize(detail::MinusTag::template ScalarExpr<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand1, typename Operand2>
#ifndef DOXYGEN
    detail::BinaryOpExpression< 
         Operand1, Operand2,
         typename detail::MinusTag::template ExprExpr<Operand1,Operand2>::BinaryFunction
    >
#else
    <unspecified-expression-type>
#endif
    operator -(Expression<Operand1> const & operand1, Expression<Operand2> const & operand2) {
        return vectorize(
            typename detail::MinusTag::template ExprExpr<Operand1,Operand2>::BinaryFunction(),
            operand1,
            operand2
        );
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::MultipliesTag::template ExprScalar<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator *(Expression<Operand> const & operand, Scalar const & scalar) {
        return vectorize(detail::MultipliesTag::template ExprScalar<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::MultipliesTag::template ScalarExpr<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator *(Scalar const & scalar, Expression<Operand> const & operand) {
        return vectorize(detail::MultipliesTag::template ScalarExpr<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand1, typename Operand2>
#ifndef DOXYGEN
    detail::BinaryOpExpression< 
         Operand1, Operand2,
         typename detail::MultipliesTag::template ExprExpr<Operand1,Operand2>::BinaryFunction
    >
#else
    <unspecified-expression-type>
#endif
    operator *(Expression<Operand1> const & operand1, Expression<Operand2> const & operand2) {
        return vectorize(
            typename detail::MultipliesTag::template ExprExpr<Operand1,Operand2>::BinaryFunction(),
            operand1,
            operand2
        );
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::DividesTag::template ExprScalar<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator /(Expression<Operand> const & operand, Scalar const & scalar) {
        return vectorize(detail::DividesTag::template ExprScalar<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::DividesTag::template ScalarExpr<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator /(Scalar const & scalar, Expression<Operand> const & operand) {
        return vectorize(detail::DividesTag::template ScalarExpr<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand1, typename Operand2>
#ifndef DOXYGEN
    detail::BinaryOpExpression< 
         Operand1, Operand2,
         typename detail::DividesTag::template ExprExpr<Operand1,Operand2>::BinaryFunction
    >
#else
    <unspecified-expression-type>
#endif
    operator /(Expression<Operand1> const & operand1, Expression<Operand2> const & operand2) {
        return vectorize(
            typename detail::DividesTag::template ExprExpr<Operand1,Operand2>::BinaryFunction(),
            operand1,
            operand2
        );
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::ModulusTag::template ExprScalar<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator %(Expression<Operand> const & operand, Scalar const & scalar) {
        return vectorize(detail::ModulusTag::template ExprScalar<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::ModulusTag::template ScalarExpr<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator %(Scalar const & scalar, Expression<Operand> const & operand) {
        return vectorize(detail::ModulusTag::template ScalarExpr<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand1, typename Operand2>
#ifndef DOXYGEN
    detail::BinaryOpExpression< 
         Operand1, Operand2,
         typename detail::ModulusTag::template ExprExpr<Operand1,Operand2>::BinaryFunction
    >
#else
    <unspecified-expression-type>
#endif
    operator %(Expression<Operand1> const & operand1, Expression<Operand2> const & operand2) {
        return vectorize(
            typename detail::ModulusTag::template ExprExpr<Operand1,Operand2>::BinaryFunction(),
            operand1,
            operand2
        );
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::BitwiseXorTag::template ExprScalar<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator ^(Expression<Operand> const & operand, Scalar const & scalar) {
        return vectorize(detail::BitwiseXorTag::template ExprScalar<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::BitwiseXorTag::template ScalarExpr<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator ^(Scalar const & scalar, Expression<Operand> const & operand) {
        return vectorize(detail::BitwiseXorTag::template ScalarExpr<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand1, typename Operand2>
#ifndef DOXYGEN
    detail::BinaryOpExpression< 
         Operand1, Operand2,
         typename detail::BitwiseXorTag::template ExprExpr<Operand1,Operand2>::BinaryFunction
    >
#else
    <unspecified-expression-type>
#endif
    operator ^(Expression<Operand1> const & operand1, Expression<Operand2> const & operand2) {
        return vectorize(
            typename detail::BitwiseXorTag::template ExprExpr<Operand1,Operand2>::BinaryFunction(),
            operand1,
            operand2
        );
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::BitwiseOrTag::template ExprScalar<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator |(Expression<Operand> const & operand, Scalar const & scalar) {
        return vectorize(detail::BitwiseOrTag::template ExprScalar<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::BitwiseOrTag::template ScalarExpr<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator |(Scalar const & scalar, Expression<Operand> const & operand) {
        return vectorize(detail::BitwiseOrTag::template ScalarExpr<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand1, typename Operand2>
#ifndef DOXYGEN
    detail::BinaryOpExpression< 
         Operand1, Operand2,
         typename detail::BitwiseOrTag::template ExprExpr<Operand1,Operand2>::BinaryFunction
    >
#else
    <unspecified-expression-type>
#endif
    operator |(Expression<Operand1> const & operand1, Expression<Operand2> const & operand2) {
        return vectorize(
            typename detail::BitwiseOrTag::template ExprExpr<Operand1,Operand2>::BinaryFunction(),
            operand1,
            operand2
        );
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::BitwiseAndTag::template ExprScalar<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator |(Expression<Operand> const & operand, Scalar const & scalar) {
        return vectorize(detail::BitwiseAndTag::template ExprScalar<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::BitwiseAndTag::template ScalarExpr<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator |(Scalar const & scalar, Expression<Operand> const & operand) {
        return vectorize(detail::BitwiseAndTag::template ScalarExpr<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand1, typename Operand2>
#ifndef DOXYGEN
    detail::BinaryOpExpression< 
         Operand1, Operand2,
         typename detail::BitwiseAndTag::template ExprExpr<Operand1,Operand2>::BinaryFunction
    >
#else
    <unspecified-expression-type>
#endif
    operator |(Expression<Operand1> const & operand1, Expression<Operand2> const & operand2) {
        return vectorize(
            typename detail::BitwiseAndTag::template ExprExpr<Operand1,Operand2>::BinaryFunction(),
            operand1,
            operand2
        );
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::BitwiseLeftShiftTag::template ExprScalar<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator <<(Expression<Operand> const & operand, Scalar const & scalar) {
        return vectorize(detail::BitwiseLeftShiftTag::template ExprScalar<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::BitwiseLeftShiftTag::template ScalarExpr<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator <<(Scalar const & scalar, Expression<Operand> const & operand) {
        return vectorize(detail::BitwiseLeftShiftTag::template ScalarExpr<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand1, typename Operand2>
#ifndef DOXYGEN
    detail::BinaryOpExpression< 
         Operand1, Operand2,
         typename detail::BitwiseLeftShiftTag::template ExprExpr<Operand1,Operand2>::BinaryFunction
    >
#else
    <unspecified-expression-type>
#endif
    operator <<(Expression<Operand1> const & operand1, Expression<Operand2> const & operand2) {
        return vectorize(
            typename detail::BitwiseLeftShiftTag::template ExprExpr<Operand1,Operand2>::BinaryFunction(),
            operand1,
            operand2
        );
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::BitwiseRightShiftTag::template ExprScalar<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator >>(Expression<Operand> const & operand, Scalar const & scalar) {
        return vectorize(detail::BitwiseRightShiftTag::template ExprScalar<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::BitwiseRightShiftTag::template ScalarExpr<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator >>(Scalar const & scalar, Expression<Operand> const & operand) {
        return vectorize(detail::BitwiseRightShiftTag::template ScalarExpr<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand1, typename Operand2>
#ifndef DOXYGEN
    detail::BinaryOpExpression< 
         Operand1, Operand2,
         typename detail::BitwiseRightShiftTag::template ExprExpr<Operand1,Operand2>::BinaryFunction
    >
#else
    <unspecified-expression-type>
#endif
    operator >>(Expression<Operand1> const & operand1, Expression<Operand2> const & operand2) {
        return vectorize(
            typename detail::BitwiseRightShiftTag::template ExprExpr<Operand1,Operand2>::BinaryFunction(),
            operand1,
            operand2
        );
    }


    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::EqualToTag::template ExprScalar<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator ==(Expression<Operand> const & operand, Scalar const & scalar) {
        return vectorize(detail::EqualToTag::template ExprScalar<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::EqualToTag::template ScalarExpr<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator ==(Scalar const & scalar, Expression<Operand> const & operand) {
        return vectorize(detail::EqualToTag::template ScalarExpr<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand1, typename Operand2>
#ifndef DOXYGEN
    detail::BinaryOpExpression< 
         Operand1, Operand2,
         typename detail::EqualToTag::template ExprExpr<Operand1,Operand2>::BinaryFunction
    >
#else
    <unspecified-expression-type>
#endif
    operator ==(Expression<Operand1> const & operand1, Expression<Operand2> const & operand2) {
        return vectorize(
            typename detail::EqualToTag::template ExprExpr<Operand1,Operand2>::BinaryFunction(),
            operand1,
            operand2
        );
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::NotEqualToTag::template ExprScalar<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator !=(Expression<Operand> const & operand, Scalar const & scalar) {
        return vectorize(detail::NotEqualToTag::template ExprScalar<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::NotEqualToTag::template ScalarExpr<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator !=(Scalar const & scalar, Expression<Operand> const & operand) {
        return vectorize(detail::NotEqualToTag::template ScalarExpr<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand1, typename Operand2>
#ifndef DOXYGEN
    detail::BinaryOpExpression< 
         Operand1, Operand2,
         typename detail::NotEqualToTag::template ExprExpr<Operand1,Operand2>::BinaryFunction
    >
#else
    <unspecified-expression-type>
#endif
    operator !=(Expression<Operand1> const & operand1, Expression<Operand2> const & operand2) {
        return vectorize(
            typename detail::NotEqualToTag::template ExprExpr<Operand1,Operand2>::BinaryFunction(),
            operand1,
            operand2
        );
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::LessTag::template ExprScalar<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator <(Expression<Operand> const & operand, Scalar const & scalar) {
        return vectorize(detail::LessTag::template ExprScalar<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::LessTag::template ScalarExpr<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator <(Scalar const & scalar, Expression<Operand> const & operand) {
        return vectorize(detail::LessTag::template ScalarExpr<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand1, typename Operand2>
#ifndef DOXYGEN
    detail::BinaryOpExpression< 
         Operand1, Operand2,
         typename detail::LessTag::template ExprExpr<Operand1,Operand2>::BinaryFunction
    >
#else
    <unspecified-expression-type>
#endif
    operator <(Expression<Operand1> const & operand1, Expression<Operand2> const & operand2) {
        return vectorize(
            typename detail::LessTag::template ExprExpr<Operand1,Operand2>::BinaryFunction(),
            operand1,
            operand2
        );
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::GreaterTag::template ExprScalar<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator >(Expression<Operand> const & operand, Scalar const & scalar) {
        return vectorize(detail::GreaterTag::template ExprScalar<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::GreaterTag::template ScalarExpr<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator >(Scalar const & scalar, Expression<Operand> const & operand) {
        return vectorize(detail::GreaterTag::template ScalarExpr<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand1, typename Operand2>
#ifndef DOXYGEN
    detail::BinaryOpExpression< 
         Operand1, Operand2,
         typename detail::GreaterTag::template ExprExpr<Operand1,Operand2>::BinaryFunction
    >
#else
    <unspecified-expression-type>
#endif
    operator >(Expression<Operand1> const & operand1, Expression<Operand2> const & operand2) {
        return vectorize(
            typename detail::GreaterTag::template ExprExpr<Operand1,Operand2>::BinaryFunction(),
            operand1,
            operand2
        );
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::LessEqualTag::template ExprScalar<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator <=(Expression<Operand> const & operand, Scalar const & scalar) {
        return vectorize(detail::LessEqualTag::template ExprScalar<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::LessEqualTag::template ScalarExpr<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator <=(Scalar const & scalar, Expression<Operand> const & operand) {
        return vectorize(detail::LessEqualTag::template ScalarExpr<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand1, typename Operand2>
#ifndef DOXYGEN
    detail::BinaryOpExpression< 
         Operand1, Operand2,
         typename detail::LessEqualTag::template ExprExpr<Operand1,Operand2>::BinaryFunction
    >
#else
    <unspecified-expression-type>
#endif
    operator <=(Expression<Operand1> const & operand1, Expression<Operand2> const & operand2) {
        return vectorize(
            typename detail::LessEqualTag::template ExprExpr<Operand1,Operand2>::BinaryFunction(),
            operand1,
            operand2
        );
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::GreaterEqualTag::template ExprScalar<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator >=(Expression<Operand> const & operand, Scalar const & scalar) {
        return vectorize(detail::GreaterEqualTag::template ExprScalar<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::GreaterEqualTag::template ScalarExpr<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator >=(Scalar const & scalar, Expression<Operand> const & operand) {
        return vectorize(detail::GreaterEqualTag::template ScalarExpr<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand1, typename Operand2>
#ifndef DOXYGEN
    detail::BinaryOpExpression< 
         Operand1, Operand2,
         typename detail::GreaterEqualTag::template ExprExpr<Operand1,Operand2>::BinaryFunction
    >
#else
    <unspecified-expression-type>
#endif
    operator >=(Expression<Operand1> const & operand1, Expression<Operand2> const & operand2) {
        return vectorize(
            typename detail::GreaterEqualTag::template ExprExpr<Operand1,Operand2>::BinaryFunction(),
            operand1,
            operand2
        );
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::LogicalAnd::template ExprScalar<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator &&(Expression<Operand> const & operand, Scalar const & scalar) {
        return vectorize(detail::LogicalAnd::template ExprScalar<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::LogicalAnd::template ScalarExpr<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator &&(Scalar const & scalar, Expression<Operand> const & operand) {
        return vectorize(detail::LogicalAnd::template ScalarExpr<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand1, typename Operand2>
#ifndef DOXYGEN
    detail::BinaryOpExpression< 
         Operand1, Operand2,
         typename detail::LogicalAnd::template ExprExpr<Operand1,Operand2>::BinaryFunction
    >
#else
    <unspecified-expression-type>
#endif
    operator &&(Expression<Operand1> const & operand1, Expression<Operand2> const & operand2) {
        return vectorize(
            typename detail::LogicalAnd::template ExprExpr<Operand1,Operand2>::BinaryFunction(),
            operand1,
            operand2
        );
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::LogicalOr::template ExprScalar<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator ||(Expression<Operand> const & operand, Scalar const & scalar) {
        return vectorize(detail::LogicalOr::template ExprScalar<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand, typename Scalar>
#ifndef DOXYGEN
    typename boost::enable_if<
        typename detail::ExpressionTraits<Scalar>::IsScalar,
        detail::UnaryOpExpression< Operand, typename detail::LogicalOr::template ScalarExpr<Operand,Scalar>::Bound >
    >::type
#else
    <unspecified-expression-type>
#endif
    operator ||(Scalar const & scalar, Expression<Operand> const & operand) {
        return vectorize(detail::LogicalOr::template ScalarExpr<Operand,Scalar>::bind(scalar),operand);
    }

    template <typename Operand1, typename Operand2>
#ifndef DOXYGEN
    detail::BinaryOpExpression< 
         Operand1, Operand2,
         typename detail::LogicalOr::template ExprExpr<Operand1,Operand2>::BinaryFunction
    >
#else
    <unspecified-expression-type>
#endif
    operator ||(Expression<Operand1> const & operand1, Expression<Operand2> const & operand2) {
        return vectorize(
            typename detail::LogicalOr::template ExprExpr<Operand1,Operand2>::BinaryFunction(),
            operand1,
            operand2
        );
    }


    template <typename Operand>
#ifndef DOXYGEN
    detail::UnaryOpExpression< Operand, std::negate<typename detail::ExpressionTraits<Operand>::Element> >
#else
    <unspecified-expression-type>
#endif
    operator -(Expression<Operand> const & operand) {
        return vectorize(std::negate<typename detail::ExpressionTraits<Operand>::Element>(),operand);
    }

    template <typename Operand>
#ifndef DOXYGEN
    detail::UnaryOpExpression< Operand, std::logical_not<typename detail::ExpressionTraits<Operand>::Element> >
#else
    <unspecified-expression-type>
#endif
    operator !(Expression<Operand> const & operand) {
        return vectorize(std::logical_not<typename detail::ExpressionTraits<Operand>::Element>(),operand);
    }

    template <typename Operand>
#ifndef DOXYGEN
    detail::UnaryOpExpression< Operand, detail::BitwiseNot<typename detail::ExpressionTraits<Operand>::Element> >
#else
    <unspecified-expression-type>
#endif
    operator ~(Expression<Operand> const & operand) {
        return vectorize(detail::BitwiseNot<typename detail::ExpressionTraits<Operand>::Element>(),operand);
    }
/// @}

template <typename Scalar>
inline typename boost::enable_if<typename detail::ExpressionTraits<Scalar>::IsScalar, bool>::type
any(Scalar const & scalar) {
    return bool(scalar);
}

/**
 *  \brief Return true if any of the elements of the given expression are true.
 *
 *  \ingroup MainGroup
 */
template <typename Derived>
inline bool
any(Expression<Derived> const & expr) {
    typename Derived::Iterator const i_end = expr.end();
    for (typename Derived::Iterator i = expr.begin(); i != i_end; ++i) {
        if (any(*i)) return true;
    }
    return false;
}

template <typename Scalar>
inline typename boost::enable_if<typename detail::ExpressionTraits<Scalar>::IsScalar, bool>::type
all(Scalar const & scalar) {
    return bool(scalar);
}

/**
 *  \brief Return true if all of the elements of the given expression are true.
 *
 *  \ingroup MainGroup
 */
template <typename Derived>
inline bool
all(Expression<Derived> const & expr) {
    typename Derived::Iterator const i_end = expr.end();
    for (typename Derived::Iterator i = expr.begin(); i != i_end; ++i) {
        if (!all(*i)) return false;
    }
    return true;
}

template <typename Scalar>
inline typename boost::enable_if<typename detail::ExpressionTraits<Scalar>::IsScalar, Scalar>::type
sum(Scalar const & scalar) { return scalar; }


/**
 *  \brief Return the sum of all elements of the given expression.
 *
 *  \ingroup MainGroup
 */
template <typename Derived>
inline typename Derived::Element
sum(Expression<Derived> const & expr) {
    typename Derived::Iterator const i_end = expr.end();
    typename Derived::Element total = static_cast<typename Derived::Element>(0);
    for (typename Derived::Iterator i = expr.begin(); i != i_end; ++i) {
        total += sum(*i);
    }
    return total;
}


} // namespace ndarray

#endif // !NDARRAY_operators_hpp_INCLUDED