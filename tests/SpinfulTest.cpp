//
// Created by Andy on 6/17/2026.
//

// ReSharper disable CppTemplateArgumentsCanBeDeduced
#include <cmath>
#include <complex>
#include <functional>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/SparseCore>

#include <SecChem/Spinful.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

namespace
{
	using Catch::Approx;

	using Eigen::MatrixXd;
	using Eigen::VectorXd;

	// Expression/view zoo of spec §5.1 (types are only probed, never named exactly).
	using SumExpr = decltype(std::declval<const MatrixXd&>() + std::declval<const MatrixXd&>());
	using ProductExpr = decltype(std::declval<const MatrixXd&>() * std::declval<const MatrixXd&>());
	using ScaledExpr = decltype(2.0 * std::declval<const MatrixXd&>());
	using BlockExpr = decltype(std::declval<MatrixXd&>().block(0, 0, 1, 1));
	using DiagonalExpr = decltype(std::declval<MatrixXd&>().diagonal());
	using AsDiagonalExpr = decltype(std::declval<const VectorXd&>().asDiagonal());
	using ConstMatrixXd = const MatrixXd;

	struct NotEigen
	{
		double Value;
	};

	// Mock foreign component family (spec §5.3): one opt-in declaration for the whole
	// family via the expression base class. Values are plain doubles so expectations
	// stay exact. `+` is the family's lazy operator (ForeignSum nodes nest by value);
	// `-` and `*` are eager owner results; compound `+=` accepts any family member.
	struct ExprBase
	{
		static constexpr bool NestByValue = true;
	};

	struct ForeignMatrix : ExprBase
	{
		using Scalar = double;
		using PlainObject = ForeignMatrix;
		static constexpr bool NestByValue = false;

		ForeignMatrix() = default;
		explicit ForeignMatrix(const double value) : m_Value{value}
		{
		}

		double Value() const
		{
			return m_Value;
		}

		template <typename TRhs, typename = std::enable_if_t<std::is_base_of_v<ExprBase, TRhs>>>
		ForeignMatrix& operator+=(const TRhs& rhs)
		{
			m_Value += rhs.Value();
			return *this;
		}

		double m_Value = 0.0;
	};

	template <typename TLhs, typename TRhs>
	struct ForeignSum : ExprBase
	{
		using Scalar = std::common_type_t<typename TLhs::Scalar, typename TRhs::Scalar>;
		using PlainObject = typename TLhs::PlainObject;

		ForeignSum(const TLhs& lhs, const TRhs& rhs) : m_Lhs{lhs}, m_Rhs{rhs}
		{
		}

		double Value() const
		{
			return m_Lhs.Value() + m_Rhs.Value();
		}

		// Implicit conversion to the family owner — Spinful's block-assignment paths
		// evaluate through this.
		operator typename TLhs::PlainObject() const
		{
			return typename TLhs::PlainObject{Value()};
		}

		TLhs m_Lhs;
		TRhs m_Rhs;
	};

	// Family-generic sum: any two members (owner or node) nest into a ForeignSum.
	template <typename TLhs,
	          typename TRhs,
	          typename = std::enable_if_t<std::is_base_of_v<ExprBase, TLhs> && std::is_base_of_v<ExprBase, TRhs>>>
	ForeignSum<TLhs, TRhs> operator+(const TLhs& lhs, const TRhs& rhs)
	{
		return ForeignSum<TLhs, TRhs>{lhs, rhs};
	}

	inline ForeignMatrix operator-(const ForeignMatrix& lhs, const ForeignMatrix& rhs)
	{
		return ForeignMatrix{lhs.Value() - rhs.Value()};
	}

	inline ForeignMatrix operator*(const ForeignMatrix& lhs, const ForeignMatrix& rhs)
	{
		return ForeignMatrix{lhs.Value() * rhs.Value()};
	}

	inline ForeignMatrix operator*(const double scalar, const ForeignMatrix& matrix)
	{
		return ForeignMatrix{scalar * matrix.Value()};
	}

	inline ForeignMatrix operator*(const ForeignMatrix& matrix, const double scalar)
	{
		return ForeignMatrix{matrix.Value() * scalar};
	}

	// A family whose base carries no NestByValue policy at all.
	struct PolicylessBase
	{
	};

	struct PolicylessNode : PolicylessBase
	{
		using Scalar = float;
		using PlainObject = PolicylessNode;
	};

	// Upcasts to the CRTP base: with a base-typed argument the plain copy/move
	// overloads are not viable, so the expression constructor / expression assignment
	// are selected — the way to exercise the assignment lattice owner→owner before
	// expression nodes exist (P2).
	template <typename E>
	const SecChem::SpinfulExpr<E>& AsExpression(const SecChem::SpinfulExpr<E>& expression)
	{
		return expression;
	}

	// Per-block value check for any spinful expression — owner or node: logical Alpha
	// and Beta against manual references.
	template <typename E>
	void CheckBlocks(const SecChem::SpinfulExpr<E>& expression,
	                 const MatrixXd& expectedAlpha,
	                 const MatrixXd& expectedBeta)
	{
		CHECK((expression[SecChem::Spin::Alpha].array() == expectedAlpha.array()).all());
		CHECK((expression[SecChem::Spin::Beta].array() == expectedBeta.array()).all());
	}

	// Declaration-only SFINAE probes for the free operators (spec §9 test 6): C++17
	// forbids lambdas in unevaluated contexts, so is_invocable drives named functors.
	struct ProbeAdd
	{
		template <typename A, typename B>
		auto operator()(const A& a, const B& b) const -> decltype(a + b);
	};

	struct ProbeSub
	{
		template <typename A, typename B>
		auto operator()(const A& a, const B& b) const -> decltype(a - b);
	};

	struct ProbeMul
	{
		template <typename A, typename B>
		auto operator()(const A& a, const B& b) const -> decltype(a * b);
	};

	struct ProbeDiv
	{
		template <typename A, typename B>
		auto operator()(const A& a, const B& b) const -> decltype(a / b);
	};

	struct ProbeCompoundAdd
	{
		template <typename A, typename B>
		auto operator()(A& a, const B& b) const -> decltype(a += b);
	};

	struct ProbeCompoundMul
	{
		template <typename A, typename B>
		auto operator()(A& a, const B& b) const -> decltype(a *= b);
	};

	// Inter operands are SFINAE-constrained to real components (spec §6.3): invalid
	// operand types have no matching operator...
	static_assert(!std::is_invocable_v<ProbeAdd, const SecChem::Spinful<MatrixXd>&, const std::string&>);
	static_assert(!std::is_invocable_v<ProbeAdd, const std::string&, const SecChem::Spinful<MatrixXd>&>);
	static_assert(!std::is_invocable_v<ProbeSub, const SecChem::Spinful<MatrixXd>&, const std::vector<double>&>);
	static_assert(!std::is_invocable_v<ProbeMul, const std::vector<double>&, const SecChem::Spinful<MatrixXd>&>);
	static_assert(!std::is_invocable_v<ProbeDiv, const SecChem::Spinful<MatrixXd>&, const std::string&>);

	// ... valid operands do match in both directions (a scalar operand constructs the
	// inter node for any component kind — evaluating matrix ± scalar would hit Eigen's
	// missing operators, availability-by-call, spec §6.5) ...
	// static_assert(!std::is_invocable_v<ProbeAdd, const SecChem::Spinful<MatrixXd>&, const double&>);
	// static_assert(!std::is_invocable_v<ProbeAdd, const double&, const SecChem::Spinful<MatrixXd>&>);
	static_assert(std::is_invocable_v<ProbeAdd, const MatrixXd&, const SecChem::Spinful<MatrixXd>&>);
	static_assert(std::is_invocable_v<ProbeAdd, const SecChem::Spinful<MatrixXd>&, const MatrixXd&>);
	static_assert(std::is_invocable_v<ProbeSub, const SecChem::Spinful<MatrixXd>&, const MatrixXd&>);
	static_assert(std::is_invocable_v<ProbeSub, const MatrixXd&, const SecChem::Spinful<MatrixXd>&>);
	static_assert(std::is_invocable_v<ProbeMul, const SecChem::Spinful<MatrixXd>&, const double&>);
	static_assert(std::is_invocable_v<ProbeMul, const double&, const SecChem::Spinful<MatrixXd>&>);
	static_assert(std::is_invocable_v<ProbeDiv, const SecChem::Spinful<MatrixXd>&, const double&>);
	// static_assert(!std::is_invocable_v<ProbeDiv, const double&, const SecChem::Spinful<MatrixXd>&>);
	static_assert(std::is_invocable_v<ProbeDiv, const double&, const SecChem::Spinful<double>&>);

	// ... and spinful ⊕ spinful stays on the intra path.
	static_assert(std::is_invocable_v<ProbeAdd, const SecChem::Spinful<MatrixXd>&, const SecChem::Spinful<MatrixXd>&>);

	// Compound scalar assignment follows the same doctrine: matrix ± scalar is
	// meaningless, no scalar operator+=/-= overloads exist — but scaling stays.
	static_assert(!std::is_invocable_v<ProbeCompoundAdd, SecChem::Spinful<MatrixXd>&, const double&>);
	static_assert(std::is_invocable_v<ProbeCompoundMul, SecChem::Spinful<MatrixXd>&, const double&>);
	static_assert(std::is_invocable_v<ProbeCompoundAdd, SecChem::Spinful<double>&, const double&>);

	// EigenBase derivation (the Trace adapter's scalar guard, spec §6.6): every Eigen
	// object — owner, view, or expression, sparse included — derives from
	// Eigen::EigenBase; scalars and foreign families never do.
	static_assert(SecChem::Detail::IsEigenType<MatrixXd>);
	static_assert(SecChem::Detail::IsEigenType<Eigen::SparseMatrix<double>>);
	static_assert(SecChem::Detail::IsEigenType<Eigen::Map<MatrixXd>>);
	static_assert(SecChem::Detail::IsEigenType<SumExpr>);
	static_assert(!SecChem::Detail::IsEigenType<double>);
	static_assert(!SecChem::Detail::IsEigenType<std::complex<double>>);
	static_assert(!SecChem::Detail::IsEigenType<Int32>);

	// The Trace adapter dispatches by component family (spec §6.6): Eigen components
	// route to trace(), scalars and complex to the identity — and an Eigen component
	// without trace() (the sparse types) finds no viable overload, so Trace() is
	// ill-formed at its call site.
	static_assert(std::is_invocable_v<SecChem::Detail::TraceOf, const double&, int>);
	static_assert(std::is_invocable_v<SecChem::Detail::TraceOf, const std::complex<double>&, int>);
	static_assert(std::is_invocable_v<SecChem::Detail::TraceOf, const MatrixXd&, int>);
	static_assert(!std::is_invocable_v<SecChem::Detail::TraceOf, const Eigen::SparseMatrix<double>&, int>);

	// Shipped scalar specializations are trivial: Scalar = Plain = T, NestByValue = true.
	template <typename T>
	constexpr bool HasTrivialScalarTraits = std::is_same_v<typename SecChem::Detail::ComponentTraits<T>::Scalar, T>
	                                        && std::is_same_v<typename SecChem::Detail::ComponentTraits<T>::Plain, T>
	                                        && SecChem::Detail::ComponentTraits<T>::NestByValue;

	// Detection table of spec §5.1: everything Eigen-ish is detected by the probe...
	static_assert(SecChem::Detail::IsEigenType<MatrixXd>);
	static_assert(SecChem::Detail::IsEigenType<VectorXd>);
	static_assert(SecChem::Detail::IsEigenType<Eigen::ArrayXXd>);
	static_assert(SecChem::Detail::IsEigenType<Eigen::Transpose<MatrixXd>>);
	static_assert(SecChem::Detail::IsEigenType<BlockExpr>);
	static_assert(SecChem::Detail::IsEigenType<DiagonalExpr>);
	static_assert(SecChem::Detail::IsEigenType<Eigen::Map<MatrixXd>>);
	static_assert(SecChem::Detail::IsEigenType<Eigen::SparseMatrix<double>>);
	static_assert(SecChem::Detail::IsEigenType<ConstMatrixXd>);
	static_assert(SecChem::Detail::IsEigenType<SumExpr>);
	static_assert(SecChem::Detail::IsEigenType<ProductExpr>);
	static_assert(SecChem::Detail::IsEigenType<ScaledExpr>);
	static_assert(SecChem::Detail::IsEigenType<AsDiagonalExpr>);

	// ... and nothing else is.
	static_assert(!SecChem::Detail::IsEigenType<double>);
	static_assert(!SecChem::Detail::IsEigenType<float>);
	static_assert(!SecChem::Detail::IsEigenType<std::complex<double>>);
	static_assert(!SecChem::Detail::IsEigenType<Int32>);
	static_assert(!SecChem::Detail::IsEigenType<NotEigen>);
	static_assert(!SecChem::Detail::IsEigenType<ForeignMatrix>);

	// NestByValue reproduces Eigen's own nesting policy: owners false...
	static_assert(!SecChem::Detail::ComponentTraits<MatrixXd>::NestByValue);
	static_assert(!SecChem::Detail::ComponentTraits<VectorXd>::NestByValue);
	static_assert(!SecChem::Detail::ComponentTraits<Eigen::ArrayXXd>::NestByValue);
	static_assert(!SecChem::Detail::ComponentTraits<Eigen::SparseMatrix<double>>::NestByValue);
	static_assert(!SecChem::Detail::ComponentTraits<ConstMatrixXd>::NestByValue);

	// ... expressions/views true.
	static_assert(SecChem::Detail::ComponentTraits<Eigen::Transpose<MatrixXd>>::NestByValue);
	static_assert(SecChem::Detail::ComponentTraits<BlockExpr>::NestByValue);
	static_assert(SecChem::Detail::ComponentTraits<DiagonalExpr>::NestByValue);
	static_assert(SecChem::Detail::ComponentTraits<Eigen::Map<MatrixXd>>::NestByValue);
	static_assert(SecChem::Detail::ComponentTraits<SumExpr>::NestByValue);
	static_assert(SecChem::Detail::ComponentTraits<ProductExpr>::NestByValue);
	static_assert(SecChem::Detail::ComponentTraits<ScaledExpr>::NestByValue);
	static_assert(SecChem::Detail::ComponentTraits<AsDiagonalExpr>::NestByValue);

	// Plain types: materialized matrix for owners and expressions...
	static_assert(std::is_same_v<SecChem::Detail::ComponentTraits<MatrixXd>::Plain, MatrixXd>);
	static_assert(std::is_same_v<SecChem::Detail::ComponentTraits<ConstMatrixXd>::Plain, MatrixXd>);
	static_assert(std::is_same_v<SecChem::Detail::ComponentTraits<SumExpr>::Plain, MatrixXd>);
	static_assert(std::is_same_v<SecChem::Detail::ComponentTraits<ProductExpr>::Plain, MatrixXd>);

	// ... sparse owners bridge through plain_matrix_type's sparse specialization,
	// preserving the storage order ...
	static_assert(std::is_same_v<SecChem::Detail::ComponentTraits<Eigen::SparseMatrix<double>>::Plain,
	                             Eigen::SparseMatrix<double>>);
	static_assert(std::is_same_v<SecChem::Detail::ComponentTraits<Eigen::SparseMatrix<double, Eigen::RowMajor>>::Plain,
	                             Eigen::SparseMatrix<double, Eigen::RowMajor>>);

	// ... and mapped sparse views carry the owning Plain with by-value nesting — the
	// mapping is what copies, the same policy as dense views.
	static_assert(std::is_same_v<SecChem::Detail::ComponentTraits<Eigen::Map<Eigen::SparseMatrix<double>>>::Plain,
	                             Eigen::SparseMatrix<double>>);
	static_assert(SecChem::Detail::ComponentTraits<Eigen::Map<Eigen::SparseMatrix<double>>>::NestByValue);
	static_assert(std::is_same_v<SecChem::Detail::ComponentTraits<Eigen::Map<const Eigen::SparseMatrix<double>>>::Plain,
	                             Eigen::SparseMatrix<double>>);
	static_assert(SecChem::Detail::ComponentTraits<Eigen::Map<const Eigen::SparseMatrix<double>>>::NestByValue);

	// ... the honest diagonal type for asDiagonal (assigns to MatrixXd when materialized)...
	static_assert(std::is_same_v<SecChem::Detail::ComponentTraits<AsDiagonalExpr>::Plain,
	                             Eigen::DiagonalMatrix<double, Eigen::Dynamic>>);

	// ... and the scalars go through the shipped specializations.
	static_assert(HasTrivialScalarTraits<Int8>);
	static_assert(HasTrivialScalarTraits<Int16>);
	static_assert(HasTrivialScalarTraits<Int32>);
	static_assert(HasTrivialScalarTraits<Int64>);
	static_assert(HasTrivialScalarTraits<float>);
	static_assert(HasTrivialScalarTraits<double>);
	static_assert(HasTrivialScalarTraits<std::complex<float>>);
	static_assert(HasTrivialScalarTraits<std::complex<double>>);

	// Layer 2 resolves through inheritance: an owner's own member shadows the family
	// default; nested nodes inherit it at arbitrary depth; families without a policy
	// get the true default.
	static_assert(std::is_same_v<SecChem::Detail::ComponentTraits<ForeignMatrix>::Scalar, double>);
	static_assert(std::is_same_v<SecChem::Detail::ComponentTraits<ForeignMatrix>::Plain, ForeignMatrix>);
	static_assert(!SecChem::Detail::ComponentTraits<ForeignMatrix>::NestByValue);
	static_assert(SecChem::Detail::ComponentTraits<ForeignSum<ForeignMatrix, ForeignMatrix>>::NestByValue);
	static_assert(std::is_same_v<SecChem::Detail::ComponentTraits<ForeignSum<ForeignMatrix, ForeignMatrix>>::Plain,
	                             ForeignMatrix>);
	static_assert(SecChem::Detail::ComponentTraits<
	              ForeignSum<ForeignMatrix, ForeignSum<ForeignMatrix, ForeignMatrix>>>::NestByValue);
	static_assert(SecChem::Detail::ComponentTraits<PolicylessNode>::NestByValue);
}  // namespace

// The single-component constructor is explicit: no implicit raw-component conversion.
static_assert(!std::is_convertible_v<const MatrixXd&, SecChem::Spinful<MatrixXd>>);

TEST_CASE("Spinful construction")
{
	SECTION("default construction is Empty with public count 1")
	{
		const SecChem::Spinful<MatrixXd> s;

		CHECK(s.IsEmpty());
		CHECK_FALSE(s.IsRestricted());
		CHECK(s.IndependentComponentCount() == 1);
		CHECK(s.Alpha().rows() == 0);
		CHECK(s.Alpha().cols() == 0);
	}

	SECTION("single-component construction is restricted")
	{
		const MatrixXd a = MatrixXd::Constant(2, 3, 1.0);
		const SecChem::Spinful<MatrixXd> s(a);

		CHECK_FALSE(s.IsEmpty());
		CHECK(s.IsRestricted());
		CHECK(s.IndependentComponentCount() == 1);
		CHECK((s.Alpha().array() == a.array()).all());
	}

	SECTION("two-component construction is unrestricted")
	{
		const MatrixXd a = MatrixXd::Constant(2, 2, 1.0);
		const MatrixXd b = MatrixXd::Constant(2, 2, 2.0);
		const SecChem::Spinful<MatrixXd> s{a, b};

		CHECK_FALSE(s.IsEmpty());
		CHECK_FALSE(s.IsRestricted());
		CHECK(s.IndependentComponentCount() == 2);
		CHECK((s.Alpha().array() == a.array()).all());
		CHECK((s.Beta().array() == b.array()).all());
	}

	SECTION("two-component construction: pairs and the replicated spelling")
	{
		const MatrixXd a = MatrixXd::Constant(2, 2, 3.0);

		const SecChem::Spinful<MatrixXd> unrestricted{a, a};
		CHECK(unrestricted.IndependentComponentCount() == 2);
		CHECK_FALSE(unrestricted.IsRestricted());
		CHECK((unrestricted.Alpha().array() == a.array()).all());
		CHECK((unrestricted.Beta().array() == a.array()).all());

		const SecChem::Spinful<MatrixXd> restricted(a);
		CHECK(restricted.IndependentComponentCount() == 1);
		CHECK(restricted.IsRestricted());

		// braces and parentheses agree on the two-component form (D28): the int
		// converts to the component scalar, so both spellings build the pair (0.5, 2.0)
		const SecChem::Spinful<double> parenPair(0.5, 2);
		CHECK(parenPair.IndependentComponentCount() == 2);
		CHECK(parenPair.Alpha() == 0.5);
		CHECK(parenPair.Beta() == 2.0);

		const SecChem::Spinful<double> valueList{0.5, 2};
		CHECK(valueList.IndependentComponentCount() == 2);
		CHECK(valueList.Alpha() == 0.5);
		CHECK(valueList.Beta() == 2.0);

		// the pair constructor stays implicit: the assignment brace form keeps compiling
		SecChem::Spinful<double> assigned;
		assigned = {1.0, 2.0};
		CHECK(assigned.Beta() == 2.0);
	}

	SECTION("move-in construction moves each component")
	{
		MatrixXd a = MatrixXd::Constant(2, 2, 1.0);
		MatrixXd b = MatrixXd::Constant(2, 2, 2.0);
		const double* const alphaData = a.data();
		const double* const betaData = b.data();

		const SecChem::Spinful<MatrixXd> s{std::move(a), std::move(b)};
		CHECK(s.IndependentComponentCount() == 2);
		CHECK(s.Alpha().data() == alphaData);
		CHECK(s.Beta().data() == betaData);
		CHECK(a.size() == 0);
		CHECK(b.size() == 0);

		// mixed argument categories: the lvalue copies, the rvalue moves
		MatrixXd c = MatrixXd::Constant(2, 2, 3.0);
		MatrixXd d = MatrixXd::Constant(2, 2, 4.0);
		const double* const betaDataOfD = d.data();
		const SecChem::Spinful<MatrixXd> mixed{c, std::move(d)};
		CHECK(mixed.Alpha().data() != c.data());
		CHECK((mixed.Alpha().array() == c.array()).all());
		CHECK(mixed.Beta().data() == betaDataOfD);
		CHECK(d.size() == 0);
	}

	SECTION("move preserves the count; the moved-from object keeps its count")
	{
		const MatrixXd a = MatrixXd::Constant(2, 2, 1.0);
		const MatrixXd b = MatrixXd::Constant(2, 2, 2.0);
		SecChem::Spinful<MatrixXd> s{a, b};

		const SecChem::Spinful<MatrixXd> moved(std::move(s));
		CHECK(moved.IndependentComponentCount() == 2);
		CHECK(s.IndependentComponentCount() == 2);
		CHECK_FALSE(s.IsEmpty());
	}

	SECTION("scalar default construction holds the value component")
	{
		const SecChem::Spinful<double> s;

		CHECK(s.IsEmpty());
		CHECK(s.Alpha() == 0.0);
	}

	SECTION("matrix factories never produce Empty")
	{
		const auto zero = SecChem::Spinful<MatrixXd>::Zero(2, 3);
		CHECK(zero.IsRestricted());
		CHECK(zero.Alpha().isZero());

		const auto zero2 = SecChem::Spinful<MatrixXd>::Zero(2, 3, 2);
		CHECK(zero2.IndependentComponentCount() == 2);
		CHECK(zero2.Alpha().isZero());
		CHECK(zero2.Beta().isZero());

		const auto ones = SecChem::Spinful<MatrixXd>::Ones(2, 3);
		CHECK(ones.IsRestricted());
		CHECK(ones.Alpha().isOnes());

		const auto ones2 = SecChem::Spinful<MatrixXd>::Ones(2, 3, 2);
		CHECK(ones2.IndependentComponentCount() == 2);
		CHECK(ones2.Alpha().isOnes());
		CHECK(ones2.Beta().isOnes());

		const auto identity = SecChem::Spinful<MatrixXd>::Identity(3);
		CHECK(identity.IsRestricted());
		CHECK(identity.Alpha().isIdentity());

		const auto identity2 = SecChem::Spinful<MatrixXd>::Identity(3, 2);
		CHECK(identity2.IndependentComponentCount() == 2);
		CHECK(identity2.Alpha().isIdentity());
		CHECK(identity2.Beta().isIdentity());
	}

	SECTION("scalar factories never produce Empty")
	{
		const auto zero = SecChem::Spinful<double>::Zero();
		CHECK(zero.IsRestricted());
		CHECK(zero.Alpha() == 0.0);

		const auto zero2 = SecChem::Spinful<double>::Zero(2);
		CHECK(zero2.IndependentComponentCount() == 2);
		CHECK(zero2.Alpha() == 0.0);
		CHECK(zero2.Beta() == 0.0);

		const auto one = SecChem::Spinful<double>::One();
		CHECK(one.IsRestricted());
		CHECK(one.Alpha() == 1.0);

		const auto one2 = SecChem::Spinful<double>::One(2);
		CHECK(one2.IndependentComponentCount() == 2);
		CHECK(one2.Alpha() == 1.0);
		CHECK(one2.Beta() == 1.0);
	}
}

TEST_CASE("Spinful broadcast reads")
{
	const MatrixXd a = MatrixXd::Constant(2, 2, 1.0);
	const MatrixXd b = MatrixXd::Constant(2, 2, 2.0);

	SECTION("restricted: Beta is the same object as Alpha")
	{
		SecChem::Spinful<MatrixXd> s(a);

		CHECK(&s[SecChem::Spin::Beta] == &s[SecChem::Spin::Alpha]);
		CHECK(&s.Beta() == &s.Alpha());

		const SecChem::Spinful<MatrixXd>& cs = s;
		CHECK(&cs[SecChem::Spin::Beta] == &cs[SecChem::Spin::Alpha]);
		CHECK(&cs.Beta() == &cs.Alpha());
	}

	SECTION("unrestricted: the blocks are distinct objects with their own values")
	{
		const SecChem::Spinful<MatrixXd> s{a, b};

		CHECK(&s[SecChem::Spin::Beta] != &s[SecChem::Spin::Alpha]);
		CHECK(&s.Beta() != &s.Alpha());
		CHECK((s[SecChem::Spin::Alpha].array() == a.array()).all());
		CHECK((s[SecChem::Spin::Beta].array() == b.array()).all());
	}

	SECTION("Empty broadcasts like restricted over default-constructed components")
	{
		const SecChem::Spinful<MatrixXd> s;

		CHECK(&s[SecChem::Spin::Beta] == &s[SecChem::Spin::Alpha]);
		CHECK(s.Alpha().rows() == 0);
	}
}

TEST_CASE("Spinful broadcast writes")
{
	const MatrixXd a = MatrixXd::Constant(2, 2, 1.0);
	const MatrixXd m = MatrixXd::Constant(2, 2, 7.0);

	SECTION("restricted: writing through [Beta] hits the shared block; count unchanged")
	{
		SecChem::Spinful<MatrixXd> s(a);
		s[SecChem::Spin::Beta] = m;

		CHECK((s.Alpha().array() == m.array()).all());
		CHECK((s.Beta().array() == m.array()).all());
		CHECK(s.IsRestricted());
		CHECK(s.IndependentComponentCount() == 1);
	}

	SECTION("restricted: writing through Beta() hits the shared block")
	{
		SecChem::Spinful<MatrixXd> s(a);
		s.Beta() = m;

		CHECK((s.Alpha().array() == m.array()).all());
		CHECK(s.IsRestricted());
	}

	SECTION("unrestricted: writing through [Beta] hits only the Beta block")
	{
		SecChem::Spinful<MatrixXd> s{a, a};
		s[SecChem::Spin::Beta] = m;

		CHECK((s.Alpha().array() == a.array()).all());
		CHECK((s.Beta().array() == m.array()).all());
		CHECK(s.IndependentComponentCount() == 2);
	}

	SECTION("count-generic assignment loop over AllSpins is legal and stays restricted")
	{
		SecChem::Spinful<MatrixXd> s(a);
		for (const SecChem::Spin spin : SecChem::AllSpins)
		{
			s[spin] = m;
		}

		// Both iterations write the one shared block.
		CHECK((s.Alpha().array() == m.array()).all());
		CHECK(s.IsRestricted());
	}

	SECTION("count-generic mutation loop over AllSpins doubles the shared block twice")
	{
		SecChem::Spinful<MatrixXd> s(a);
		for (const SecChem::Spin spin : SecChem::AllSpins)
		{
			s[spin] *= 2;
		}

		// Alpha and Beta address the same storage, so *= applies twice: 1 * 2 * 2.
		CHECK((s.Alpha().array() == Eigen::MatrixXd::Constant(2, 2, 4.0).array()).all());
		CHECK(s.IsRestricted());
	}

	SECTION("non-const write through [spin] exits Empty to restricted")
	{
		SecChem::Spinful<MatrixXd> s;
		s[SecChem::Spin::Alpha] = m;

		CHECK_FALSE(s.IsEmpty());
		CHECK(s.IsRestricted());
		CHECK(s.IndependentComponentCount() == 1);
		CHECK((s[SecChem::Spin::Beta].array() == m.array()).all());
	}

	SECTION("even a non-const read through [spin] exits Empty (harmless)")
	{
		SecChem::Spinful<MatrixXd> s;
		(void)s[SecChem::Spin::Alpha];

		CHECK_FALSE(s.IsEmpty());
		CHECK(s.IsRestricted());
	}

	SECTION("non-const Alpha() exits Empty; the write lands in the shared block")
	{
		SecChem::Spinful<MatrixXd> s;
		s.Alpha() = m;

		CHECK_FALSE(s.IsEmpty());
		CHECK(s.IsRestricted());
		CHECK((s.Alpha().array() == m.array()).all());
		CHECK((s.Beta().array() == m.array()).all());
	}

	SECTION("non-const Beta() exits Empty and aliases the Alpha block")
	{
		SecChem::Spinful<MatrixXd> s;
		s.Beta() = m;

		CHECK_FALSE(s.IsEmpty());
		CHECK(s.IsRestricted());
		CHECK(&s.Beta() == &s.Alpha());
	}

	SECTION("non-const AtIndex exits Empty; index 1 broadcasts to slot 0")
	{
		SecChem::Spinful<MatrixXd> s;
		s.AtIndex(1) = m;

		CHECK_FALSE(s.IsEmpty());
		CHECK(s.IsRestricted());
		CHECK(&s.AtIndex(1) == &s.AtIndex(0));
		CHECK((s.Alpha().array() == m.array()).all());
	}

	SECTION("const access never exits Empty")
	{
		SecChem::Spinful<MatrixXd> s;
		const SecChem::Spinful<MatrixXd>& cs = s;

		(void)cs[SecChem::Spin::Alpha];
		(void)cs.Alpha();
		(void)cs.Beta();
		(void)cs.AtIndex(1);

		CHECK(s.IsEmpty());
	}
}

TEST_CASE("Spinful empty state and assignment lattice")
{
	using SpinfulMatrix = SecChem::Spinful<MatrixXd>;

	const MatrixXd a = MatrixXd::Constant(2, 2, 1.0);
	const MatrixXd b = MatrixXd::Constant(2, 2, 2.0);

	SECTION("expression construction from a restricted source adopts count 1")
	{
		SpinfulMatrix source(a);
		SpinfulMatrix target(AsExpression(source));

		CHECK(target.IsRestricted());
		CHECK((target.Alpha().array() == a.array()).all());
	}

	SECTION("expression construction from an unrestricted source adopts count 2")
	{
		SpinfulMatrix source{a, b};
		SpinfulMatrix target(AsExpression(source));

		CHECK(target.IndependentComponentCount() == 2);
		CHECK((target.Alpha().array() == a.array()).all());
		CHECK((target.Beta().array() == b.array()).all());
	}

	SECTION("first expression assignment adopts the source count from Empty")
	{
		SpinfulMatrix restrictedSource(a);
		SpinfulMatrix unrestrictedSource{a, b};

		SpinfulMatrix adoptOne;
		adoptOne = AsExpression(restrictedSource);
		CHECK(adoptOne.IsRestricted());
		CHECK_FALSE(adoptOne.IsEmpty());

		SpinfulMatrix adoptTwo;
		adoptTwo = AsExpression(unrestrictedSource);
		CHECK(adoptTwo.IndependentComponentCount() == 2);
		CHECK_FALSE(adoptTwo.IsEmpty());
	}

	SECTION("a component write exits Empty to restricted; a subsequent 1<-2 assignment then throws")
	{
		SpinfulMatrix target;
		target[SecChem::Spin::Beta] = a;
		CHECK(target.IsRestricted());

		SpinfulMatrix source{a, b};
		REQUIRE_THROWS_AS(target = AsExpression(source), SecUtility::InvalidOperationException);
	}

	SECTION("copy and copy-assignment from an Empty object stay Empty")
	{
		SpinfulMatrix empty;

		SpinfulMatrix copied(empty);
		CHECK(copied.IsEmpty());

		SpinfulMatrix assigned;
		assigned = empty;
		CHECK(assigned.IsEmpty());
	}

	SECTION("1<-1 expression assignment copies the single block")
	{
		SpinfulMatrix source(b);
		SpinfulMatrix target(a);

		target = AsExpression(source);

		CHECK(target.IsRestricted());
		CHECK((target.Alpha().array() == b.array()).all());
	}

	SECTION("2<-2 expression assignment copies both blocks")
	{
		SpinfulMatrix source{b, a};
		SpinfulMatrix target{a, b};

		target = AsExpression(source);

		CHECK((target.Alpha().array() == b.array()).all());
		CHECK((target.Beta().array() == a.array()).all());
	}

	SECTION("2<-1 expression assignment broadcasts the single block into both")
	{
		SpinfulMatrix source(a);
		SpinfulMatrix target{a, b};

		target = AsExpression(source);

		CHECK(target.IndependentComponentCount() == 2);
		CHECK((target.Alpha().array() == a.array()).all());
		CHECK((target.Beta().array() == a.array()).all());
	}

	SECTION("1<-2 on a non-empty restricted target throws naming both counts")
	{
		SpinfulMatrix target(a);
		SpinfulMatrix source{a, b};

		REQUIRE_THROWS_MATCHES(target = AsExpression(source),
		                       SecUtility::InvalidOperationException,
		                       Catch::Matchers::MessageMatches(Catch::Matchers::ContainsSubstring("count 1")));
		REQUIRE_THROWS_MATCHES(target = AsExpression(source),
		                       SecUtility::InvalidOperationException,
		                       Catch::Matchers::MessageMatches(Catch::Matchers::ContainsSubstring("count 2")));
	}

	SECTION("object copy/move assignment replaces the count wholesale")
	{
		SpinfulMatrix source{a, b};

		SpinfulMatrix copyAssigned;
		copyAssigned = source;
		CHECK(copyAssigned.IndependentComponentCount() == 2);
		CHECK_FALSE(copyAssigned.IsEmpty());

		SpinfulMatrix moveAssigned;
		moveAssigned = std::move(source);
		CHECK(moveAssigned.IndependentComponentCount() == 2);
	}

	SECTION("compound assignment on Empty adopts the expression count (scalar components)")
	{
		SecChem::Spinful<double> target;
		target += SecChem::Spinful<double>{1.0, 2.0};

		CHECK(target.IndependentComponentCount() == 2);
		CHECK(target[SecChem::Spin::Alpha] == 1.0);
		CHECK(target[SecChem::Spin::Beta] == 2.0);

		SecChem::Spinful<double> adoptedOne;
		adoptedOne += SecChem::Spinful<double>{3.0};

		CHECK(adoptedOne.IsRestricted());
		CHECK(adoptedOne[SecChem::Spin::Alpha] == 3.0);
	}

	SECTION("compound assignment of a raw value on Empty exits to restricted")
	{
		SecChem::Spinful<double> target;
		target += 2.0;

		CHECK(target.IsRestricted());
		CHECK(target[SecChem::Spin::Alpha] == 2.0);
	}

	SECTION("compound 2<-1 broadcasts the single source block")
	{
		SpinfulMatrix target{a, b};
		target += SpinfulMatrix{b};

		CHECK(target.IndependentComponentCount() == 2);
		CHECK((target.Alpha().array() == (a + b).array()).all());
		CHECK((target.Beta().array() == (b + b).array()).all());
	}

	SECTION("compound 1<-2 on a non-empty restricted target throws naming both counts")
	{
		SpinfulMatrix target(a);
		SpinfulMatrix source{a, b};

		REQUIRE_THROWS_MATCHES(target += source,
		                       SecUtility::InvalidOperationException,
		                       Catch::Matchers::MessageMatches(Catch::Matchers::ContainsSubstring("count 1")));
		REQUIRE_THROWS_MATCHES(target += source,
		                       SecUtility::InvalidOperationException,
		                       Catch::Matchers::MessageMatches(Catch::Matchers::ContainsSubstring("count 2")));
	}
}

TEST_CASE("Spinful intra ops")
{
	const MatrixXd one = MatrixXd::Constant(2, 2, 1.0);
	const MatrixXd two = MatrixXd::Constant(2, 2, 2.0);
	const MatrixXd three = MatrixXd::Constant(2, 2, 3.0);
	const MatrixXd four = MatrixXd::Constant(2, 2, 4.0);
	const MatrixXd five = MatrixXd::Constant(2, 2, 5.0);
	const MatrixXd six = MatrixXd::Constant(2, 2, 6.0);

	const SecChem::Spinful<MatrixXd> r1(one);
	const SecChem::Spinful<MatrixXd> r2(two);
	const SecChem::Spinful<MatrixXd> u1{three, four};
	const SecChem::Spinful<MatrixXd> u2{five, six};

	SECTION("addition: restricted + restricted joins to count 1, both reads broadcast")
	{
		const auto sum = r1 + r2;

		REQUIRE(sum.IndependentComponentCount() == 1);
		CheckBlocks(sum, one + two, one + two);
	}

	SECTION("addition: restricted + unrestricted joins to count 2, the restricted block broadcasts")
	{
		const auto sum = r1 + u1;

		REQUIRE(sum.IndependentComponentCount() == 2);
		CheckBlocks(sum, one + three, one + four);
	}

	SECTION("addition: unrestricted + restricted mirrors the operand order")
	{
		const auto sum = u1 + r1;

		REQUIRE(sum.IndependentComponentCount() == 2);
		CheckBlocks(sum, three + one, four + one);
	}

	SECTION("addition: unrestricted + unrestricted stays count 2")
	{
		const auto sum = u1 + u2;

		REQUIRE(sum.IndependentComponentCount() == 2);
		CheckBlocks(sum, three + five, four + six);
	}

	SECTION("subtraction: restricted - restricted joins to count 1")
	{
		const auto difference = r1 - r2;

		REQUIRE(difference.IndependentComponentCount() == 1);
		CheckBlocks(difference, one - two, one - two);
	}

	SECTION("subtraction: restricted - unrestricted joins to count 2")
	{
		const auto difference = r1 - u1;

		REQUIRE(difference.IndependentComponentCount() == 2);
		CheckBlocks(difference, one - three, one - four);
	}

	SECTION("subtraction: unrestricted - restricted mirrors the operand order")
	{
		const auto difference = u1 - r1;

		REQUIRE(difference.IndependentComponentCount() == 2);
		CheckBlocks(difference, three - one, four - one);
	}

	SECTION("subtraction: unrestricted - unrestricted stays count 2")
	{
		const auto difference = u2 - u1;

		REQUIRE(difference.IndependentComponentCount() == 2);
		CheckBlocks(difference, five - three, six - four);
	}

	// Constant 2x2 product reference: (const c) * (const d) = 2*c*d entrywise.
	SECTION("multiplication: restricted * restricted joins to count 1")
	{
		const auto product = r1 * r2;

		REQUIRE(product.IndependentComponentCount() == 1);
		CheckBlocks(product, MatrixXd::Constant(2, 2, 4.0), MatrixXd::Constant(2, 2, 4.0));
	}

	SECTION("multiplication: restricted * unrestricted joins to count 2")
	{
		const auto product = r1 * u1;

		REQUIRE(product.IndependentComponentCount() == 2);
		CheckBlocks(product, MatrixXd::Constant(2, 2, 6.0), MatrixXd::Constant(2, 2, 8.0));
	}

	SECTION("multiplication: unrestricted * restricted mirrors the operand order")
	{
		const auto product = u1 * r1;

		REQUIRE(product.IndependentComponentCount() == 2);
		CheckBlocks(product, MatrixXd::Constant(2, 2, 6.0), MatrixXd::Constant(2, 2, 8.0));
	}

	SECTION("multiplication: unrestricted * unrestricted stays count 2")
	{
		const auto product = u1 * u2;

		REQUIRE(product.IndependentComponentCount() == 2);
		CheckBlocks(product, MatrixXd::Constant(2, 2, 30.0), MatrixXd::Constant(2, 2, 48.0));
	}

	SECTION("unary minus negates per block and keeps the operand count")
	{
		const auto negated = -u1;

		REQUIRE(negated.IndependentComponentCount() == 2);
		CheckBlocks(negated, -three, -four);

		const auto negatedRestricted = -r1;

		REQUIRE(negatedRestricted.IndependentComponentCount() == 1);
		CheckBlocks(negatedRestricted, -one, -one);
	}

	SECTION("sum nodes stay lazy; Evaluate materializes an owner")
	{
		const auto sum = r1 + u1;

		static_assert(!std::is_same_v<std::decay_t<decltype(sum)>, SecChem::Spinful<MatrixXd>>);

		const auto materialized = sum.Evaluate();

		static_assert(std::is_same_v<std::decay_t<decltype(materialized)>, SecChem::Spinful<MatrixXd>>);
		REQUIRE(materialized.IndependentComponentCount() == 2);
		CheckBlocks(materialized, one + three, one + four);
	}

	SECTION("expression construction consumes node values")
	{
		const SecChem::Spinful<MatrixXd> sum = r1 + u1;

		REQUIRE(sum.IndependentComponentCount() == 2);
		CheckBlocks(sum, one + three, one + four);
	}

	SECTION("division: scalar components divide per block")
	{
		const SecChem::Spinful<double> x(3.0);
		const SecChem::Spinful<double> y(4.0);
		const SecChem::Spinful<double> u{2.0, 8.0};

		const auto quotient = x / u;

		REQUIRE(quotient.IndependentComponentCount() == 2);
		CHECK(quotient[SecChem::Spin::Alpha] == 1.5);
		CHECK(quotient[SecChem::Spin::Beta] == 0.375);

		const auto quotientRestricted = x / y;

		REQUIRE(quotientRestricted.IndependentComponentCount() == 1);
		CHECK(quotientRestricted[SecChem::Spin::Alpha] == 0.75);
		CHECK(quotientRestricted[SecChem::Spin::Beta] == 0.75);

		const auto reversed = u / x;

		REQUIRE(reversed.IndependentComponentCount() == 2);
		CHECK(reversed[SecChem::Spin::Alpha] == 2.0 / 3.0);
		CHECK(reversed[SecChem::Spin::Beta] == 8.0 / 3.0);

		static_assert(std::is_same_v<SecChem::Traits<decltype(x / u)>::Evaluated,
		                             SecChem::Spinful<double>>);
	}
}

TEST_CASE("Spinful inter ops")
{
	const MatrixXd a1 = MatrixXd::Constant(2, 2, 1.0);
	const MatrixXd b1 = MatrixXd::Constant(2, 2, 2.0);
	const MatrixXd a2 = MatrixXd::Constant(2, 2, 3.0);
	const MatrixXd b2 = MatrixXd::Constant(2, 2, 4.0);
	const MatrixXd m = MatrixXd::Constant(2, 2, 6.0);
	const MatrixXd c2 = MatrixXd::Constant(2, 2, 2.0);

	const SecChem::Spinful<MatrixXd> r(a1);
	const SecChem::Spinful<MatrixXd> u{a2, b2};

	// Matrix ± scalar is meaningless — a scalar is a 1-by-1 matrix — and Eigen defines
	// no such operators. The free inter nodes construct generically but fail at
	// component evaluation for matrix components (availability-by-call, spec §6.5);
	// additive scalar arithmetic is exercised on scalar components below, and the
	// scalar compound operator+=/-= overloads do not exist at all.
	SECTION("scalar multiplication on the right reaches every logical component")
	{
		const auto product = u * 5.0;

		REQUIRE(product.IndependentComponentCount() == 2);
		CheckBlocks(product, MatrixXd::Constant(2, 2, 15.0), MatrixXd::Constant(2, 2, 20.0));

		const auto restrictedProduct = r * 2.0;

		REQUIRE(restrictedProduct.IndependentComponentCount() == 1);
		CheckBlocks(restrictedProduct, MatrixXd::Constant(2, 2, 2.0), MatrixXd::Constant(2, 2, 2.0));
	}

	SECTION("scalar on the left (the 2 * s form)")
	{
		const auto product = 2 * u;

		REQUIRE(product.IndependentComponentCount() == 2);
		CheckBlocks(product, MatrixXd::Constant(2, 2, 6.0), MatrixXd::Constant(2, 2, 8.0));

		const auto scaled = 2.5 * r;

		REQUIRE(scaled.IndependentComponentCount() == 1);
		CheckBlocks(scaled, MatrixXd::Constant(2, 2, 2.5), MatrixXd::Constant(2, 2, 2.5));
	}

	SECTION("division by a scalar")
	{
		REQUIRE((u / 2.0).IndependentComponentCount() == 2);
		CheckBlocks(u / 2.0, MatrixXd::Constant(2, 2, 1.5), MatrixXd::Constant(2, 2, 2.0));
	}

	SECTION("plain matrix on the right applies to every logical component")
	{
		const auto sum = u + m;

		REQUIRE(sum.IndependentComponentCount() == 2);
		CheckBlocks(sum, a2 + m, b2 + m);

		const auto restrictedSum = r + m;

		REQUIRE(restrictedSum.IndependentComponentCount() == 1);
		CheckBlocks(restrictedSum, a1 + m, a1 + m);
	}

	SECTION("plain matrix on the left")
	{
		const auto sum = m + u;

		REQUIRE(sum.IndependentComponentCount() == 2);
		CheckBlocks(sum, m + a2, m + b2);

		const auto restrictedSum = m + r;

		REQUIRE(restrictedSum.IndependentComponentCount() == 1);
		CheckBlocks(restrictedSum, m + a1, m + a1);
	}

	SECTION("scalar components take the inter path in both directions")
	{
		const SecChem::Spinful<double> s(3.0);
		const SecChem::Spinful<double> t{1.0, 2.0};

		REQUIRE((2.0 * t).IndependentComponentCount() == 2);
		CHECK((2.0 * t)[SecChem::Spin::Alpha] == 2.0);
		CHECK((2.0 * t)[SecChem::Spin::Beta] == 4.0);
		CHECK((s + 1.0)[SecChem::Spin::Alpha] == 4.0);
		CHECK((1.0 + s)[SecChem::Spin::Beta] == 4.0);
		CHECK((t - 1.0)[SecChem::Spin::Alpha] == 0.0);
		CHECK((1.0 - t)[SecChem::Spin::Beta] == -1.0);
		CHECK((2.0 / s)[SecChem::Spin::Alpha] == 2.0 / 3.0);
		CHECK((s * 2)[SecChem::Spin::Alpha] == 6.0);
	}

	SECTION("compound assignment with an expression rhs obeys the lattice")
	{
		SecChem::Spinful<MatrixXd> sum{a1, b1};
		sum += u;
		CheckBlocks(sum, a1 + a2, b1 + b2);

		SecChem::Spinful<MatrixXd> difference{a1, b1};
		difference -= u;
		CheckBlocks(difference, a1 - a2, b1 - b2);

		// 2<-1 compound broadcast: both blocks compound with the single source block.
		SecChem::Spinful<MatrixXd> product{a1, b1};
		product *= SecChem::Spinful<MatrixXd>{c2};

		REQUIRE(product.IndependentComponentCount() == 2);
		CheckBlocks(product, a1 * c2, b1 * c2);
	}

	SECTION("compound assignment with a raw component")
	{
		SecChem::Spinful<MatrixXd> target{a1, b1};
		target += m;
		CheckBlocks(target, a1 + m, b1 + m);

		SecChem::Spinful<MatrixXd> difference{a1, b1};
		difference -= m;
		CheckBlocks(difference, a1 - m, b1 - m);

		SecChem::Spinful<MatrixXd> restricted(a1);
		restricted *= c2;

		REQUIRE(restricted.IsRestricted());
		CheckBlocks(restricted, a1 * c2, a1 * c2);
	}

	SECTION("compound scaling by a scalar on matrix components")
	{
		SecChem::Spinful<MatrixXd> target{a1, b1};
		target *= 2.0;
		CheckBlocks(target, MatrixXd::Constant(2, 2, 2.0), MatrixXd::Constant(2, 2, 4.0));

		target /= 4.0;
		CheckBlocks(target, MatrixXd::Constant(2, 2, 0.5), MatrixXd::Constant(2, 2, 1.0));
	}

	SECTION("scalar components compound too")
	{
		SecChem::Spinful<double> s(3.0);
		s += 1.5;
		CHECK(s[SecChem::Spin::Alpha] == 4.5);
		s -= 0.5;
		s *= 2.0;
		CHECK(s[SecChem::Spin::Alpha] == 8.0);
		s /= 4;
		CHECK(s[SecChem::Spin::Alpha] == 2.0);

		// An int argument took the scalar overload above; a double argument selects the
		// raw-component overload (the scalar one is disabled when Scalar == TSpinless).
		s /= 2.0;
		CHECK(s[SecChem::Spin::Alpha] == 1.0);

		SecChem::Spinful<double> t{1.0, 2.0};
		t += 10.0;
		CHECK(t[SecChem::Spin::Beta] == 12.0);

		SecChem::Spinful<double> d{8.0, 4.0};
		d /= SecChem::Spinful<double>{2.0, 4.0};
		CHECK(d[SecChem::Spin::Alpha] == 4.0);
		CHECK(d[SecChem::Spin::Beta] == 1.0);
	}
}

TEST_CASE("Spinful chained expressions")
{
	const MatrixXd one = MatrixXd::Constant(2, 2, 1.0);
	const MatrixXd two = MatrixXd::Constant(2, 2, 2.0);
	const MatrixXd three = MatrixXd::Constant(2, 2, 3.0);
	const MatrixXd four = MatrixXd::Constant(2, 2, 4.0);
	const MatrixXd five = MatrixXd::Constant(2, 2, 5.0);
	const MatrixXd six = MatrixXd::Constant(2, 2, 6.0);

	const SecChem::Spinful<MatrixXd> a{one, two};
	const SecChem::Spinful<MatrixXd> b{three, four};
	const SecChem::Spinful<MatrixXd> c{five, six};

	const MatrixXd expectedAlpha = (one + three) * five;
	const MatrixXd expectedBeta = (two + four) * six;

	SECTION("(a + b) * c computes per block over the node tree")
	{
		const auto product = (a + b) * c;

		REQUIRE(product.IndependentComponentCount() == 2);
		CheckBlocks(product, expectedAlpha, expectedBeta);
	}

	SECTION("named intermediates carry the node tree by value and stay safe")
	{
		const auto sum = a + b;
		const auto product = sum * c;

		static_assert(std::is_same_v<std::decay_t<decltype(sum)>,
		                             SecChem::SpinfulIntraBinaryOp<std::plus<>,
		                                                           SecChem::Spinful<MatrixXd>,
		                                                           SecChem::Spinful<MatrixXd>>>);
		REQUIRE(product.IndependentComponentCount() == 2);
		CheckBlocks(product, expectedAlpha, expectedBeta);
	}

	SECTION("depth-three chains keep computing per block")
	{
		const auto deep = (a + b) * c + a;

		REQUIRE(deep.IndependentComponentCount() == 2);
		CheckBlocks(deep, expectedAlpha + one, expectedBeta + two);
	}

	SECTION("Evaluate materializes the chain into its Evaluated owner")
	{
		const auto materialized = ((a + b) * c).Evaluate();

		static_assert(std::is_same_v<std::decay_t<decltype(materialized)>, SecChem::Spinful<MatrixXd>>);
		REQUIRE(materialized.IndependentComponentCount() == 2);
		CheckBlocks(materialized, expectedAlpha, expectedBeta);
	}

	SECTION("Evaluated and Scalar types come from the per-block Plain (spec §7.7)")
	{
		const SecChem::Spinful<VectorXd> vector{VectorXd::Constant(2, 2.0), VectorXd::Constant(2, 3.0)};
		const auto matrixTimesVector = a * vector;

		static_assert(std::is_same_v<SecChem::Traits<decltype(a * vector)>::Evaluated,
		                             SecChem::Spinful<VectorXd>>);
		REQUIRE(matrixTimesVector.IndependentComponentCount() == 2);
		CHECK((matrixTimesVector[SecChem::Spin::Alpha].array() == Eigen::VectorXd::Constant(2, 4.0).array()).all());
		CHECK((matrixTimesVector[SecChem::Spin::Beta].array() == Eigen::VectorXd::Constant(2, 12.0).array()).all());

		static_assert(std::is_same_v<SecChem::Traits<decltype((a + b) * c)>::Evaluated,
		                             SecChem::Spinful<MatrixXd>>);
		static_assert(std::is_same_v<SecChem::Traits<decltype((a + b) * c)>::Scalar, double>);
	}

	SECTION("inter nodes are directly constructible; their free operators arrive in P2.2")
	{
		const SecChem::SpinfulInterBinaryOpComponentRhs<std::plus<>, SecChem::Spinful<MatrixXd>, MatrixXd>
		        spinfulPlusComponent(a, three);

		REQUIRE(spinfulPlusComponent.IndependentComponentCount() == 2);
		CheckBlocks(spinfulPlusComponent, one + three, two + three);

		const SecChem::SpinfulInterBinaryOpComponentLhs<std::plus<>, MatrixXd, SecChem::Spinful<MatrixXd>>
		        componentPlusSpinful{three, a};

		REQUIRE(componentPlusSpinful.IndependentComponentCount() == 2);
		CheckBlocks(componentPlusSpinful, three + one, three + two);

		const SecChem::Spinful<MatrixXd> restricted(one);
		const SecChem::SpinfulInterBinaryOpComponentLhs<std::multiplies<>, MatrixXd, SecChem::Spinful<MatrixXd>>
		        restrictedScaled(two, restricted);

		REQUIRE(restrictedScaled.IndependentComponentCount() == 1);
		CheckBlocks(restrictedScaled, two * one, two * one);
	}
}

TEST_CASE("Spinful foreign family opt-in")
{
	const ForeignMatrix fa{2.0};
	const ForeignMatrix fb{3.0};
	const SecChem::Spinful<ForeignMatrix> ab{fa, fb};

	SECTION("intra ops evaluate through the family")
	{
		const SecChem::Spinful<ForeignMatrix> sum = ab + ab;
		REQUIRE(sum.IndependentComponentCount() == 2);
		CHECK(sum[SecChem::Spin::Alpha].Value() == 4.0);
		CHECK(sum[SecChem::Spin::Beta].Value() == 6.0);

		const SecChem::Spinful<ForeignMatrix> difference = ab - ab;
		CHECK(difference[SecChem::Spin::Alpha].Value() == 0.0);
		CHECK(difference[SecChem::Spin::Beta].Value() == 0.0);

		const SecChem::Spinful<ForeignMatrix> product = ab * ab;
		CHECK(product[SecChem::Spin::Alpha].Value() == 4.0);
		CHECK(product[SecChem::Spin::Beta].Value() == 9.0);
	}

	SECTION("sum-of-sums nests to depth >= 2 inside one expression")
	{
		const SecChem::Spinful<ForeignMatrix> cd{ForeignMatrix{7.0}, ForeignMatrix{11.0}};
		const SecChem::Spinful<ForeignMatrix> deep = (ab + cd) + ab;

		static_assert(std::is_same_v<SecChem::Traits<decltype((ab + cd) + ab)>::SpinComponent,
		                             ForeignSum<ForeignSum<ForeignMatrix, ForeignMatrix>, ForeignMatrix>>);
		static_assert(std::is_same_v<SecChem::Traits<decltype((ab + cd) + ab)>::Evaluated,
		                             SecChem::Spinful<ForeignMatrix>>);

		REQUIRE(deep.IndependentComponentCount() == 2);
		CHECK(deep[SecChem::Spin::Alpha].Value() == 11.0);
		CHECK(deep[SecChem::Spin::Beta].Value() == 17.0);
	}

	SECTION("inter ops accept foreign component operands on both sides")
	{
		const SecChem::Spinful<ForeignMatrix> shifted = ab + fb;
		REQUIRE(shifted.IndependentComponentCount() == 2);
		CHECK(shifted[SecChem::Spin::Alpha].Value() == 5.0);
		CHECK(shifted[SecChem::Spin::Beta].Value() == 6.0);

		const SecChem::Spinful<ForeignMatrix> reduced = ab - fb;
		CHECK(reduced[SecChem::Spin::Alpha].Value() == -1.0);
		CHECK(reduced[SecChem::Spin::Beta].Value() == 0.0);

		const SecChem::Spinful<ForeignMatrix> scaledLhs = fa * ab;
		CHECK(scaledLhs[SecChem::Spin::Alpha].Value() == 4.0);
		CHECK(scaledLhs[SecChem::Spin::Beta].Value() == 6.0);

		const SecChem::Spinful<ForeignMatrix> scaledScalar = 2.0 * ab;
		CHECK(scaledScalar[SecChem::Spin::Alpha].Value() == 4.0);
		CHECK(scaledScalar[SecChem::Spin::Beta].Value() == 6.0);

		const SecChem::Spinful<ForeignMatrix> scaledRhs = ab * 2.0;
		CHECK(scaledRhs[SecChem::Spin::Alpha].Value() == 4.0);
		CHECK(scaledRhs[SecChem::Spin::Beta].Value() == 6.0);

		const SecChem::Spinful<ForeignMatrix> restricted{fa};
		const SecChem::Spinful<ForeignMatrix> restrictedShifted = restricted + fb;
		REQUIRE(restrictedShifted.IndependentComponentCount() == 1);
		CHECK(restrictedShifted[SecChem::Spin::Alpha].Value() == 5.0);
		CHECK(restrictedShifted[SecChem::Spin::Beta].Value() == 5.0);
	}

	SECTION("a foreign expression operand nests by value")
	{
		auto nested = ab + (fa + fb);

		static_assert(std::is_same_v<SecChem::Traits<decltype(nested)>::Evaluated,
		                             SecChem::Spinful<ForeignMatrix>>);

		// The ForeignSum temporary died with its statement; the node holds its own copy
		// (inherited NestByValue = true), so late evaluation stays valid.
		const SecChem::Spinful<ForeignMatrix> evaluated = nested.Evaluate();
		CHECK(evaluated[SecChem::Spin::Alpha].Value() == 7.0);
		CHECK(evaluated[SecChem::Spin::Beta].Value() == 8.0);
	}

	SECTION("compound assignment with foreign right-hand sides")
	{
		SecChem::Spinful<ForeignMatrix> target{fa, fb};

		target += ab;
		CHECK(target[SecChem::Spin::Alpha].Value() == 4.0);
		CHECK(target[SecChem::Spin::Beta].Value() == 6.0);

		target += fa;
		CHECK(target[SecChem::Spin::Alpha].Value() == 6.0);
		CHECK(target[SecChem::Spin::Beta].Value() == 8.0);

		const SecChem::Spinful<ForeignMatrix> restricted{ForeignMatrix{5.0}};
		target += restricted;
		CHECK(target[SecChem::Spin::Alpha].Value() == 11.0);
		CHECK(target[SecChem::Spin::Beta].Value() == 13.0);

		target += ab + ab;
		CHECK(target[SecChem::Spin::Alpha].Value() == 15.0);
		CHECK(target[SecChem::Spin::Beta].Value() == 19.0);
	}

	SECTION("compound 1<-2 with foreign components throws naming both counts")
	{
		SecChem::Spinful<ForeignMatrix> target{fa};

		REQUIRE_THROWS_MATCHES(target += ab,
		                       SecUtility::InvalidOperationException,
		                       Catch::Matchers::MessageMatches(Catch::Matchers::ContainsSubstring("count 1")));
		REQUIRE_THROWS_MATCHES(target += ab,
		                       SecUtility::InvalidOperationException,
		                       Catch::Matchers::MessageMatches(Catch::Matchers::ContainsSubstring("count 2")));
	}
}

TEST_CASE("Spinful lifted ops")
{
	// Integer-valued blocks keep every expectation exact.
	MatrixXd alpha{{1, 2}, {1, 1}};
	MatrixXd beta{{2, 1}, {3, 2}};
	// alpha << 1, 2, 1, 1;
	// beta << 2, 1, 3, 2;
	const SecChem::Spinful<MatrixXd> unrestricted{alpha, beta};
	const SecChem::Spinful<MatrixXd> restricted(alpha);

	// Spinful<double>{}.Transpose() is ill-formed at the call site — double has no
	// transpose member — and member bodies instantiate lazily (spec §6.5), so there is
	// nothing to runtime-test for that rejection.
	SECTION("Transpose lifts Eigen's transpose onto every block")
	{
		const auto transposed = unrestricted.Transpose();

		REQUIRE(transposed.IndependentComponentCount() == 2);
		CheckBlocks(transposed, alpha.transpose(), beta.transpose());

		const auto restrictedTransposed = restricted.Transpose();

		REQUIRE(restrictedTransposed.IndependentComponentCount() == 1);
		CheckBlocks(restrictedTransposed, alpha.transpose(), alpha.transpose());
	}

	SECTION("Adjoint lifts Eigen's adjoint onto every block")
	{
		Eigen::MatrixXcd complexAlpha(2, 2);
		Eigen::MatrixXcd complexBeta(2, 2);
		complexAlpha << std::complex<double>(1.0, -1.0), std::complex<double>(2.0, 0.0), std::complex<double>(0.0, 3.0),
		        std::complex<double>(-4.0, 5.0);
		complexBeta << std::complex<double>(-2.0, -3.0), std::complex<double>(4.0, 0.5), std::complex<double>(1.0, 1.0),
		        std::complex<double>(0.0, -6.0);
		const SecChem::Spinful<Eigen::MatrixXcd> complexUnrestricted{complexAlpha, complexBeta};

		const auto adjoint = complexUnrestricted.Adjoint();

		REQUIRE(adjoint.IndependentComponentCount() == 2);
		CHECK((adjoint[SecChem::Spin::Alpha].array() == complexAlpha.adjoint().array()).all());
		CHECK((adjoint[SecChem::Spin::Beta].array() == complexBeta.adjoint().array()).all());

		const SecChem::Spinful<Eigen::MatrixXcd> complexRestricted(complexAlpha);

		const auto restrictedAdjoint = complexRestricted.Adjoint();

		REQUIRE(restrictedAdjoint.IndependentComponentCount() == 1);
		CHECK((restrictedAdjoint[SecChem::Spin::Alpha].array() == complexAlpha.adjoint().array()).all());
		CHECK((restrictedAdjoint[SecChem::Spin::Beta].array() == complexAlpha.adjoint().array()).all());
	}

	SECTION("Conjugate lifts Eigen's conjugate onto every block")
	{
		Eigen::MatrixXcd complexAlpha(2, 2);
		Eigen::MatrixXcd complexBeta(2, 2);
		complexAlpha << std::complex<double>(1.0, -1.0), std::complex<double>(2.0, 0.0), std::complex<double>(0.0, 3.0),
		        std::complex<double>(-4.0, 5.0);
		complexBeta << std::complex<double>(-2.0, -3.0), std::complex<double>(4.0, 0.5), std::complex<double>(1.0, 1.0),
		        std::complex<double>(0.0, -6.0);
		const SecChem::Spinful<Eigen::MatrixXcd> complexUnrestricted{complexAlpha, complexBeta};

		const auto conjugated = complexUnrestricted.Conjugate();

		REQUIRE(conjugated.IndependentComponentCount() == 2);
		CHECK((conjugated[SecChem::Spin::Alpha].array() == complexAlpha.conjugate().array()).all());
		CHECK((conjugated[SecChem::Spin::Beta].array() == complexBeta.conjugate().array()).all());

		const SecChem::Spinful<Eigen::MatrixXcd> complexRestricted(complexAlpha);

		const auto restrictedConjugated = complexRestricted.Conjugate();

		REQUIRE(restrictedConjugated.IndependentComponentCount() == 1);
		CHECK((restrictedConjugated[SecChem::Spin::Beta].array() == complexAlpha.conjugate().array()).all());
	}

	SECTION("Inverse lifts Eigen's inverse onto every block")
	{
		const auto inverted = unrestricted.Inverse();

		REQUIRE(inverted.IndependentComponentCount() == 2);
		CheckBlocks(inverted, alpha.inverse(), beta.inverse());

		const auto restrictedInverted = restricted.Inverse();

		REQUIRE(restrictedInverted.IndependentComponentCount() == 1);
		CheckBlocks(restrictedInverted, alpha.inverse(), alpha.inverse());
	}

	SECTION("AsDiagonal lifts Eigen's asDiagonal onto every vector block")
	{
		VectorXd energyAlpha(2);
		VectorXd energyBeta(2);
		energyAlpha << 10, 100;
		energyBeta << 1000, 10000;
		const SecChem::Spinful<VectorXd> energies{energyAlpha, energyBeta};

		const auto diagonalized = energies.AsDiagonal();

		static_assert(std::is_same_v<SecChem::Traits<decltype(energies.AsDiagonal())>::Evaluated,
		                             SecChem::Spinful<Eigen::DiagonalMatrix<double, Eigen::Dynamic>>>);
		REQUIRE(diagonalized.IndependentComponentCount() == 2);
		CHECK((diagonalized[SecChem::Spin::Alpha].diagonal().array() == energyAlpha.array()).all());
		CHECK((diagonalized[SecChem::Spin::Beta].diagonal().array() == energyBeta.array()).all());

		const SecChem::Spinful<VectorXd> restrictedEnergies(energyAlpha);

		const auto restrictedDiagonalized = restrictedEnergies.AsDiagonal();

		REQUIRE(restrictedDiagonalized.IndependentComponentCount() == 1);
		CHECK((restrictedDiagonalized[SecChem::Spin::Beta].diagonal().array() == energyAlpha.array()).all());
	}

	SECTION("Diagonal lifts Eigen's diagonal onto every block")
	{
		const auto diagonals = unrestricted.Diagonal();

		REQUIRE(diagonals.IndependentComponentCount() == 2);
		CHECK((diagonals[SecChem::Spin::Alpha].array() == alpha.diagonal().array()).all());
		CHECK((diagonals[SecChem::Spin::Beta].array() == beta.diagonal().array()).all());

		const auto restrictedDiagonals = restricted.Diagonal();

		REQUIRE(restrictedDiagonals.IndependentComponentCount() == 1);
		CHECK((restrictedDiagonals[SecChem::Spin::Beta].array() == alpha.diagonal().array()).all());
	}

	SECTION("Cast converts every block to the target scalar")
	{
		const Eigen::MatrixXf floatAlpha = Eigen::MatrixXf::Constant(2, 3, 0.5f);
		const Eigen::MatrixXf floatBeta = Eigen::MatrixXf::Constant(2, 3, -1.25f);
		const SecChem::Spinful<Eigen::MatrixXf> floatPair{floatAlpha, floatBeta};

		const auto widened = floatPair.Cast<double>();
		static_assert(std::is_same_v<SecChem::Traits<decltype(floatPair.Cast<double>())>::Evaluated,
		                             SecChem::Spinful<MatrixXd>>);
		REQUIRE(widened.IndependentComponentCount() == 2);
		CHECK((widened[SecChem::Spin::Alpha].array() == floatAlpha.cast<double>().array()).all());
		CHECK((widened[SecChem::Spin::Beta].array() == floatBeta.cast<double>().array()).all());

		// Cast composes through further nodes and stays lazy until Evaluate.
		const auto doubled = (floatPair.Cast<double>() * 2.0).Evaluate();
		static_assert(std::is_same_v<std::decay_t<decltype(doubled)>, SecChem::Spinful<MatrixXd>>);
		CHECK((doubled[SecChem::Spin::Alpha].array() == (2.0 * floatAlpha.cast<double>()).array()).all());

		const SecChem::Spinful<Eigen::MatrixXf> floatRestricted(floatAlpha);
		REQUIRE(floatRestricted.Cast<double>().IndependentComponentCount() == 1);
		CHECK((floatRestricted.Cast<double>()[SecChem::Spin::Beta].array() == floatAlpha.cast<double>().array()).all());
	}

	SECTION("Cast converts scalar components through static_cast")
	{
		const SecChem::Spinful<double> energies{0.5, -1.5};

		const auto narrowed = energies.Cast<float>();
		static_assert(std::is_same_v<SecChem::Traits<decltype(energies.Cast<float>())>::Evaluated,
		                             SecChem::Spinful<float>>);
		CHECK(narrowed[SecChem::Spin::Alpha] == 0.5f);
		CHECK(narrowed[SecChem::Spin::Beta] == -1.5f);

		const SecChem::Spinful<std::complex<float>> complexPair{std::complex<float>(0.5f, 1.0f),
		                                                        std::complex<float>(-2.0f, 0.25f)};
		const auto complexWidened = complexPair.Cast<std::complex<double>>();
		static_assert(std::is_same_v<
		              SecChem::Traits<decltype(complexPair.Cast<std::complex<double>>())>::Evaluated,
		              SecChem::Spinful<std::complex<double>>>);
		CHECK(complexWidened[SecChem::Spin::Alpha] == std::complex<double>(0.5, 1.0));
		CHECK(complexWidened[SecChem::Spin::Beta] == std::complex<double>(-2.0, 0.25));

		// A same-scalar cast is the identity.
		CHECK(energies.Cast<double>().Evaluate().EqualsTo(energies));
	}

	SECTION("Re and Im extract the complex parts of every block")
	{
		Eigen::MatrixXcd complexAlpha(2, 2);
		Eigen::MatrixXcd complexBeta(2, 2);
		complexAlpha << std::complex<double>(1.0, -1.0), std::complex<double>(2.0, 0.0), std::complex<double>(0.0, 3.0),
		        std::complex<double>(-4.0, 5.0);
		complexBeta << std::complex<double>(-2.0, -3.0), std::complex<double>(4.0, 0.5), std::complex<double>(1.0, 1.0),
		        std::complex<double>(0.0, -6.0);
		const SecChem::Spinful<Eigen::MatrixXcd> complexPair{complexAlpha, complexBeta};

		const auto realParts = complexPair.Re();
		static_assert(std::is_same_v<SecChem::Traits<decltype(complexPair.Re())>::Evaluated,
		                             SecChem::Spinful<MatrixXd>>);
		REQUIRE(realParts.IndependentComponentCount() == 2);
		CHECK((realParts[SecChem::Spin::Alpha].array() == complexAlpha.real().array()).all());
		CHECK((realParts[SecChem::Spin::Beta].array() == complexBeta.real().array()).all());

		const auto imaginaryParts = complexPair.Im();
		CHECK((imaginaryParts[SecChem::Spin::Alpha].array() == complexAlpha.imag().array()).all());
		CHECK((imaginaryParts[SecChem::Spin::Beta].array() == complexBeta.imag().array()).all());

		const SecChem::Spinful<Eigen::MatrixXcd> complexRestricted(complexAlpha);
		REQUIRE(complexRestricted.Re().IndependentComponentCount() == 1);
		CHECK((complexRestricted.Re()[SecChem::Spin::Beta].array() == complexAlpha.real().array()).all());

		// Node operand: statement-scoped component expressions, as everywhere.
		CHECK(((complexPair + complexPair).Im()[SecChem::Spin::Alpha].array() == (2 * complexAlpha.imag()).array())
		              .all());
	}

	SECTION("Re is the identity and Im zero on real components")
	{
		const auto realParts = unrestricted.Re().Evaluate();
		static_assert(std::is_same_v<std::decay_t<decltype(realParts)>, SecChem::Spinful<MatrixXd>>);
		CHECK(realParts.EqualsTo(unrestricted));

		const auto imaginaryParts = unrestricted.Im().Evaluate();
		CHECK((imaginaryParts[SecChem::Spin::Alpha].array() == Eigen::MatrixXd::Zero(2, 2).array()).all());
		CHECK((imaginaryParts[SecChem::Spin::Beta].array() == Eigen::MatrixXd::Zero(2, 2).array()).all());

		const SecChem::Spinful<std::complex<double>> scalarComplex{std::complex<double>(0.5, 2.5),
		                                                           std::complex<double>(-1.0, -0.25)};
		const auto re = scalarComplex.Re().Evaluate();
		static_assert(std::is_same_v<std::decay_t<decltype(re)>, SecChem::Spinful<double>>);
		CHECK(re.Alpha() == 0.5);
		CHECK(re.Beta() == -1.0);
		CHECK(scalarComplex.Im().Evaluate().Alpha() == 2.5);

		// Arithmetic scalar components take the std::real / std::imag fallbacks.
		const SecChem::Spinful<double> realScalars{0.5, -1.5};
		const auto scalarRe = realScalars.Re().Evaluate();
		static_assert(std::is_same_v<std::decay_t<decltype(scalarRe)>, SecChem::Spinful<double>>);
		CHECK(scalarRe.Alpha() == 0.5);
		CHECK(scalarRe.Beta() == -1.5);
		CHECK(realScalars.Im().Evaluate().Alpha() == 0.0);
	}

	SECTION("Componentwise lifts an arbitrary callable onto every block")
	{
		const double factor = 3.0;

		// The callable is decay-copied and invoked at AtIndex time — the capture is
		// stored, not applied eagerly.
		const auto scaled = unrestricted.Componentwise([factor](const auto& block) { return factor * block; });

		REQUIRE(scaled.IndependentComponentCount() == 2);
		CheckBlocks(scaled, 3.0 * alpha, 3.0 * beta);

		const auto restrictedScaled = restricted.Componentwise([factor](const auto& block) { return factor * block; });

		REQUIRE(restrictedScaled.IndependentComponentCount() == 1);
		CheckBlocks(restrictedScaled, 3.0 * alpha, 3.0 * alpha);

		// A callable may morph the component type: matrix blocks to scalars.
		const auto totals = unrestricted.Componentwise([](const auto& block) { return block.sum(); }).Evaluate();

		static_assert(std::is_same_v<std::decay_t<decltype(totals)>, SecChem::Spinful<double>>);
		REQUIRE(totals.IndependentComponentCount() == 2);
		CHECK(totals[SecChem::Spin::Alpha] == 5.0);  // 1 + 2 + 1 + 1
		CHECK(totals[SecChem::Spin::Beta] == 8.0);   // 2 + 1 + 3 + 2
	}

	SECTION("Lifted ops compose through expression nodes")
	{
		const auto sumThenTranspose = (unrestricted + unrestricted).Transpose();

		REQUIRE(sumThenTranspose.IndependentComponentCount() == 2);
		CheckBlocks(sumThenTranspose, (alpha + alpha).transpose(), (beta + beta).transpose());

		const auto doubleTransposed = unrestricted.Transpose().Transpose();

		REQUIRE(doubleTransposed.IndependentComponentCount() == 2);
		CheckBlocks(doubleTransposed, alpha, beta);
	}

	SECTION("acceptance-bar composition evaluates; Trace arrives in P3.2")
	{
		VectorXd energyAlpha(2);
		VectorXd energyBeta(2);
		energyAlpha << 10, 100;
		energyBeta << 1000, 10000;
		const SecChem::Spinful<VectorXd> energies{energyAlpha, energyBeta};

		const auto composed = (unrestricted.Transpose() * energies.AsDiagonal()).Evaluate();

		// The transposed operand carries Eigen's RowMajorBit, so the materialized
		// owner is a row-major Matrix — plain_matrix_type propagates the flag. Pin the
		// scalar, not the storage order.
		static_assert(std::is_same_v<SecChem::Traits<decltype(unrestricted.Transpose()
		                                                                         * energies.AsDiagonal())>::Scalar,
		                             double>);
		REQUIRE(composed.IndependentComponentCount() == 2);

		// Manual per-block trace stands in for Trace() until P3.2 (spec §6.5 bar).
		const double alphaTrace = composed[SecChem::Spin::Alpha](0, 0) + composed[SecChem::Spin::Alpha](1, 1);
		const double betaTrace = composed[SecChem::Spin::Beta](0, 0) + composed[SecChem::Spin::Beta](1, 1);

		CHECK(alphaTrace == 110.0);   // 1*10 + 1*100
		CHECK(betaTrace == 22000.0);  // 2*1000 + 2*10000
	}
}

TEST_CASE("Spinful reductions")
{
	MatrixXd alpha(2, 2);
	MatrixXd beta(2, 2);
	alpha << 1, 2, 1, 1;
	beta << 2, 1, 3, 2;
	const SecChem::Spinful<MatrixXd> density{alpha, beta};
	const SecChem::Spinful<MatrixXd> restricted(alpha);

	SECTION("Unrestricted folds combine both components")
	{
		CHECK(density.Trace() == 6.0);                        // 2 + 4
		CHECK(density.Determinant() == Catch::Approx(-1.0));  // det(A) * det(B) = -1 * 1
		CHECK(density.Sum() == 13.0);                         // 5 + 8
		CHECK(density.SquaredNorm() == 25.0);                 // 7 + 18
		CHECK(density.Norm() == 5.0);                         // sqrt(25) - exact
		CHECK(density.MaxCoeff() == 3.0);
		CHECK(density.MinCoeff() == 1.0);
		CHECK(density.IsFinite());
	}

	SECTION("Restricted folds double the single component")
	{
		CHECK(restricted.Trace() == 4.0);                       // 2 * 2
		CHECK(restricted.Determinant() == Catch::Approx(1.0));  // det(A)^2 = (-1)^2
		CHECK(restricted.Sum() == 10.0);                        // 2 * 5
		CHECK(restricted.SquaredNorm() == 14.0);                // 2 * 7
		CHECK(restricted.Norm() == Approx(std::sqrt(14.0)));    // sqrt(2) * ||alpha||
		CHECK(restricted.MaxCoeff() == 2.0);                    // undoubled
		CHECK(restricted.MinCoeff() == 1.0);                    // undoubled
		CHECK(restricted.IsFinite());
	}

	SECTION("Determinant is the block-diagonal product; SpinResolved members stay per-slot")
	{
		const auto spinResolved = density.SpinResolvedDeterminant();
		CHECK(spinResolved.IndependentComponentCount() == 2);
		CHECK(spinResolved.Alpha() == Catch::Approx(-1.0));  // det(A)
		CHECK(spinResolved.Beta() == Catch::Approx(1.0));    // det(B)

		const auto spinResolvedTrace = density.SpinResolvedTrace();
		CHECK(spinResolvedTrace.IndependentComponentCount() == 2);
		CHECK(spinResolvedTrace.Alpha() == 2.0);  // tr(A)
		CHECK(spinResolvedTrace.Beta() == 4.0);   // tr(B)

		const auto restrictedResolved = restricted.SpinResolvedDeterminant();
		CHECK(restrictedResolved.IndependentComponentCount() == 1);
		CHECK(restrictedResolved.IsRestricted());
		CHECK(restrictedResolved.Alpha() == Catch::Approx(-1.0));  // the shared block's determinant

		const auto restrictedResolvedTrace = restricted.SpinResolvedTrace();
		CHECK(restrictedResolvedTrace.IndependentComponentCount() == 1);
		CHECK(restrictedResolvedTrace.Alpha() == 2.0);  // the shared block's trace

		// nodes take the members: det(2A) = 4*det(A), tr(2A) = 2*tr(A)
		CHECK((density + density).Determinant() == Catch::Approx(-16.0));  // -4 * 4
		const auto nodeResolved = (density + density).SpinResolvedDeterminant();
		CHECK(nodeResolved.IndependentComponentCount() == 2);
		CHECK(nodeResolved.Alpha() == Catch::Approx(-4.0));
		CHECK(nodeResolved.Beta() == Catch::Approx(4.0));
		const auto nodeResolvedTrace = (density + density).SpinResolvedTrace();
		CHECK(nodeResolvedTrace.Alpha() == 4.0);
		CHECK(nodeResolvedTrace.Beta() == 8.0);

		// Empty reads as restricted: the 0-by-0 determinant is the empty product 1,
		// the 0-by-0 trace is 0
		const SecChem::Spinful<MatrixXd> empty;
		CHECK(empty.Determinant() == Catch::Approx(1.0));
		const auto emptyResolved = empty.SpinResolvedDeterminant();
		CHECK(emptyResolved.IndependentComponentCount() == 1);
		CHECK(emptyResolved.Alpha() == Catch::Approx(1.0));
		const auto emptyResolvedTrace = empty.SpinResolvedTrace();
		CHECK(emptyResolvedTrace.IndependentComponentCount() == 1);
		CHECK(emptyResolvedTrace.Alpha() == 0.0);
	}

	SECTION("LpNorm folds the per-component norms for each order")
	{
		// |alpha|_1 = 5, |alpha|_inf = 2; |beta|_1 = 8, |beta|_inf = 3.
		CHECK(density.LpNorm<1>() == 13.0);               // 5 + 8
		CHECK(restricted.LpNorm<1>() == 10.0);            // 2 * 5
		CHECK(density.LpNorm<Eigen::Infinity>() == 3.0);  // undoubled, like MaxCoeff
		CHECK(restricted.LpNorm<Eigen::Infinity>() == 2.0);

		CHECK(density.LpNorm<2>() == Approx(5.0));  // agrees with Norm()
		CHECK(density.LpNorm<2>() == Approx(density.Norm()));
		CHECK(restricted.LpNorm<2>() == Approx(std::sqrt(14.0)));         // 2^(1/2) * ||alpha||_2
		CHECK(density.LpNorm<3>() == Approx(std::pow(55.0, 1.0 / 3.0)));  // 11 + 44
	}

	SECTION("LpNorm takes the scalar fallback, nodes, and the zero of Empty")
	{
		const SecChem::Spinful<double> scalarPair{0.5, -1.5};
		const SecChem::Spinful<double> scalarRestricted(0.5);

		CHECK(scalarPair.LpNorm<1>() == 2.0);
		CHECK(scalarPair.LpNorm<2>() == Approx(std::sqrt(2.5)));
		CHECK(scalarPair.LpNorm<Eigen::Infinity>() == 1.5);
		CHECK(scalarRestricted.LpNorm<1>() == 1.0);  // 2 * 0.5
		CHECK(scalarRestricted.LpNorm<Eigen::Infinity>() == 0.5);

		CHECK((density + density).LpNorm<1>() == 26.0);  // expression node operand

		const SecChem::Spinful<MatrixXd> empty;
		CHECK(empty.LpNorm<1>() == 0.0);
		CHECK(empty.LpNorm<Eigen::Infinity>() == 0.0);
	}

	SECTION("Reductions run on expression components")
	{
		CHECK((density + density).Trace() == 12.0);
		CHECK((density + density).Sum() == 26.0);
		CHECK((density + density).SquaredNorm() == 100.0);  // quadratic: 4*7 + 4*18
		CHECK((density + density).MaxCoeff() == 6.0);
		CHECK((density + density).MinCoeff() == 2.0);
		CHECK((density + density).IsFinite());
		CHECK(density.Transpose().Trace() == 6.0);        // tr(A^T) = tr(A)
		CHECK((restricted + restricted).Trace() == 8.0);  // restricted node, still doubled
	}

	SECTION("IsFinite detects non-finite components")
	{
		MatrixXd withInfinity(2, 2);
		MatrixXd withNaN(2, 2);
		withInfinity << 1.0, 2.0, std::numeric_limits<double>::infinity(), 4.0;
		withNaN << 1.0, std::numeric_limits<double>::quiet_NaN(), 3.0, 4.0;

		CHECK(SecChem::Spinful<MatrixXd>{alpha, beta}.IsFinite());
		CHECK_FALSE(SecChem::Spinful<MatrixXd>{withInfinity, beta}.IsFinite());
		CHECK_FALSE(SecChem::Spinful<MatrixXd>{alpha, withNaN}.IsFinite());
		CHECK_FALSE(SecChem::Spinful<MatrixXd>(withNaN).IsFinite());
	}

	SECTION("Scalar components reduce through the scalar fallbacks")
	{
		const SecChem::Spinful<double> scalarPair{0.5, -1.5};
		const SecChem::Spinful<double> scalarRestricted(0.5);

		CHECK(scalarPair.Trace() == -1.0);         // 0.5 + (-1.5)
		CHECK(scalarPair.Determinant() == -0.75);  // 0.5 * (-1.5) — the 1-by-1 identity
		CHECK(scalarPair.Sum() == -1.0);
		CHECK(scalarPair.SquaredNorm() == 2.5);  // 0.25 + 2.25
		CHECK(scalarPair.Norm() == Approx(std::sqrt(2.5)));
		CHECK(scalarPair.MaxCoeff() == 0.5);
		CHECK(scalarPair.MinCoeff() == -1.5);
		CHECK(scalarPair.IsFinite());

		CHECK(scalarRestricted.Trace() == 1.0);         // 2 * 0.5
		CHECK(scalarRestricted.Determinant() == 0.25);  // 0.5 * 0.5
		CHECK(scalarRestricted.Sum() == 1.0);
		CHECK(scalarRestricted.SquaredNorm() == 0.5);              // 2 * 0.25
		CHECK(scalarRestricted.Norm() == Approx(std::sqrt(0.5)));  // |a| without std::abs
		CHECK(scalarRestricted.MaxCoeff() == 0.5);
		CHECK(scalarRestricted.MinCoeff() == 0.5);
		CHECK(scalarRestricted.IsFinite());

		const auto scalarResolved = scalarPair.SpinResolvedDeterminant();
		CHECK(scalarResolved.IndependentComponentCount() == 2);
		CHECK(scalarResolved.Alpha() == 0.5);
		CHECK(scalarResolved.Beta() == -1.5);

		const auto scalarResolvedTrace = scalarPair.SpinResolvedTrace();
		CHECK(scalarResolvedTrace.Alpha() == 0.5);
		CHECK(scalarResolvedTrace.Beta() == -1.5);
	}

	SECTION("Complex scalar components")
	{
		const std::complex<double> alphaComplex(3.0, 4.0);
		const std::complex<double> betaComplex(0.0, -2.0);
		const SecChem::Spinful<std::complex<double>> complexPair{alphaComplex, betaComplex};

		CHECK(complexPair.SquaredNorm() == Approx(29.0));  // 25 + 4
		CHECK(complexPair.Norm() == Approx(std::sqrt(29.0)));
		CHECK(complexPair.Trace().real() == 3.0);  // identity fallback, then +
		CHECK(complexPair.Trace().imag() == 2.0);
		CHECK(complexPair.Sum().real() == 3.0);
		CHECK(complexPair.Sum().imag() == 2.0);
		const std::complex<double> determinantProduct(8.0, -6.0);  // (3+4i) * (-2i)
		CHECK(complexPair.Determinant() == determinantProduct);
		const auto complexTraceResolved = complexPair.SpinResolvedTrace();
		CHECK(complexTraceResolved.Alpha() == alphaComplex);
		CHECK(complexTraceResolved.Beta() == betaComplex);
		CHECK(complexPair.IsFinite());

		const SecChem::Spinful<std::complex<double>> complexNonFinite{
		        std::complex<double>(1.0, std::numeric_limits<double>::infinity()), betaComplex};
		CHECK_FALSE(complexNonFinite.IsFinite());
	}

	SECTION("Acceptance bar: band energy as a reduction")
	{
		VectorXd energyAlpha(2);
		VectorXd energyBeta(2);
		energyAlpha << 10, 100;
		energyBeta << 1000, 10000;
		const SecChem::Spinful<VectorXd> energies{energyAlpha, energyBeta};

		const auto bandEnergy = (density.Transpose() * energies.AsDiagonal()).Trace();
		CHECK(bandEnergy == 22110.0);  // 110 + 22000
	}
}

TEST_CASE("Spinful conversions")
{
	MatrixXd alpha(2, 2);
	MatrixXd beta(2, 2);
	alpha << 1, 2, 1, 1;
	beta << 2, 1, 3, 2;
	const SecChem::Spinful<MatrixXd> density{alpha, beta};
	const SecChem::Spinful<MatrixXd> restricted(alpha);

	SECTION("ToUnrestricted expands a restricted spinful")
	{
		const SecChem::Spinful<MatrixXd> expanded = restricted.ToUnrestricted();
		CHECK_FALSE(expanded.IsRestricted());
		CheckBlocks(expanded, alpha, alpha);
		CHECK(&expanded.Alpha() != &expanded.Beta());  // two independent copies
	}

	SECTION("ToUnrestricted is the identity on an unrestricted spinful")
	{
		const SecChem::Spinful<MatrixXd> again = density.ToUnrestricted();
		CHECK_FALSE(again.IsRestricted());
		CheckBlocks(again, alpha, beta);
	}

	SECTION("ToRestricted succeeds on agreeing components")
	{
		const SecChem::Spinful<MatrixXd> agreeing{alpha, alpha};
		const SecChem::Spinful<MatrixXd> result = agreeing.ToRestricted();
		CHECK(result.IsRestricted());
		CheckBlocks(result, alpha, alpha);
	}

	SECTION("ToRestricted probes the Frobenius tolerance boundary")
	{
		MatrixXd shifted(2, 2);
		shifted << 1, 2, 1 + 1e-6, 1;  // Frobenius distance 1e-6
		const SecChem::Spinful<MatrixXd> nearPair{alpha, shifted};

		const SecChem::Spinful<MatrixXd> accepted = nearPair.ToRestricted(2e-6);
		CHECK(accepted.IsRestricted());
		CheckBlocks(accepted, alpha, alpha);  // Alpha survives the restriction

		REQUIRE_THROWS_AS(nearPair.ToRestricted(0.5e-6), SecUtility::InvalidOperationException);
		REQUIRE_THROWS_AS(nearPair.ToRestricted(), SecUtility::InvalidOperationException);  // default tolerance
	}

	SECTION("ToRestricted is trivially legal on a restricted spinful")
	{
		const SecChem::Spinful<MatrixXd> result = restricted.ToRestricted();
		CHECK(result.IsRestricted());
		CheckBlocks(result, alpha, alpha);
	}

	SECTION("ToRestrictedAverage averages both components")
	{
		const SecChem::Spinful<MatrixXd> averaged = density.ToRestrictedAverage();
		CHECK(averaged.IsRestricted());
		CheckBlocks(averaged, (alpha + beta) * 0.5, (alpha + beta) * 0.5);
	}

	SECTION("ToRestrictedAverage reproduces a restricted spinful exactly")
	{
		const SecChem::Spinful<MatrixXd> averaged = restricted.ToRestrictedAverage();
		CHECK(averaged.IsRestricted());
		CheckBlocks(averaged, alpha, alpha);
	}

	SECTION("IsNearlyRestricted is the predicate of ToRestricted")
	{
		CHECK(restricted.IsNearlyRestricted());     // broadcast distance 0
		CHECK(density.IsNearlyRestricted(5.0));     // Frobenius distance sqrt(7) < 5
		CHECK_FALSE(density.IsNearlyRestricted());  // default tolerance is dummy precision

		MatrixXd shifted(2, 2);
		shifted << 1, 2, 1 + 1e-6, 1;  // Frobenius distance 1e-6
		const SecChem::Spinful<MatrixXd> nearPair{alpha, shifted};
		CHECK(nearPair.IsNearlyRestricted(2e-6));
		CHECK_FALSE(nearPair.IsNearlyRestricted(0.5e-6));
		REQUIRE_THROWS_AS(nearPair.ToRestricted(0.5e-6), SecUtility::InvalidOperationException);

		// Empty reports true: its broadcast distance is zero and ToRestricted succeeds —
		// the documented asymmetry against IsRestricted().
		const SecChem::Spinful<MatrixXd> empty;
		CHECK(empty.IsNearlyRestricted());
		CHECK_FALSE(empty.IsRestricted());
	}

	SECTION("ToSpinFlipped swaps the blocks of an unrestricted spinful")
	{
		const SecChem::Spinful<MatrixXd> flipped = density.ToSpinFlipped();
		CHECK(flipped.IndependentComponentCount() == 2);
		CheckBlocks(flipped, beta, alpha);
		CHECK(flipped.ToSpinFlipped().EqualsTo(density));  // an involution

		CHECK(restricted.ToSpinFlipped().EqualsTo(restricted));  // spin-symmetric: its own flip
		CHECK(restricted.ToSpinFlipped().IsRestricted());

		const SecChem::Spinful<MatrixXd> empty;
		CHECK(empty.ToSpinFlipped().IsEmpty());

		const SecChem::Spinful<double> scalarPair{0.5, -1.5};
		CHECK(scalarPair.ToSpinFlipped().Alpha() == -1.5);
		CHECK(scalarPair.ToSpinFlipped().Beta() == 0.5);
	}

	SECTION("FlipSpin swaps in place and returns the same object")
	{
		SecChem::Spinful<MatrixXd> flippable{alpha, beta};

		SecChem::Spinful<MatrixXd>& self = flippable.FlipSpin();
		CHECK(&self == &flippable);  // in-place, chainable
		CheckBlocks(flippable, beta, alpha);

		flippable.FlipSpin();
		CheckBlocks(flippable, alpha, beta);  // an involution

		SecChem::Spinful<MatrixXd> restrictedCopy(alpha);
		restrictedCopy.FlipSpin();  // no-op: its own flip
		CHECK(restrictedCopy.EqualsTo(restricted));
		CHECK(restrictedCopy.IsRestricted());

		SecChem::Spinful<MatrixXd> emptyFlip;
		CHECK(emptyFlip.FlipSpin().IsEmpty());

		SecChem::Spinful<double> scalarPair{0.5, -1.5};
		scalarPair.FlipSpin();
		CHECK(scalarPair.Alpha() == -1.5);
		CHECK(scalarPair.Beta() == 0.5);
	}

	SECTION("SpinSum adds the logical components, evaluated")
	{
		const MatrixXd unrestrictedSum = density.SpinSum();
		static_assert(std::is_same_v<decltype(unrestrictedSum), const MatrixXd>);  // evaluated, not an expression
		CHECK((unrestrictedSum.array() == (alpha + beta).array()).all());

		const MatrixXd restrictedSum = restricted.SpinSum();
		CHECK((restrictedSum.array() == (2 * alpha).array()).all());
	}

	SECTION("SpinSum runs on expressions")
	{
		const MatrixXd sum = (density + density).SpinSum();
		static_assert(std::is_same_v<decltype(sum), const MatrixXd>);
		CHECK((sum.array() == (2 * (alpha + beta)).array()).all());
	}

	SECTION("Scalar components restrict through the distance fallback")
	{
		const SecChem::Spinful<double> scalarPair{0.5, 0.5 + 1e-9};

		const SecChem::Spinful<double> accepted = scalarPair.ToRestricted(2e-9);
		CHECK(accepted.IsRestricted());
		CHECK(accepted.Alpha() == 0.5);

		REQUIRE_THROWS_AS(scalarPair.ToRestricted(0.5e-9), SecUtility::InvalidOperationException);

		const SecChem::Spinful<double> agreeing{0.5, 0.5};
		CHECK(agreeing.ToRestricted().Alpha() == 0.5);  // default tolerance accepts an exact pair
	}
}

TEST_CASE("Spinful iteration")
{
	MatrixXd alpha(2, 2);
	MatrixXd beta(2, 2);
	alpha << 1, 2, 1, 1;
	beta << 2, 1, 3, 2;
	const SecChem::Spinful<MatrixXd> density{alpha, beta};
	const SecChem::Spinful<MatrixXd> restricted(alpha);

	static_assert(SecChem::AllSpins[0] == SecChem::Spin::Alpha);  // iteration order
	static_assert(SecChem::AllSpins[1] == SecChem::Spin::Beta);

	SECTION("Logical components alias on a restricted spinful")
	{
		const auto logical = restricted.LogicalComponents();
		const auto* entries = logical.begin();
		CHECK(entries[0].Spin == SecChem::AllSpins[0]);
		CHECK(entries[1].Spin == SecChem::AllSpins[1]);
		CHECK(&entries[0].Component == &entries[1].Component);  // the alias, observable
		CHECK(&entries[0].Component == &restricted.Alpha());

		std::size_t spinsSeen = 0;
		for (const auto& entry : logical)
		{
			CHECK(entry.Spin == SecChem::AllSpins[spinsSeen]);
			CHECK(entry.Component(0, 0) == 1.0);
			++spinsSeen;
		}
		CHECK(spinsSeen == 2);
	}

	SECTION("Logical components are distinct on an unrestricted spinful")
	{
		const auto logical = density.LogicalComponents();
		const auto* entries = logical.begin();
		CHECK(entries[0].Spin == SecChem::Spin::Alpha);
		CHECK(entries[1].Spin == SecChem::Spin::Beta);
		CHECK(&entries[0].Component == &density.Alpha());
		CHECK(&entries[1].Component == &density.Beta());
		CHECK(&entries[0].Component != &entries[1].Component);
		CHECK((entries[0].Component.array() == alpha.array()).all());
		CHECK((entries[1].Component.array() == beta.array()).all());
	}

	SECTION("Independent components follow the stored count")
	{
		const auto unrestricted = density.IndependentComponents();
		CHECK(std::distance(unrestricted.begin(), unrestricted.end()) == 2);
		CHECK(unrestricted.begin()->Spin == SecChem::Spin::Alpha);

		const auto restrictedRange = restricted.IndependentComponents();
		CHECK(std::distance(restrictedRange.begin(), restrictedRange.end()) == 1);

		const SecChem::Spinful<MatrixXd> empty;
		const auto emptyRange = empty.IndependentComponents();
		CHECK(std::distance(emptyRange.begin(), emptyRange.end()) == 1);  // public count 1
		CHECK(empty.IsEmpty());                                           // const iteration never exits Empty
	}
}

TEST_CASE("Spinful comparison")
{
	MatrixXd alpha(2, 2);
	MatrixXd beta(2, 2);
	alpha << 1, 2, 1, 1;
	beta << 2, 1, 3, 2;
	const SecChem::Spinful<MatrixXd> restricted(alpha);
	const SecChem::Spinful<MatrixXd> expanded{alpha, alpha};
	const SecChem::Spinful<MatrixXd> density{alpha, beta};

	SECTION("A restricted spinful equals its unrestricted expansion")
	{
		CHECK(restricted.EqualsTo(expanded));  // (A) vs (A, A): both logical pairs are (A, A)
		CHECK(expanded.EqualsTo(restricted));
	}

	SECTION("Restricted spinfuls compare on the single distance")
	{
		const SecChem::Spinful<MatrixXd> same(alpha);
		CHECK(restricted.EqualsTo(same));
		const SecChem::Spinful<MatrixXd> other(beta);
		CHECK_FALSE(restricted.EqualsTo(other));
	}

	SECTION("Unequal blocks are rejected")
	{
		CHECK_FALSE(restricted.EqualsTo(density));
		CHECK_FALSE(density.EqualsTo(restricted));
	}

	SECTION("Tolerance is respected below and above")
	{
		MatrixXd shifted(2, 2);
		shifted << 1, 2, 1 + 1e-6, 1;  // Frobenius distance 1e-6 from alpha
		const SecChem::Spinful<MatrixXd> nearPair{shifted, shifted};
		CHECK(nearPair.EqualsTo(restricted, 2e-6));
		CHECK_FALSE(nearPair.EqualsTo(restricted, 0.5e-6));
		CHECK_FALSE(nearPair.EqualsTo(restricted));  // default is dummy precision, 1e-12
	}

	SECTION("The inherited operator== is exact; NotEqualsTo mirrors EqualsTo")
	{
		MatrixXd perturbed(2, 2);
		perturbed << 1, 2, 1, 1 + 1e-15;  // distance 1e-15, below any sane tolerance
		const SecChem::Spinful<MatrixXd> tinyPair{perturbed, perturbed};
		CHECK(restricted == expanded);  // bit-identical blocks
		CHECK(restricted != density);
		CHECK_FALSE(restricted == tinyPair);  // tolerance 0: a tiny perturbation differs
		CHECK(restricted.NotEqualsTo(tinyPair, 0.5e-15));
		CHECK_FALSE(restricted.NotEqualsTo(tinyPair, 2e-15));
	}

	SECTION("Empty spinfuls compare equal")
	{
		const SecChem::Spinful<MatrixXd> lhs;
		const SecChem::Spinful<MatrixXd> rhs;
		CHECK(lhs.EqualsTo(rhs));  // (0x0 - 0x0).norm() == 0 — the A2 empty case
		CHECK(lhs == rhs);
		// Empty vs non-Empty is a dimension mismatch — the documented eigen_assert
		// abort case; it aborts rather than throwing, so it is not runtime-tested.
	}

	SECTION("Scalar components compare through std::abs")
	{
		const SecChem::Spinful<double> a(0.5);
		const SecChem::Spinful<double> b{0.5, 0.5 + 2e-9};
		CHECK(a.EqualsTo(b, 3e-9));  // differing counts, scalar distance 2e-9
		CHECK_FALSE(a.EqualsTo(b, 1e-9));
	}

	SECTION("The free operators accept expressions on either side")
	{
		CHECK((expanded + expanded - expanded) == expanded);                         // node == owner
		CHECK(expanded == (expanded + expanded - expanded));                         // owner == node
		CHECK((expanded + expanded - expanded) == (restricted + restricted) * 0.5);  // node == node
		CHECK((restricted + density) != expanded);                                   // node != owner
		CHECK(expanded != (restricted + density));                                   // owner != node
	}

	SECTION("Any expression can be the comparison callee")
	{
		CHECK((expanded + expanded - expanded).EqualsTo(restricted));  // node callee, owner rhs
		CHECK((expanded + expanded - expanded).EqualsTo(expanded));

		const auto sum = restricted + density;               // (2a, a+b)
		CHECK(sum.NotEqualsTo(restricted + restricted));     // same-shell nodes through the interface
		CHECK(sum.NotEqualsTo(density + density));           // differs in Beta
		CHECK_FALSE(sum.NotEqualsTo(density + restricted));  // (2a, a+b) again
	}

	SECTION("Equality accepts mixed component types")
	{
		// Eigen mixes dynamic- and fixed-size matrices of the same scalar; mixing
		// numeric types (float vs double) is NOT allowed — Eigen demands an explicit
		// cast there, and our availability-by-call rejects it identically.
		const Eigen::Matrix2d alpha2 = alpha;
		const SecChem::Spinful<Eigen::Matrix2d> restrictedFixed(alpha2);
		CHECK(restricted == restrictedFixed);
		CHECK(restrictedFixed != density);

		CHECK(restricted.EqualsTo(restrictedFixed));                 // cross-type tolerance comparison
		CHECK(restricted.EqualsTo(expanded + expanded - expanded));  // lazy rhs, no materialization

		MatrixXd shifted(2, 2);
		shifted << 1, 2, 1 + 1e-6, 1;  // Frobenius distance 1e-6 from alpha
		const SecChem::Spinful<MatrixXd> nearPair{shifted, shifted};
		CHECK(nearPair.EqualsTo(restrictedFixed, 2e-6));  // the callee's tolerance governs
		CHECK_FALSE(nearPair.EqualsTo(restrictedFixed, 0.5e-6));

		// Scalar mixing the std library allows: complex minus real distances through std::abs.
		const SecChem::Spinful<std::complex<double>> complexPair{0.5, 0.25};
		const SecChem::Spinful<double> scalarPair(0.5);
		CHECK(scalarPair.EqualsTo(complexPair, 0.3));
		CHECK_FALSE(scalarPair.EqualsTo(complexPair, 0.2));
	}

	// Integral components: the exact routes (the free and inherited operator==/!=,
	// NotEqualsTo at tolerance 0) accept integral scalars — covered in the
	// scalar-components test case; only the tolerance-parameterized routes (EqualsTo
	// with a tolerance, Restricted) remain compile-time errors for them. A
	// static_assert failure cannot be runtime-tested, hence no SECTION here.
}

TEST_CASE("Spinful streaming")
{
	MatrixXd alpha(2, 2);
	alpha << 1, 2, 3, 4;
	MatrixXd beta(2, 2);
	beta << 5, 6, 7, 8;

	SECTION("A restricted spinful streams the shared block once")
	{
		const SecChem::Spinful<MatrixXd> restricted(alpha);
		std::ostringstream stream;
		stream << restricted;
		CHECK(stream.str() == "(Alpha|Beta):\n1 2\n3 4\n");
	}

	SECTION("An unrestricted spinful streams both blocks")
	{
		const SecChem::Spinful<MatrixXd> pair{alpha, beta};
		std::ostringstream stream;
		stream << pair;
		CHECK(stream.str() == "(Alpha):\n1 2\n3 4\n(Beta):\n5 6\n7 8\n");
	}

	SECTION("An Empty spinful streams its state, not a restricted frame")
	{
		// v2: Empty prints its own marker (the honest tri-state) instead of the v1
		// restricted frame around 0x0 blocks.
		const SecChem::Spinful<MatrixXd> empty;
		std::ostringstream stream;
		stream << empty;
		CHECK(stream.str() == "(Empty):\n");
	}

	SECTION("Lazy expressions stream unevaluated through the CRTP base")
	{
		const SecChem::Spinful<MatrixXd> restricted(alpha);
		const SecChem::Spinful<MatrixXd> pair{alpha, beta};

		std::ostringstream sumStream;
		sumStream << (pair + pair);
		CHECK(sumStream.str() == "(Alpha):\n2 4\n6 8\n(Beta):\n10 12\n14 16\n");

		std::ostringstream restrictedStream;
		restrictedStream << (restricted + restricted);
		CHECK(restrictedStream.str() == "(Alpha|Beta):\n2 4\n6 8\n");
	}

	SECTION("Scalar components stream through the built-in operator<<")
	{
		const SecChem::Spinful<double> energies{0.5, 1.25};
		std::ostringstream stream;
		stream << energies;
		CHECK(stream.str() == "(Alpha):\n0.5\n(Beta):\n1.25\n");
	}
}

TEST_CASE("Spinful scalar components")
{
	using SpinfulScalar = SecChem::Spinful<double>;
	// Int32 lives at global scope (Raw/Int.hpp default), not in SecUtility.
	using IntSpinfulInt = SecChem::Spinful<Int32>;

	SECTION("Arithmetic runs through the scalar component operators")
	{
		const SpinfulScalar restricted(3.0);
		const SpinfulScalar pair{1.0, 2.0};

		const SpinfulScalar sum = restricted + pair;  // broadcast: (3+1, 3+2)
		CHECK(sum.Alpha() == 4.0);
		CHECK(sum.Beta() == 5.0);

		const SpinfulScalar scaled = 2.0 * pair;  // inter: the component is on the left
		CHECK(scaled.Alpha() == 2.0);
		CHECK(scaled.Beta() == 4.0);

		const SpinfulScalar negated = -restricted;
		CHECK(negated.Alpha() == -3.0);
		CHECK(negated.Beta() == -3.0);  // broadcast

		const SpinfulScalar quotient = pair / SpinfulScalar(2.0);
		CHECK(quotient.Alpha() == 0.5);
		CHECK(quotient.Beta() == 1.0);
	}

	SECTION("Reductions take the scalar adapter fallbacks")
	{
		const SpinfulScalar restricted(3.0);
		const SpinfulScalar pair{1.0, 2.0};

		CHECK(restricted.Trace() == 6.0);  // logical components: 2*3
		CHECK(pair.Trace() == 3.0);        // identity adapter: 1 + 2
		CHECK(pair.Sum() == 3.0);
		CHECK(pair.SquaredNorm() == 5.0);  // 1*1 + 2*2
		CHECK(pair.Norm() == std::sqrt(5.0));
		CHECK(pair.MaxCoeff() == 2.0);
		CHECK(pair.MinCoeff() == 1.0);
		CHECK(pair.IsFinite());
		CHECK_FALSE(SpinfulScalar(std::numeric_limits<double>::infinity()).IsFinite());
	}

	SECTION("Compound assignment runs on the raw-component overloads")
	{
		const SpinfulScalar pair{1.0, 2.0};

		SpinfulScalar scaled = pair;
		scaled *= 2.0;  // TSpinless == the scalar: the component overload takes it
		CHECK(scaled.Alpha() == 2.0);
		CHECK(scaled.Beta() == 4.0);

		SpinfulScalar quotient = pair;
		quotient /= 2.0;
		CHECK(quotient.Alpha() == 0.5);
		CHECK(quotient.Beta() == 1.0);

		SpinfulScalar shifted = pair;
		shifted += 1.0;  // the scalar compound overload is disabled for TSpinless == double
		CHECK(shifted.Alpha() == 2.0);
		CHECK(shifted.Beta() == 3.0);
		shifted -= 4.0;
		CHECK(shifted.Alpha() == -2.0);
		CHECK(shifted.Beta() == -1.0);

		SpinfulScalar combined = pair;
		combined += pair + pair;  // expression compound, count 2
		CHECK(combined.Alpha() == 3.0);
		CHECK(combined.Beta() == 6.0);
	}

	SECTION("Integral scalar components stay integral")
	{
		const IntSpinfulInt pair{Int32{2}, Int32{3}};

		const IntSpinfulInt sum = pair + IntSpinfulInt(Int32{5});  // broadcast: (7, 8)
		CHECK(sum.Alpha() == 7);
		CHECK(sum.Beta() == 8);

		IntSpinfulInt product = pair;
		product *= Int32{4};
		CHECK(product.Alpha() == 8);
		CHECK(product.Beta() == 12);

		CHECK(pair.Sum() == 5);
		CHECK(pair.SquaredNorm() == 13);  // 2*2 + 3*3 in integer arithmetic
		CHECK(pair.MaxCoeff() == 3);
		CHECK(pair.MinCoeff() == 2);
	}

	SECTION("Integral scalar components compare exactly")
	{
		const IntSpinfulInt samePair{Int32{2}, Int32{3}};
		const IntSpinfulInt pair{Int32{2}, Int32{3}};
		const IntSpinfulInt other{Int32{2}, Int32{4}};
		const IntSpinfulInt single(Int32{2});
		const IntSpinfulInt expansion{Int32{2}, Int32{2}};

		CHECK(pair == samePair);
		CHECK_FALSE(pair == other);
		CHECK(pair != other);
		CHECK_FALSE(pair != samePair);
		CHECK(single == expansion);       // restricted equals its unrestricted expansion
		CHECK(single.NotEqualsTo(pair));  // the interface route at tolerance 0 (dummy precision of int is 0)
		CHECK_FALSE(single.NotEqualsTo(expansion));
	}

	SECTION("Scalar factories branch by count")
	{
		const SpinfulScalar zero = SpinfulScalar::Zero();
		CHECK(zero.IsRestricted());
		CHECK(zero.Alpha() == 0.0);

		const SpinfulScalar zeroPair = SpinfulScalar::Zero(2);
		CHECK(zeroPair.IndependentComponentCount() == 2);
		CHECK(zeroPair.Beta() == 0.0);

		const SpinfulScalar one = SpinfulScalar::One(2);
		CHECK(one.Alpha() == 1.0);
		CHECK(one.Beta() == 1.0);

		const IntSpinfulInt intOne = IntSpinfulInt::One();
		CHECK(intOne.Alpha() == 1);
	}

	SECTION("Conversions run through std::abs")
	{
		const SpinfulScalar pair{1.0, 3.0};

		const SpinfulScalar expanded = SpinfulScalar(2.0).ToUnrestricted();
		CHECK(expanded.IndependentComponentCount() == 2);
		CHECK(expanded.Alpha() == 2.0);
		CHECK(expanded.Beta() == 2.0);
		CHECK(pair.ToUnrestricted().Beta() == 3.0);  // identity on the unrestricted

		const SpinfulScalar averaged = pair.ToRestrictedAverage();
		CHECK(averaged.IsRestricted());
		CHECK(averaged.Alpha() == 2.0);

		CHECK(pair.ToRestricted(2.0).Alpha() == 1.0);  // |1-3| = 2 at the boundary
		REQUIRE_THROWS_AS(pair.ToRestricted(1.5), SecUtility::InvalidOperationException);
		const SpinfulScalar nearAgreeing{0.5, 0.5 + 1e-13};  // inside the 1e-12 default
		CHECK(nearAgreeing.ToRestricted().Alpha() == 0.5);

		CHECK(pair.SpinSum() == 4.0);
		CHECK(SpinfulScalar(2.0).SpinSum() == 4.0);  // restricted doubles

		// Integral components: ToUnrestricted and SpinSum are available; ToRestricted,
		// ToRestrictedAverage and IsNearlyRestricted static_assert on the floating-point
		// requirement — a compile-time rejection, not runtime-testable here.
		const IntSpinfulInt intPair{Int32{2}, Int32{5}};
		CHECK(intPair.ToUnrestricted().Beta() == 5);
		CHECK(intPair.SpinSum() == 7);
	}

	SECTION("Empty adoption semantics on scalars")
	{
		SpinfulScalar empty;
		CHECK(std::as_const(empty).IsEmpty());       // non-const component access exits Empty
		CHECK(std::as_const(empty).Alpha() == 0.0);  // value-initialized storage, no indeterminacy

		empty = SpinfulScalar{1.0, 2.0} + SpinfulScalar(5.0);  // expression: Empty adopts count 2
		CHECK(empty.IndependentComponentCount() == 2);
		CHECK(empty.Alpha() == 6.0);
		CHECK(empty.Beta() == 7.0);

		SpinfulScalar adopter;
		adopter += SpinfulScalar(1.0);  // compound adoption exits Empty restricted
		CHECK(adopter.IsRestricted());
		CHECK(adopter.Alpha() == 1.0);

		SpinfulScalar exited;
		exited.Alpha() = 7.0;  // non-const component access exits Empty
		CHECK(exited.IsRestricted());
		CHECK(exited.Beta() == 7.0);  // broadcast

		IntSpinfulInt emptyInt;
		CHECK(emptyInt.IsEmpty());
		CHECK(emptyInt.Alpha() == 0);
	}
}

// Non-owning spinfuls are component compositions: Spinful<Map<T>> views foreign
// memory per block — separate buffers, one contiguous [A|B] buffer, strided — where
// Eigen's Map defines the view semantics (alias on copy, write-through assignment).
// There is no top-level Map<Spinful<T>>: Spinful defines no scalar-buffer layout to
// map over, and the restricted/unrestricted count is object state, not bytes.
TEST_CASE("Spinful mapped components")
{
	using MappedMatrix = Eigen::Map<MatrixXd>;
	using ConstMappedMatrix = Eigen::Map<const MatrixXd>;
	using MappedSpinful = SecChem::Spinful<MappedMatrix>;
	using ConstMappedSpinful = SecChem::Spinful<ConstMappedMatrix>;

	MatrixXd alpha(2, 2);
	alpha << 1, 2, 3, 4;
	MatrixXd beta(2, 2);
	beta << 5, 6, 7, 8;

	// A foreign [A|B] buffer — the contiguous layout a QC driver or file mapping hands over.
	std::vector<double> buffer(alpha.size() + beta.size());
	Eigen::Map<MatrixXd>(buffer.data(), 2, 2) = alpha;
	Eigen::Map<MatrixXd>(buffer.data() + alpha.size(), 2, 2) = beta;
	double* const pAlpha = buffer.data();
	double* const pBeta = buffer.data() + alpha.size();

	SECTION("An unrestricted view aliases the mapped buffers")
	{
		MappedSpinful density{MappedMatrix{pAlpha, 2, 2}, MappedMatrix{pBeta, 2, 2}};

		CHECK(density.IndependentComponentCount() == 2);
		CHECK(density.Alpha().data() == pAlpha);
		CHECK(density.Beta().data() == pBeta);
		CheckBlocks(density, alpha, beta);

		density.Alpha()(0, 0) = 9.0;  // component access writes through the view
		CHECK(pAlpha[0] == 9.0);

		const MappedSpinful copy = density;  // copies the mapping, not the data
		CHECK(copy.Alpha().data() == pAlpha);

		const MappedSpinful flipped = density.ToSpinFlipped();  // a permutation of views
		CHECK(flipped.Alpha().data() == pBeta);
		CHECK(flipped.Beta().data() == pAlpha);
	}

	SECTION("A restricted view broadcasts one mapped block")
	{
		MappedSpinful restricted{MappedMatrix{pAlpha, 2, 2}};

		CHECK(restricted.IsRestricted());
		CHECK(restricted.Alpha().data() == pAlpha);
		CHECK(restricted.Beta().data() == pAlpha);  // broadcast: both spins view the same memory
		CheckBlocks(restricted, alpha, alpha);
	}

	SECTION("Views are lazy-expression operands")
	{
		const MappedSpinful density{MappedMatrix{pAlpha, 2, 2}, MappedMatrix{pBeta, 2, 2}};
		const SecChem::Spinful<MatrixXd> other{beta, alpha};

		CheckBlocks(density + other, alpha + beta, beta + alpha);

		// A node over mapped components materializes to an owning spinful: the
		// component Plain of Map<T> is T.
		const auto materialized = (density + other).Evaluate();
		static_assert(std::is_same_v<decltype(materialized), const SecChem::Spinful<MatrixXd>>);
		CheckBlocks(materialized, alpha + beta, beta + alpha);
	}

	SECTION("Expression assignment writes through the view")
	{
		MappedSpinful target{MappedMatrix{pAlpha, 2, 2}, MappedMatrix{pBeta, 2, 2}};
		const SecChem::Spinful<MatrixXd> source{beta, alpha};

		target = source;
		CHECK((Eigen::Map<const MatrixXd>(pAlpha, 2, 2).array() == beta.array()).all());
		CHECK((Eigen::Map<const MatrixXd>(pBeta, 2, 2).array() == alpha.array()).all());

		target += source;  // compound assignment writes through as well
		CHECK((Eigen::Map<const MatrixXd>(pAlpha, 2, 2).array() == (2 * beta).array()).all());
		CHECK((Eigen::Map<const MatrixXd>(pBeta, 2, 2).array() == (2 * alpha).array()).all());
	}

	SECTION("The assignment lattice holds for views")
	{
		MappedSpinful restricted{MappedMatrix{pAlpha, 2, 2}};
		const SecChem::Spinful<MatrixXd> unrestricted{alpha, beta};

		REQUIRE_THROWS_AS(restricted = unrestricted, SecUtility::InvalidOperationException);
	}

	SECTION("Const-mapped views are read-only operands")
	{
		const ConstMappedSpinful density{ConstMappedMatrix{pAlpha, 2, 2}, ConstMappedMatrix{pBeta, 2, 2}};

		CHECK(density.EqualsTo(SecChem::Spinful<MatrixXd>{alpha, beta}));

		CheckBlocks(density + density, alpha + alpha, beta + beta);

		// Materializing a view into an owner detaches: the owner copies the data.
		const SecChem::Spinful<MatrixXd> owned = density;
		CheckBlocks(owned, alpha, beta);
		CHECK(owned.Alpha().data() != pAlpha);
	}
}

TEST_CASE("Spinful reference views")
{
	using RefMatrix = std::reference_wrapper<MatrixXd>;
	using ConstRefMatrix = std::reference_wrapper<const MatrixXd>;
	using RefSpinful = SecChem::Spinful<RefMatrix>;
	using ConstRefSpinful = SecChem::Spinful<ConstRefMatrix>;

	// The wrapper delegates to the referenced component's traits (D30, spec §5.3).
	static_assert(std::is_same_v<SecChem::Detail::ComponentTraits<RefMatrix>::Scalar, double>);
	static_assert(std::is_same_v<SecChem::Detail::ComponentTraits<RefMatrix>::Plain, MatrixXd>);
	static_assert(std::is_same_v<SecChem::Detail::ComponentTraits<ConstRefMatrix>::Plain, MatrixXd>);
	static_assert(SecChem::Detail::is_component<RefMatrix>::value);
	static_assert(std::is_same_v<SecChem::Detail::ComponentOf<RefMatrix>, MatrixXd>);
	static_assert(std::is_same_v<SecChem::Detail::ComponentOf<ConstRefMatrix>, const MatrixXd>);
	static_assert(std::is_same_v<RefSpinful::SpinComponent, MatrixXd>);
	static_assert(std::is_same_v<SecChem::Traits<RefSpinful>::Evaluated, RefSpinful>);

	MatrixXd alpha(2, 2);
	alpha << 1, 2, 3, 4;
	MatrixXd beta(2, 2);
	beta << 5, 6, 7, 8;

	SECTION("std::ref and std::cref deduce unrestricted views")
	{
		SecChem::Spinful density{std::ref(alpha), std::ref(beta)};
		static_assert(std::is_same_v<decltype(density), RefSpinful>);

		CHECK(density.IndependentComponentCount() == 2);
		CHECK(!density.IsRestricted());
		CHECK(!density.IsEmpty());
		CHECK(density.Alpha().data() == alpha.data());
		CHECK(density.Beta().data() == beta.data());
		CheckBlocks(density, alpha, beta);

		SecChem::Spinful constView{std::cref(alpha), std::cref(beta)};
		static_assert(std::is_same_v<decltype(constView), ConstRefSpinful>);
		CheckBlocks(constView, alpha, beta);
	}

	SECTION("Component access writes through the view")
	{
		SecChem::Spinful density{std::ref(alpha), std::ref(beta)};

		density.Alpha()(0, 0) = 9.0;  // component access writes through the view
		CHECK(alpha(0, 0) == 9.0);

		const RefSpinful copy = density;  // copies the handles, not the data
		CHECK(copy.Alpha().data() == alpha.data());

		const RefSpinful flipped = density.ToSpinFlipped();  // a permutation of views
		CHECK(flipped.Alpha().data() == beta.data());
		CHECK(flipped.Beta().data() == alpha.data());

		RefSpinful mutableFlipped = density.ToSpinFlipped();
		mutableFlipped.FlipSpin();  // swapping handles swaps the referenced blocks
		CHECK(mutableFlipped.Alpha().data() == alpha.data());

		const RefSpinful evaluated = density.Evaluate();  // a view evaluates to itself
		CHECK(evaluated.Alpha().data() == alpha.data());
	}

	SECTION("Reads and reductions run on the referenced blocks")
	{
		SecChem::Spinful density{std::ref(alpha), std::ref(beta)};

		CHECK(density.Trace() == 18.0);  // tr(A) + tr(B) = 5 + 13
		CHECK(density.Determinant() == Approx(4.0));  // det(A) * det(B) = (-2) * (-2)

		const SecChem::Spinful<double> resolved = density.SpinResolvedTrace();
		CHECK(resolved.Alpha() == 5.0);
		CHECK(resolved.Beta() == 13.0);

		int logicalCount = 0;
		double traceSum = 0.0;
		for (const auto& entry : density.LogicalComponents())
		{
			++logicalCount;
			traceSum += entry.Component.trace();
		}
		CHECK(logicalCount == 2);
		CHECK(traceSum == 18.0);
	}

	SECTION("Views are lazy-expression operands")
	{
		const SecChem::Spinful density{std::ref(alpha), std::ref(beta)};
		const SecChem::Spinful<MatrixXd> other{beta, alpha};

		CheckBlocks(density + other, alpha + beta, beta + alpha);

		// A node over view components materializes to an owning spinful: the
		// component Plain of reference_wrapper<T> is T.
		const auto materialized = (density + other).Evaluate();
		static_assert(std::is_same_v<decltype(materialized), const SecChem::Spinful<MatrixXd>>);
		CheckBlocks(materialized, alpha + beta, beta + alpha);

		// Materializing a view into an owner detaches: the owner copies the data.
		const SecChem::Spinful<MatrixXd> owned = density;
		CheckBlocks(owned, alpha, beta);
		CHECK(owned.Alpha().data() != alpha.data());
	}

	SECTION("Expression assignment and fills write through the view")
	{
		// The expectations must be detached: alpha and beta are the target's own
		// storage, so every write-through overwrites the pristine values too.
		const MatrixXd expectedAlpha = alpha;
		const MatrixXd expectedBeta = beta;
		SecChem::Spinful target{std::ref(alpha), std::ref(beta)};
		const SecChem::Spinful<MatrixXd> source{expectedBeta, expectedAlpha};

		target = source;
		CHECK((alpha.array() == expectedBeta.array()).all());
		CHECK((beta.array() == expectedAlpha.array()).all());

		target += source;  // compound assignment writes through as well
		CHECK((alpha.array() == (2 * expectedBeta).array()).all());
		CHECK((beta.array() == (2 * expectedAlpha).array()).all());

		target.SetConstant(3.0);  // in-place fills write through too
		CHECK((alpha.array() == 3.0).all());
		CHECK((beta.array() == 3.0).all());
	}

	SECTION("Const views are read-only operands")
	{
		const SecChem::Spinful constView{std::cref(alpha), std::cref(beta)};

		CHECK(constView.EqualsTo(SecChem::Spinful<MatrixXd>{alpha, beta}));

		CheckBlocks(constView + constView, alpha + alpha, beta + beta);

		const SecChem::Spinful<MatrixXd> owned = constView;
		CheckBlocks(owned, alpha, beta);
		CHECK(owned.Alpha().data() != alpha.data());
	}

	SECTION("Scalar views")
	{
		double alphaValue = 1.0;
		double betaValue = 2.0;
		SecChem::Spinful scalarView{std::ref(alphaValue), std::ref(betaValue)};
		static_assert(std::is_same_v<decltype(scalarView), SecChem::Spinful<std::reference_wrapper<double>>>);

		CHECK(scalarView.Trace() == 3.0);  // the 1-by-1 identity: tr(d) = d

		scalarView.Alpha() = 10.0;  // writes through to the referenced scalar
		CHECK(alphaValue == 10.0);

		scalarView.SetConstant(7.0);
		CHECK(alphaValue == 7.0);
		CHECK(betaValue == 7.0);

		scalarView *= 2.0;  // scalar scaling writes through
		CHECK(alphaValue == 14.0);
		CHECK(betaValue == 14.0);
	}
}

TEST_CASE("Spinful sparse components")
{
	using SparseMatrix = Eigen::SparseMatrix<double>;
	using SparseSpinful = SecChem::Spinful<SparseMatrix>;

	SparseMatrix alpha(2, 2);
	alpha.insert(0, 0) = 1.0;
	alpha.insert(1, 1) = 2.0;
	SparseMatrix beta(2, 2);
	beta.insert(0, 0) = 3.0;
	beta.insert(1, 1) = 5.0;

	SECTION("Sparse blocks compose through the storage layer")
	{
		const SparseSpinful density{alpha, beta};

		CHECK(density.IndependentComponentCount() == 2);
		CHECK(density.Alpha().isApprox(alpha));
		CHECK(density.Beta().isApprox(beta));

		const SparseSpinful flipped = density.ToSpinFlipped();
		CHECK(flipped.Alpha().isApprox(beta));
		CHECK(flipped.Beta().isApprox(alpha));
	}

	SECTION("A restricted sparse spinful broadcasts one block")
	{
		const SparseSpinful restricted{alpha};

		CHECK(restricted.IsRestricted());
		CHECK(restricted.Alpha().isApprox(alpha));
		CHECK(restricted.Beta().isApprox(alpha));  // broadcast: both spins read one block

		const SparseSpinful expanded = restricted.ToUnrestricted();
		CHECK(expanded.IndependentComponentCount() == 2);
		CHECK(expanded.Alpha().isApprox(alpha));
		CHECK(expanded.Beta().isApprox(alpha));
	}

	SECTION("Sparse reductions fold over the logical components")
	{
		const SparseSpinful density{alpha, beta};
		const SparseSpinful restricted{alpha};

		CHECK(density.Sum() == Approx(11.0));          // (1+2) + (3+5)
		CHECK(restricted.Sum() == Approx(6.0));        // 2 * (1+2)
		CHECK(density.SquaredNorm() == Approx(39.0));  // (1+4) + (9+25)
		CHECK(density.Norm() == Approx(std::sqrt(39.0)));
	}

	SECTION("Sparse arithmetic and assignment evaluate through Eigen")
	{
		const SparseSpinful density{alpha, beta};

		const SparseSpinful doubled = (density + density).Evaluate();
		static_assert(std::is_same_v<decltype(doubled), const SparseSpinful>);
		CHECK(doubled.Alpha().isApprox(2 * alpha));
		CHECK(doubled.Beta().isApprox(2 * beta));

		const auto scaled = (2.0 * density).Evaluate();
		CHECK(scaled.Alpha().isApprox(2 * alpha));
		CHECK(scaled.Beta().isApprox(2 * beta));

		SparseSpinful target{alpha};
		target += SparseSpinful{alpha};  // 1←1: stays restricted, one stored block
		CHECK(target.IsRestricted());
		CHECK(target.Alpha().isApprox(2 * alpha));
	}

	SECTION("Sparse comparisons and conversions use the distance rule")
	{
		const SparseSpinful density{alpha, beta};

		CHECK(density == density);  // exact equality: distance 0
		CHECK(density.EqualsTo(SparseSpinful{alpha, beta}));
		CHECK(!density.EqualsTo(SparseSpinful{beta, alpha}));

		const SparseSpinful pair{alpha, alpha};  // storage-unrestricted, A == B
		CHECK(pair.IsNearlyRestricted());
		CHECK(pair.ToRestricted().IsRestricted());

		const SparseSpinful average = density.ToRestrictedAverage();
		const SparseMatrix expected = (alpha + beta) * 0.5;
		CHECK(average.Alpha().isApprox(expected));
	}
}

// A mapped sparse spinful is a source operand, not an output buffer: Eigen's sparse
// assignment rebuilds the storage pattern (resize/reserve/insertBackByOuterInner),
// which a fixed foreign buffer cannot host — expression and compound assignment
// through the view are unavailable (availability-by-call, D24's sparse exception), and
// coeffRef writes only existing nonzeros.
TEST_CASE("Spinful mapped sparse components")
{
	using SparseMatrix = Eigen::SparseMatrix<double>;
	using MappedSparse = Eigen::Map<SparseMatrix>;
	using ConstMappedSparse = Eigen::Map<const SparseMatrix>;
	using MappedSparseSpinful = SecChem::Spinful<MappedSparse>;

	// Foreign [A|B] CSC buffers — the compressed layout a sparse driver hands over:
	// A = diag(1, 2), B = diag(3, 5).
	std::vector<int> outerA{0, 1, 2}, innerA{0, 1};
	std::vector<double> valuesA{1.0, 2.0};
	std::vector<int> outerB{0, 1, 2}, innerB{0, 1};
	std::vector<double> valuesB{3.0, 5.0};

	SparseMatrix alpha(2, 2);
	alpha.insert(0, 0) = 1.0;
	alpha.insert(1, 1) = 2.0;
	SparseMatrix beta(2, 2);
	beta.insert(0, 0) = 3.0;
	beta.insert(1, 1) = 5.0;

	SECTION("An unrestricted view aliases the mapped buffers")
	{
		MappedSparseSpinful density{MappedSparse{2, 2, 2, outerA.data(), innerA.data(), valuesA.data()},
		                            MappedSparse{2, 2, 2, outerB.data(), innerB.data(), valuesB.data()}};

		CHECK(density.IndependentComponentCount() == 2);
		CHECK(density.Alpha().outerIndexPtr() == outerA.data());
		CHECK(density.Alpha().valuePtr() == valuesA.data());
		CHECK(density.Beta().outerIndexPtr() == outerB.data());
		CHECK(density.Beta().valuePtr() == valuesB.data());
		CHECK(density.Alpha().coeff(1, 1) == 2.0);

		// Overwriting an existing nonzero is the entire writable surface of a sparse
		// view — and it writes through to the raw buffer.
		density.Alpha().coeffRef(0, 0) = 9.0;
		CHECK(valuesA[0] == 9.0);

		const MappedSparseSpinful copy = density;  // copies the mapping, not the data
		CHECK(copy.Alpha().valuePtr() == valuesA.data());

		const MappedSparseSpinful flipped = density.ToSpinFlipped();  // a permutation of views
		CHECK(flipped.Alpha().valuePtr() == valuesB.data());
		CHECK(flipped.Beta().valuePtr() == valuesA.data());
	}

	SECTION("A restricted view broadcasts one mapped block")
	{
		MappedSparseSpinful restricted{MappedSparse{2, 2, 2, outerA.data(), innerA.data(), valuesA.data()}};

		CHECK(restricted.IsRestricted());
		CHECK(restricted.Alpha().outerIndexPtr() == outerA.data());
		CHECK(restricted.Beta().outerIndexPtr() == outerA.data());  // broadcast: both spins view one buffer

		const MappedSparseSpinful expanded = restricted.ToUnrestricted();
		CHECK(expanded.IndependentComponentCount() == 2);
		CHECK(expanded.Alpha().outerIndexPtr() == outerA.data());
		CHECK(expanded.Beta().outerIndexPtr() == outerA.data());
	}

	SECTION("Views are lazy-expression operands")
	{
		const MappedSparseSpinful density{MappedSparse{2, 2, 2, outerA.data(), innerA.data(), valuesA.data()},
		                                  MappedSparse{2, 2, 2, outerB.data(), innerB.data(), valuesB.data()}};
		const SecChem::Spinful<SparseMatrix> owned{alpha, beta};

		// A node over mixed view/owner components materializes to the owning type.
		const auto materialized = (density + owned).Evaluate();
		static_assert(std::is_same_v<decltype(materialized), const SecChem::Spinful<SparseMatrix>>);
		CHECK(materialized.Alpha().isApprox(2 * alpha));
		CHECK(materialized.Beta().isApprox(2 * beta));

		const auto scaled = (2.0 * density).Evaluate();
		CHECK(scaled.Alpha().isApprox(2 * alpha));
		CHECK(scaled.Beta().isApprox(2 * beta));

		// Materializing a view into an owner detaches: the owner copies the data.
		const SecChem::Spinful<SparseMatrix> detached = density;
		CHECK(detached.Alpha().isApprox(alpha));
		CHECK(detached.Alpha().valuePtr() != valuesA.data());
	}

	SECTION("Const-mapped views are read-only operands (2)")
	{
		const SecChem::Spinful<ConstMappedSparse> density{
		        ConstMappedSparse{2, 2, 2, outerA.data(), innerA.data(), valuesA.data()},
		        ConstMappedSparse{2, 2, 2, outerB.data(), innerB.data(), valuesB.data()}};

		CHECK(density.Alpha().coeff(0, 0) == 1.0);
		CHECK(density.Alpha().coeff(1, 0) == 0.0);
		CHECK(density.EqualsTo(SecChem::Spinful<SparseMatrix>{alpha, beta}));
	}

	SECTION("Const-mapped views are read-only operands (1)")
	{
		const SecChem::Spinful<ConstMappedSparse> density{
		        ConstMappedSparse{2, 2, 2, outerA.data(), innerA.data(), valuesA.data()}};

		CHECK(density.Alpha().coeff(0, 0) == 1.0);
		CHECK(density.Alpha().coeff(1, 0) == 0.0);
		CHECK(density.EqualsTo(SecChem::Spinful<SparseMatrix>{alpha, alpha}));
	}

	SECTION("Empty construction holds null sparse views")
	{
		const MappedSparseSpinful empty;

		CHECK(empty.IsEmpty());
		CHECK(!empty.IsRestricted());
		CHECK(empty.Alpha().rows() == 0);
		CHECK(empty.Alpha().cols() == 0);
		CHECK(empty.Alpha().nonZeros() == 0);

		// Shape queries are the null view's whole surface: Eigen's sparse reductions
		// treat 0×0 as uninitialized ("you are using a non initialized matrix",
		// SparseRedux.h), so sum(), norm() — and with them the distance-rule
		// comparisons — assert on an Empty mapped-sparse spinful rather than
		// yielding 0. Availability-by-call: the dense-view Empty story (zero
		// reductions, valid comparison operand) does not carry over to sparse.
	}
}

TEST_CASE("Spinful eager per-component application")
{
	using SpinfulMatrix = SecChem::Spinful<MatrixXd>;
	using SpinfulVector = SecChem::Spinful<VectorXd>;

	const MatrixXd alpha{{2, 1}, {1, 3}};  // eigenvalues (5 ± √5) / 2, trace 5
	const MatrixXd beta{{6, 1}, {1, 2}};   // eigenvalues 4 ± √5, trace 8
	const MatrixXd identity = MatrixXd::Identity(2, 2);

	const auto eigendecompose = [](const MatrixXd& block)
	{
		Eigen::SelfAdjointEigenSolver<MatrixXd> solver{block};
		return std::make_tuple(std::move(solver).eigenvalues(), std::move(solver).eigenvectors());
	};

	SECTION("Self-adjoint EVD round trip — restricted")
	{
		const SpinfulMatrix fock{alpha};

		const auto result = SecChem::ApplyPerIndependentComponent(eigendecompose, fock);
		const auto& eps = std::get<0>(result);
		const auto& coefficients = std::get<1>(result);

		static_assert(std::is_same_v<std::decay_t<decltype(eps)>, SpinfulVector>);
		static_assert(std::is_same_v<std::decay_t<decltype(coefficients)>, SpinfulMatrix>);

		CHECK(eps.IndependentComponentCount() == 1);
		CHECK(coefficients.IndependentComponentCount() == 1);

		CHECK(eps.Alpha()(0) == Approx((5.0 - std::sqrt(5.0)) / 2.0));
		CHECK(eps.Alpha()(1) == Approx((5.0 + std::sqrt(5.0)) / 2.0));

		const MatrixXd reconstructed =
		        coefficients.Alpha() * eps.Alpha().asDiagonal() * coefficients.Alpha().transpose();
		CHECK(reconstructed.isApprox(alpha));
	}

	SECTION("Self-adjoint EVD round trip — unrestricted")
	{
		const SpinfulMatrix fock{alpha, beta};

		const auto result = SecChem::ApplyPerIndependentComponent(eigendecompose, fock);
		const auto& eps = std::get<0>(result);
		const auto& coefficients = std::get<1>(result);

		CHECK(eps.IndependentComponentCount() == 2);
		CHECK(coefficients.IndependentComponentCount() == 2);

		CHECK(eps.Alpha()(0) == Approx((5.0 - std::sqrt(5.0)) / 2.0));
		CHECK(eps.Beta()(0) == Approx(4.0 - std::sqrt(5.0)));
		CHECK(eps.Beta()(1) == Approx(4.0 + std::sqrt(5.0)));

		const MatrixXd reconstructedAlpha =
		        coefficients.Alpha() * eps.Alpha().asDiagonal() * coefficients.Alpha().transpose();
		const MatrixXd reconstructedBeta =
		        coefficients.Beta() * eps.Beta().asDiagonal() * coefficients.Beta().transpose();
		CHECK(reconstructedAlpha.isApprox(alpha));
		CHECK(reconstructedBeta.isApprox(beta));
	}

	SECTION("SVD three-tuple unzip — unrestricted")
	{
		const MatrixXd wide{{1, 2}, {0, 0}, {1, 3}};
		const MatrixXd other{{2, 0}, {1, 1}, {1, 0}};
		const SpinfulMatrix pair{wide, other};

		const auto result = SecChem::ApplyPerIndependentComponent(
		        [](const MatrixXd& block)
		        {
			        Eigen::JacobiSVD<MatrixXd> solver{block, Eigen::ComputeThinU | Eigen::ComputeThinV};
			        return std::tuple{solver.matrixU(), solver.singularValues(), solver.matrixV()};
		        },
		        pair);

		CHECK(std::get<0>(result).IndependentComponentCount() == 2);
		CHECK(std::get<1>(result).Beta().size() == 2);
		CHECK(std::get<2>(result).Alpha().cols() == 2);  // thin V of a 2×3 block

		const MatrixXd reconstructed = std::get<0>(result).Alpha() * std::get<1>(result).Alpha().asDiagonal()
		                               * std::get<2>(result).Alpha().transpose();
		CHECK(reconstructed.isApprox(wide));
	}

	SECTION("Generalized EVD — restricted overlap, unrestricted Fock")
	{
		const SpinfulMatrix overlap{identity};  // restricted: broadcasts to both spins
		const SpinfulMatrix fock{alpha, beta};

		const auto result = SecChem::ApplyPerIndependentComponent(
		        [](const MatrixXd& f, const MatrixXd& s)
		        {
			        Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXd> solver{f, s};
			        return std::tuple{solver.eigenvalues(), solver.eigenvectors()};
		        },
		        fock,
		        overlap);

		const auto& eps = std::get<0>(result);
		const auto& coefficients = std::get<1>(result);

		CHECK(eps.IndependentComponentCount() == 2);  // the unified count of {2, 1}
		CHECK(coefficients.IndependentComponentCount() == 2);

		// S = I reduces the generalized problem to the ordinary one: same eigenvalues
		// as the plain EVD of each block, S-orthogonality is plain orthogonality.
		CHECK(eps.Alpha()(1) == Approx((5.0 + std::sqrt(5.0)) / 2.0));
		CHECK(eps.Beta()(0) == Approx(4.0 - std::sqrt(5.0)));
		CHECK((fock.Alpha() * coefficients.Alpha()).isApprox(coefficients.Alpha() * eps.Alpha().asDiagonal()));
		CHECK((fock.Beta() * coefficients.Beta()).isApprox(coefficients.Beta() * eps.Beta().asDiagonal()));
		CHECK((coefficients.Alpha().transpose() * identity * coefficients.Alpha()).isApprox(identity));
	}

	SECTION("Restricted fusion — the callable runs exactly once per slot")
	{
		int restrictedCalls = 0;
		const SpinfulMatrix restricted{alpha};
		SecChem::ApplyPerIndependentComponent(
		        [&restrictedCalls](const MatrixXd& block)
		        {
			        ++restrictedCalls;
			        return block.trace();
		        },
		        restricted);
		CHECK(restrictedCalls == 1);

		int unrestrictedCalls = 0;
		const SpinfulMatrix fock{alpha, beta};
		SecChem::ApplyPerIndependentComponent(
		        [&unrestrictedCalls](const MatrixXd& block)
		        {
			        ++unrestrictedCalls;
			        return block.trace();
		        },
		        fock);
		CHECK(unrestrictedCalls == 2);

		int mixedCalls = 0;
		const SpinfulMatrix overlap{identity};
		const auto traces = SecChem::ApplyPerIndependentComponent(
		        [&mixedCalls](const MatrixXd& f, const MatrixXd& s)
		        {
			        ++mixedCalls;
			        return f.trace() + s.trace();
		        },
		        fock,
		        overlap);
		CHECK(mixedCalls == 2);  // the restricted overlap broadcasts, the count joins to 2
		CHECK(traces.Alpha() == 7.0);
		CHECK(traces.Beta() == 10.0);
		static_assert(std::is_same_v<std::decay_t<decltype(traces)>, SecChem::Spinful<double>>);
	}

	SECTION("Non-tuple return — a single Spinful")
	{
		const SpinfulMatrix fock{alpha, beta};
		const auto trace =
		        SecChem::ApplyPerIndependentComponent([](const MatrixXd& block) { return block.trace(); }, fock);

		static_assert(std::is_same_v<std::decay_t<decltype(trace)>, SecChem::Spinful<double>>);
		CHECK(trace.IndependentComponentCount() == 2);
		CHECK(trace.Alpha() == 5.0);
		CHECK(trace.Beta() == 8.0);
	}

	SECTION("Expression operand — eager, no dangling")
	{
		const SpinfulMatrix fock{alpha, beta};
		const SpinfulMatrix shift{identity};

		const auto trace = SecChem::ApplyPerIndependentComponent([](const MatrixXd& block) { return block.trace(); },
		                                                         fock + shift);
		CHECK(trace.IndependentComponentCount() == 2);  // the node's join of {2, 1}
		CHECK(trace.Alpha() == 7.0);
		CHECK(trace.Beta() == 10.0);
	}

	SECTION("Empty operand reads as restricted")
	{
		SpinfulMatrix empty;
		int calls = 0;

		const auto result = SecChem::ApplyPerIndependentComponent(
		        [&calls](const MatrixXd& block)
		        {
			        ++calls;
			        return block;
		        },
		        empty);

		CHECK(calls == 1);
		CHECK(empty.IsEmpty());  // const reads do not exit the Empty state

		static_assert(std::is_same_v<std::decay_t<decltype(result)>, SpinfulMatrix>);
		CHECK(result.IndependentComponentCount() == 1);
		CHECK(!result.IsEmpty());  // outputs are never Empty — restricted 0×0
		CHECK(result.Alpha().rows() == 0);
		CHECK(result.Alpha().cols() == 0);
	}

	SECTION("Scalar components")
	{
		using SpinfulScalar = SecChem::Spinful<double>;
		const SpinfulScalar energies{1.5, 2.5};

		const auto doubled = SecChem::ApplyPerIndependentComponent([](const double x) { return 2.0 * x; }, energies);
		static_assert(std::is_same_v<std::decay_t<decltype(doubled)>, SpinfulScalar>);
		CHECK(doubled.IndependentComponentCount() == 2);
		CHECK(doubled.Alpha() == 3.0);
		CHECK(doubled.Beta() == 5.0);

		const SpinfulScalar single{1.5};
		const auto square =
		        SecChem::ApplyPerIndependentComponent([](const double x) { return std::tuple{x, x * x}; }, single);
		static_assert(std::is_same_v<std::decay_t<decltype(square)>, std::tuple<SpinfulScalar, SpinfulScalar>>);
		CHECK(std::get<0>(square).IndependentComponentCount() == 1);
		CHECK(std::get<0>(square).Alpha() == 1.5);
		CHECK(std::get<1>(square).Alpha() == 2.25);
	}
}

TEST_CASE("Spinful in-place fill mutators")
{
	using SpinfulMatrix = SecChem::Spinful<MatrixXd>;
	using SpinfulScalar = SecChem::Spinful<double>;
	using ComplexScalar = std::complex<double>;
	using SpinfulComplex = SecChem::Spinful<ComplexScalar>;
	using MappedMatrix = Eigen::Map<MatrixXd>;
	using MappedSpinful = SecChem::Spinful<MappedMatrix>;

	SECTION("Dense matrix fills preserve count and shapes and chain")
	{
		SpinfulMatrix s{MatrixXd::Constant(2, 3, 5.0), MatrixXd::Constant(2, 3, 5.0)};

		REQUIRE(&s.SetZero() == &s);
		CHECK(s.IndependentComponentCount() == 2);
		CHECK(s.Alpha().rows() == 2);
		CHECK(s.Alpha().cols() == 3);
		CHECK(s.Alpha().isZero());
		CHECK(s.Beta().isZero());

		s.SetOnes();
		CHECK((s.Alpha().array() == 1.0).all());
		CHECK((s.Beta().array() == 1.0).all());

		s.SetConstant(3.5);
		CHECK((s.Alpha().array() == 3.5).all());
		CHECK((s.Beta().array() == 3.5).all());

		s.SetRandom();
		CHECK(s.Alpha().allFinite());
		CHECK(s.Beta().allFinite());
		CHECK((s.Alpha().array() >= -1.0).all());
		CHECK((s.Alpha().array() <= 1.0).all());
		CHECK((s.Beta().array() >= -1.0).all());
		CHECK((s.Beta().array() <= 1.0).all());
		CHECK(s.IndependentComponentCount() == 2);
		CHECK(s.Alpha().rows() == 2);
		CHECK(s.Alpha().cols() == 3);
	}

	SECTION("A restricted fill reads back through the broadcast")
	{
		SpinfulMatrix s{MatrixXd::Constant(2, 2, 5.0)};

		s.SetConstant(-1.25);
		CHECK(s.IsRestricted());
		CHECK(s.IndependentComponentCount() == 1);
		CHECK(s.Alpha().rows() == 2);
		CHECK((s.Alpha().array() == -1.25).all());
		CHECK((s.Beta().array() == -1.25).all());  // Beta() aliases Alpha() when restricted

		s.SetZero();
		CHECK(s.Alpha().isZero());
	}

	SECTION("An Empty owner fills its fillers and stays Empty")
	{
		SpinfulMatrix s;

		s.SetZero();
		CHECK(s.IsEmpty());

		s.SetConstant(2.0);
		CHECK(s.IsEmpty());
	}

	SECTION("Scalar components fill by assignment")
	{
		SpinfulScalar s{1.5, 2.5};

		s.SetZero();
		CHECK(s.Alpha() == 0.0);
		CHECK(s.Beta() == 0.0);

		s.SetOnes();
		CHECK(s.Alpha() == 1.0);
		CHECK(s.Beta() == 1.0);

		s.SetConstant(2.5);
		CHECK(s.Alpha() == 2.5);
		CHECK(s.Beta() == 2.5);
		CHECK(s.IndependentComponentCount() == 2);

		s.SetRandom();
		CHECK(s.Alpha() >= -1.0);
		CHECK(s.Alpha() <= 1.0);
		CHECK(s.Beta() >= -1.0);
		CHECK(s.Beta() <= 1.0);
	}

	SECTION("Complex scalar components")
	{
		SpinfulComplex s{ComplexScalar{1.0, 1.0}, ComplexScalar{2.0, 2.0}};

		const ComplexScalar value{1.0, 2.0};
		s.SetConstant(value);
		CHECK(s.Alpha() == value);
		CHECK(s.Beta() == value);

		s.SetZero();
		CHECK(s.Alpha() == ComplexScalar{});
		CHECK(s.Beta() == ComplexScalar{});

		s.SetOnes();
		CHECK(s.Alpha().real() == 1.0);
		CHECK(s.Alpha().imag() == 0.0);

		s.SetRandom();
		CHECK(std::isfinite(s.Alpha().real()));
		CHECK(std::isfinite(s.Alpha().imag()));
		CHECK(std::isfinite(s.Beta().real()));
		CHECK(std::isfinite(s.Beta().imag()));
	}

	SECTION("Mapped components fill write-through")
	{
		std::vector<double> alphaBuf(4, 1.0);
		std::vector<double> betaBuf(4, 2.0);
		MappedSpinful view{MappedMatrix{alphaBuf.data(), 2, 2}, MappedMatrix{betaBuf.data(), 2, 2}};

		view.SetZero();
		CHECK(alphaBuf[0] == 0.0);
		CHECK(alphaBuf[3] == 0.0);
		CHECK(betaBuf[2] == 0.0);

		view.SetConstant(7.0);
		CHECK(alphaBuf[2] == 7.0);
		CHECK(betaBuf[0] == 7.0);

		view.SetOnes();
		CHECK(alphaBuf[1] == 1.0);
		CHECK(betaBuf[3] == 1.0);
	}
}

TEST_CASE("Spinful component traits (compile-time)")
{
	CHECK(HasTrivialScalarTraits<double>);
	CHECK(HasTrivialScalarTraits<std::complex<double>>);
	CHECK(SecChem::Detail::IsEigenType<SumExpr>);
	CHECK(!SecChem::Detail::IsEigenType<NotEigen>);
	CHECK(!SecChem::Detail::ComponentTraits<ForeignMatrix>::NestByValue);
	CHECK(SecChem::Detail::ComponentTraits<ForeignSum<ForeignMatrix, ForeignMatrix>>::NestByValue);
}
