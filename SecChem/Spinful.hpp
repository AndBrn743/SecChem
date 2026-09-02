// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Andy Brown

// ReSharper disable CppRedundantTypenameKeyword
// ReSharper disable CppTemplateArgumentsCanBeDeduced
#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <ostream>
#include <tuple>
#include <type_traits>
#include <utility>

#include <Eigen/Core>
#include <SecUtility/Diagnostic/Exception.hpp>
#include <SecUtility/Math/Core.hpp>
#include <SecUtility/Raw/Int.hpp>

#include <SecChem/Spin.hpp>
#include <SecChem/Utility/IEquatableWithTolerance.hpp>

namespace Eigen
{
	template <typename TScalar, int TOptions, typename TStorageIndex>
	class SparseMatrix;
}  // namespace Eigen

namespace SecChem
{
	using namespace SecUtility;

	template <typename TSpinless>
	class Spinful;

	template <typename Derived>
	class SpinfulExpr;

	template <typename TOp, typename TOperand>
	class SpinfulUnaryOp;
}  // namespace SecChem

namespace SecChem::Detail
{
	template <typename T>
	inline constexpr bool IsSpinfulType = std::is_base_of_v<SpinfulExpr<std::remove_reference_t<std::remove_cv_t<T>>>,
	                                                        std::remove_reference_t<std::remove_cv_t<T>>>;

	template <typename T>
	inline constexpr bool IsEigenType =
	        std::is_base_of_v<Eigen::EigenBase<std::remove_reference_t<std::remove_cv_t<T>>>,
	                          std::remove_reference_t<std::remove_cv_t<T>>>;

	template <typename>
	struct is_spinful_owner : std::false_type
	{
		/* NO CODE */
	};

	template <typename TSpinless>
	struct is_spinful_owner<Spinful<TSpinless>> : std::true_type
	{
		/* NO CODE */
	};

	constexpr int BroadcastIndex(const int index, const int count) noexcept
	{
		assert(count == 1 || count == 2);
		// ReSharper disable once CppRedundantParentheses
		return index & (count - 1);
	}

	template <typename T>
	struct DefaultValueOf
	{
		T operator()() const
		{
			return T{};
		}
	};

	template <typename TPlainObjectType, int MapOptions, typename TStrideType>
	struct DefaultValueOf<Eigen::Map<TPlainObjectType, MapOptions, TStrideType>>
	{
		Eigen::Map<TPlainObjectType, MapOptions, TStrideType> operator()() const
		{
			return {nullptr, 0, 0};
		}
	};

	template <typename TScalar, int TMatOptions, typename TStorageIndex, int TMapOptions, typename TStrideType>
	struct DefaultValueOf<
	        Eigen::Map<Eigen::SparseMatrix<TScalar, TMatOptions, TStorageIndex>, TMapOptions, TStrideType>>
	{
		Eigen::Map<Eigen::SparseMatrix<TScalar, TMatOptions, TStorageIndex>, TMapOptions, TStrideType> operator()()
		        const
		{
			static const TStorageIndex nullOuterIndex[1] = {0};
			return {0, 0, 0, const_cast<TStorageIndex*>(nullOuterIndex), nullptr, nullptr};
		}
	};

	template <typename TScalar, int TMatOptions, typename TStorageIndex, int TMapOptions, typename TStrideType>
	struct DefaultValueOf<
	        Eigen::Map<const Eigen::SparseMatrix<TScalar, TMatOptions, TStorageIndex>, TMapOptions, TStrideType>>
	{
		Eigen::Map<const Eigen::SparseMatrix<TScalar, TMatOptions, TStorageIndex>, TMapOptions, TStrideType>
		operator()() const
		{
			static const TStorageIndex nullOuterIndex[1] = {0};
			return {0, 0, 0, nullOuterIndex, nullptr, nullptr};
		}
	};

	template <typename TComponent>
	TComponent& AsComponent(TComponent& component) noexcept
	{
		return component;
	}

	template <typename TReferenced>
	TReferenced& AsComponent(const std::reference_wrapper<TReferenced>& component) noexcept
	{
		return component.get();
	}

	template <typename TReferenced>
	TReferenced& AsComponent(std::reference_wrapper<TReferenced>& component) noexcept
	{
		return component.get();
	}

	template <typename TSpinless>
	using ComponentOf = std::remove_reference_t<decltype(AsComponent(std::declval<TSpinless&>()))>;

	template <typename, typename = void>
	struct has_nest_by_value_member : std::false_type
	{
		/* NO CODE */
	};

	template <typename T>
	struct has_nest_by_value_member<T, std::void_t<decltype(T::NestByValue)>> : std::true_type
	{
		/* NO CODE */
	};

	template <typename T>
	constexpr bool ResolveNestByValue() noexcept
	{
		if constexpr (has_nest_by_value_member<T>::value)
		{
			return T::NestByValue;
		}
		else
		{
			return true;
		}
	}

	template <typename, typename = void>
	struct has_set_zero_rows_cols : std::false_type
	{
		/* NO CODE */
	};

	template <typename T>
	struct has_set_zero_rows_cols<T,
	                              std::void_t<decltype(std::declval<T&>().setZero(std::declval<Eigen::Index>(),
	                                                                              std::declval<Eigen::Index>()))>>
	    : std::true_type
	{
		/* NO CODE */
	};

	template <typename, typename = void>
	struct is_self_describing_component : std::false_type
	{
		/* NO CODE */
	};

	template <typename T>
	struct is_self_describing_component<T, std::void_t<typename T::Scalar, typename T::PlainObject>> : std::true_type
	{
		/* NO CODE */
	};

	template <typename T, typename = void>
	struct ComponentTraits;

	template <typename T>
	struct ComponentTraits<T, std::enable_if_t<IsEigenType<T>>>
	{
		using Scalar = typename Eigen::internal::traits<T>::Scalar;
		using Plain = typename Eigen::internal::plain_matrix_type<T>::type;
		static constexpr bool NestByValue = !(Eigen::internal::traits<T>::Flags & Eigen::NestByRefBit);
	};

	template <typename T>
	struct ComponentTraits<T, std::enable_if_t<!Detail::IsEigenType<T> && is_self_describing_component<T>::value>>
	{
		using Scalar = typename T::Scalar;
		using Plain = typename T::PlainObject;
		static constexpr bool NestByValue = Detail::ResolveNestByValue<T>();
	};

	template <typename T>
	struct ComponentTraits<std::reference_wrapper<T>> : ComponentTraits<T>
	{
		/* NO CODE */
	};

	template <>
	struct ComponentTraits<Int8>
	{
		using Scalar = Int8;
		using Plain = Int8;
		static constexpr bool NestByValue = true;
	};

	template <>
	struct ComponentTraits<Int16>
	{
		using Scalar = Int16;
		using Plain = Int16;
		static constexpr bool NestByValue = true;
	};

	template <>
	struct ComponentTraits<Int32>
	{
		using Scalar = Int32;
		using Plain = Int32;
		static constexpr bool NestByValue = true;
	};

	template <>
	struct ComponentTraits<Int64>
	{
		using Scalar = Int64;
		using Plain = Int64;
		static constexpr bool NestByValue = true;
	};

	template <>
	struct ComponentTraits<float>
	{
		using Scalar = float;
		using Plain = float;
		static constexpr bool NestByValue = true;
	};

	template <>
	struct ComponentTraits<double>
	{
		using Scalar = double;
		using Plain = double;
		static constexpr bool NestByValue = true;
	};

	template <>
	struct ComponentTraits<std::complex<float>>
	{
		using Scalar = std::complex<float>;
		using Plain = std::complex<float>;
		static constexpr bool NestByValue = true;
	};

	template <>
	struct ComponentTraits<std::complex<double>>
	{
		using Scalar = std::complex<double>;
		using Plain = std::complex<double>;
		static constexpr bool NestByValue = true;
	};
}  // namespace SecChem::Detail

template <typename TSpinless>
struct SecUtility::Traits<SecChem::Spinful<TSpinless>>
{
	using RealScalar = typename Eigen::NumTraits<typename SecChem::Detail::ComponentTraits<TSpinless>::Scalar>::Real;
	static constexpr auto DefaultEqualityComparisonTolerance = Eigen::NumTraits<RealScalar>::dummy_precision();
	using Scalar = typename SecChem::Detail::ComponentTraits<TSpinless>::Scalar;
	using SpinComponent = SecChem::Detail::ComponentOf<TSpinless>;
	using Evaluated = SecChem::Spinful<TSpinless>;
	static constexpr bool NestByValue = false;
};

template <typename TDerived>
struct SecUtility::Traits<SecChem::SpinfulExpr<TDerived>>
{
	using RealScalar = typename Eigen::NumTraits<typename Traits<TDerived>::Scalar>::Real;
	static constexpr auto DefaultEqualityComparisonTolerance = Eigen::NumTraits<RealScalar>::dummy_precision();
};

namespace SecChem::Detail
{
	template <typename, typename = void>
	struct is_component : std::false_type
	{
		/* NO CODE */
	};

	template <typename T>
	struct is_component<T, std::void_t<typename ComponentTraits<T>::Scalar>> : std::true_type
	{
		/* NO CODE */
	};

	template <typename E>
	using SpinfulOperandStorage = std::conditional_t<Traits<E>::NestByValue, E, const E&>;

	template <typename T>
	using ComponentOperandStorage = std::conditional_t<ComponentTraits<T>::NestByValue, T, const T&>;

	struct TraceOf
	{
		template <typename TComponent>
		auto operator()(const TComponent& component, const int) const -> decltype(component.trace())
		{
			return component.trace();
		}

		template <typename TComponent, typename = std::enable_if_t<!IsEigenType<TComponent>>>
		TComponent operator()(const TComponent& component, const long) const
		{
			return component;
		}
	};

	struct DeterminantOf
	{
		template <typename TComponent>
		auto operator()(const TComponent& component, const int) const -> decltype(component.determinant())
		{
			return component.determinant();
		}

		template <typename TComponent, typename = std::enable_if_t<!IsEigenType<TComponent>>>
		TComponent operator()(const TComponent& component, const long) const
		{
			return component;
		}
	};

	struct SumOf
	{
		template <typename TComponent>
		auto operator()(const TComponent& component, const int) const -> decltype(component.sum())
		{
			return component.sum();
		}

		template <typename TComponent>
		TComponent operator()(const TComponent& component, const long) const
		{
			return component;
		}
	};

	struct SquaredNormOf
	{
		template <typename TComponent>
		auto operator()(const TComponent& component, const int) const -> decltype(component.squaredNorm())
		{
			return component.squaredNorm();
		}

		template <typename TComponent>
		auto operator()(const TComponent& component, const long) const
		{
			if constexpr (std::is_arithmetic_v<TComponent>)
			{
				return component * component;
			}
			else
			{
				return std::norm(component);
			}
		}
	};

	struct MaxCoeffOf
	{
		template <typename TComponent>
		auto operator()(const TComponent& component, const int) const -> decltype(component.maxCoeff())
		{
			return component.maxCoeff();
		}

		template <typename TComponent>
		TComponent operator()(const TComponent& component, const long) const
		{
			return component;
		}
	};

	struct MinCoeffOf
	{
		template <typename TComponent>
		auto operator()(const TComponent& component, const int) const -> decltype(component.minCoeff())
		{
			return component.minCoeff();
		}

		template <typename TComponent>
		TComponent operator()(const TComponent& component, const long) const
		{
			return component;
		}
	};

	struct AllFiniteOf
	{
		template <typename TComponent>
		auto operator()(const TComponent& component, const int) const -> decltype(component.allFinite())
		{
			return component.allFinite();
		}

		template <typename TComponent>
		bool operator()(const TComponent& component, const long) const
		{
			if constexpr (std::is_arithmetic_v<TComponent>)
			{
				return std::isfinite(component);
			}
			else
			{
				return std::isfinite(component.real()) && std::isfinite(component.imag());
			}
		}
	};

	template <int P>
	struct LpNormOf
	{
		template <typename TComponent>
		auto operator()(const TComponent& component, const int) const
		        -> decltype(std::pow(component.template lpNorm<P>(), P))
		{
			return std::pow(component.template lpNorm<P>(), P);
		}

		template <typename TComponent>
		auto operator()(const TComponent& component, const long) const
		{
			return std::pow(std::abs(component), P);
		}
	};

	template <>
	struct LpNormOf<Eigen::Infinity>
	{
		template <typename TComponent>
		auto operator()(const TComponent& component, const int) const
		        -> decltype(component.template lpNorm<Eigen::Infinity>())
		{
			return component.template lpNorm<Eigen::Infinity>();
		}

		template <typename TComponent>
		auto operator()(const TComponent& component, const long) const
		{
			return std::abs(component);
		}
	};

	struct DistanceOf
	{
		template <typename TLhs, typename TRhs>
		auto operator()(const TLhs& lhs, const TRhs& rhs, const int) const -> decltype((lhs - rhs).norm())
		{
			return (lhs - rhs).norm();
		}

		template <typename TLhs, typename TRhs>
		auto operator()(const TLhs& lhs, const TRhs& rhs, const long) const
		{
			return std::abs(lhs - rhs);
		}
	};

	template <typename TScalar>
	struct CastTo
	{
		template <typename TComponent>
		auto operator()(const TComponent& component, const int) const -> decltype(component.template cast<TScalar>())
		{
			return component.template cast<TScalar>();
		}

		template <typename TComponent>
		TScalar operator()(const TComponent& component, const long) const
		{
			return static_cast<TScalar>(component);
		}
	};

	struct RealOf
	{
		template <typename TComponent>
		auto operator()(const TComponent& component, const int) const -> decltype(component.real())
		{
			return component.real();
		}

		template <typename TComponent>
		auto operator()(const TComponent& component, const long) const -> decltype(std::real(component))
		{
			return std::real(component);
		}
	};

	struct ImagOf
	{
		template <typename TComponent>
		auto operator()(const TComponent& component, const int) const -> decltype(component.imag())
		{
			return component.imag();
		}

		template <typename TComponent>
		auto operator()(const TComponent& component, const long) const -> decltype(std::imag(component))
		{
			return std::imag(component);
		}
	};

	template <typename TComponent>
	struct SpinComponentRef
	{
		// ReSharper disable once CppRedundantElaboratedTypeSpecifier
		enum Spin Spin;
		const TComponent& Component;
	};

	template <typename TComponent>
	class SpinComponentRange
	{
	public:
		using value_type = SpinComponentRef<TComponent>;
		using const_iterator = const value_type*;

		SpinComponentRange(const TComponent& alpha, const TComponent& beta, const int size) noexcept
		    : m_Entries{value_type{Spin::Alpha, alpha}, value_type{Spin::Beta, beta}}, m_Size(size)
		{
			/* NO CODE */
		}

		const_iterator begin() const noexcept
		{
			return m_Entries.data();
		}

		const_iterator end() const noexcept
		{
			return m_Entries.data() + m_Size;
		}

	private:
		std::array<value_type, 2> m_Entries;
		int m_Size;
	};

	template <typename TLhs, typename TRhs, typename TTolerance>
	bool IsLogicallyEqual(const SpinfulExpr<TLhs>& lhs, const SpinfulExpr<TRhs>& rhs, TTolerance tolerance) noexcept;

	template <typename>
	struct is_std_tuple : std::false_type
	{
		/* NO CODE */
	};

	template <typename... Ts>
	struct is_std_tuple<std::tuple<Ts...>> : std::true_type
	{
		/* NO CODE */
	};

	template <std::size_t... Indices, typename SpinDegeneratedResult>
	auto MakeSpinfuls(std::index_sequence<Indices...>, SpinDegeneratedResult spinDegeneratedResult)
	{
		return std::tuple{Spinful<std::tuple_element_t<Indices, SpinDegeneratedResult>>(
		        std::get<Indices>(std::move(spinDegeneratedResult)))...};
	}

	template <std::size_t... Indices, typename AlphaSpinResult, typename BetaSpinResult>
	auto MakeSpinfuls(std::index_sequence<Indices...>, AlphaSpinResult alphaSpinResult, BetaSpinResult betaSpinResult)
	{
		return std::tuple{Spinful<std::tuple_element_t<Indices, AlphaSpinResult>>(
		        std::get<Indices>(std::move(alphaSpinResult)), std::get<Indices>(std::move(betaSpinResult)))...};
	}

	template <typename SpinDegeneratedResult>
	auto AssembleSpinfulResults(SpinDegeneratedResult spinDegeneratedResult)
	{
		if constexpr (is_std_tuple<SpinDegeneratedResult>::value)
		{
			return MakeSpinfuls(std::make_index_sequence<std::tuple_size_v<SpinDegeneratedResult>>{},
			                    std::move(spinDegeneratedResult));
		}
		else
		{
			return Spinful<SpinDegeneratedResult>(std::move(spinDegeneratedResult));
		}
	}

	template <typename AlphaSpinResult, typename BetaSpinResult>
	auto AssembleSpinfulResults(AlphaSpinResult alphaSpinResult, BetaSpinResult betaSpinResult)
	{
		static_assert(std::is_same_v<AlphaSpinResult, BetaSpinResult>,
		              "ApplyPerIndependentComponent's callable must return the same type at both spin slots");

		if constexpr (is_std_tuple<AlphaSpinResult>::value)
		{
			return MakeSpinfuls(std::make_index_sequence<std::tuple_size_v<AlphaSpinResult>>{},
			                    std::move(alphaSpinResult),
			                    std::move(betaSpinResult));
		}
		else
		{
			return Spinful<AlphaSpinResult>(std::move(alphaSpinResult), std::move(betaSpinResult));
		}
	}
}  // namespace SecChem::Detail

namespace SecChem
{
	template <typename Derived>
	class SpinfulExpr : public IEquatableWithTolerance<SpinfulExpr<Derived>>
	{
	public:
		Derived& AsDerived() noexcept
		{
			return static_cast<Derived&>(*this);
		}

		const Derived& AsDerived() const noexcept
		{
			return static_cast<const Derived&>(*this);
		}

		using ToleranceScalar = std::decay_t<decltype(Traits<SpinfulExpr>::DefaultEqualityComparisonTolerance)>;

		/* CRTP VIRTUAL */ int IndependentComponentCount() const noexcept
		{
			const int count = AsDerived().IndependentComponentCount();
			assert(count == 1 || count == 2);
			return count;
		}

		/* CRTP VIRTUAL */ decltype(auto) AtIndex(const int index) noexcept
		{
			assert(index >= 0 && index < 2);
			return AsDerived().AtIndex(index);
		}

		/* CRTP VIRTUAL */ decltype(auto) AtIndex(const int index) const noexcept
		{
			assert(index >= 0 && index < 2);
			return AsDerived().AtIndex(index);
		}

		decltype(auto) operator[](const Spin spin) noexcept
		{
			return AsDerived().AtIndex(static_cast<int>(spin));
		}

		decltype(auto) operator[](const Spin spin) const noexcept
		{
			return AsDerived().AtIndex(static_cast<int>(spin));
		}

		typename Traits<Derived>::Evaluated Evaluate() const
		{
			return static_cast<typename Traits<Derived>::Evaluated>(AsDerived());
		}

		template <typename F>
		auto Componentwise(F&& f) const
		{
			return SpinfulUnaryOp<std::decay_t<F>, Derived>{AsDerived(), f};
		}

		auto Transpose() const
		{
			return Componentwise([](const auto& component) { return component.transpose(); });
		}

		auto Adjoint() const
		{
			return Componentwise([](const auto& component) { return component.adjoint(); });
		}

		auto Conjugate() const
		{
			return Componentwise([](const auto& component) { return component.conjugate(); });
		}

		auto Inverse() const
		{
			return Componentwise([](const auto& component) { return component.inverse(); });
		}

		auto AsDiagonal() const
		{
			return Componentwise([](const auto& component) { return component.asDiagonal(); });
		}

		auto Diagonal() const
		{
			return Componentwise([](const auto& component) { return component.diagonal(); });
		}

		template <typename TScalar>
		auto Cast() const
		{
			return Componentwise([](const auto& component) { return Detail::CastTo<TScalar>{}(component, 0); });
		}

		auto Re() const
		{
			return Componentwise([](const auto& component) { return Detail::RealOf{}(component, 0); });
		}

		auto Im() const
		{
			return Componentwise([](const auto& component) { return Detail::ImagOf{}(component, 0); });
		}

		auto Trace() const
		{
			return FoldComponents(Detail::TraceOf{}, std::plus<>{});
		}

		auto Determinant() const
		{
			return FoldComponents(Detail::DeterminantOf{}, std::multiplies<>{});
		}

		auto Sum() const
		{
			return FoldComponents(Detail::SumOf{}, std::plus<>{});
		}

		auto SquaredNorm() const
		{
			return FoldComponents(Detail::SquaredNormOf{}, std::plus<>{});
		}

		auto Norm() const
		{
			return std::sqrt(SquaredNorm());
		}

		template <int P>
		auto LpNorm() const
		{
			if constexpr (P == Eigen::Infinity)
			{
				return FoldComponents(Detail::LpNormOf<P>{}, Math::Max);
			}
			return std::pow(FoldComponents(Detail::LpNormOf<P>{}, std::plus<>{}), 1.0 / P);
		}

		auto MaxCoeff() const
		{
			return FoldComponents(Detail::MaxCoeffOf{}, Math::Max);
		}

		auto MinCoeff() const
		{
			return FoldComponents(Detail::MinCoeffOf{}, Math::Min);
		}

		bool IsFinite() const
		{
			return FoldComponents(Detail::AllFiniteOf{}, std::logical_and<>{});
		}

		typename Detail::ComponentTraits<typename Traits<Derived>::SpinComponent>::Plain SpinSum() const
		{
			using TPlain = typename Detail::ComponentTraits<typename Traits<Derived>::SpinComponent>::Plain;
			TPlain sum = AsDerived().AtIndex(0);
			if (AsDerived().IndependentComponentCount() == 1)
			{
				return typename Traits<Derived>::Scalar{2} * sum;
			}
			return sum + TPlain(AsDerived().AtIndex(1));
		}

		Spinful<typename Traits<Derived>::Scalar> SpinResolvedTrace() const
		{
			return SpinResolvedReduction(Detail::TraceOf{});
		}

		Spinful<typename Traits<Derived>::Scalar> SpinResolvedDeterminant() const
		{
			return SpinResolvedReduction(Detail::DeterminantOf{});
		}

		template <typename TOther>
		bool EqualsTo(const SpinfulExpr<TOther>& other,
		              const ToleranceScalar tolerance =
		                      Traits<SpinfulExpr>::DefaultEqualityComparisonTolerance) const noexcept
		{
			static_assert(std::is_floating_point_v<ToleranceScalar>,
			              "EqualsTo requires a floating-point component scalar");
			return Detail::IsLogicallyEqual(AsDerived(), other.AsDerived(), tolerance);
		}

	private:
		friend Derived;
		friend IEquatableWithTolerance<SpinfulExpr>;

		bool EqualsTo_Impl(const SpinfulExpr& other, const ToleranceScalar tolerance) const noexcept
		{
			static_assert(std::is_arithmetic_v<ToleranceScalar>,
			              "Spinful equality requires an arithmetic component scalar");
			return Detail::IsLogicallyEqual(AsDerived(), other.AsDerived(), tolerance);
		}

		SpinfulExpr() noexcept = default;
		SpinfulExpr(const SpinfulExpr&) noexcept = default;
		SpinfulExpr& operator=(const SpinfulExpr&) noexcept = default;
		SpinfulExpr(SpinfulExpr&&) noexcept = default;
		SpinfulExpr& operator=(SpinfulExpr&&) noexcept = default;
		~SpinfulExpr() noexcept = default;

		template <typename TPerComponent, typename TCombine>
		auto FoldComponents(TPerComponent perComponent, TCombine combine) const
		{
			if (AsDerived().IndependentComponentCount() == 1)
			{
				const auto x = perComponent(AsDerived().AtIndex(0), 0);
				return combine(x, x);
			}

			return combine(perComponent(AsDerived().AtIndex(0), 0), perComponent(AsDerived().AtIndex(1), 0));
		}

		template <typename TReductionOp>
		Spinful<typename Traits<Derived>::Scalar> SpinResolvedReduction(TReductionOp reductionOp) const
		{
			if (AsDerived().IndependentComponentCount() == 2)
			{
				return Spinful{reductionOp(AsDerived().AtIndex(0), 0), reductionOp(AsDerived().AtIndex(1), 0)};
			}
			return Spinful(reductionOp(AsDerived().AtIndex(0), 0));
		}
	};

	namespace Detail
	{
		template <typename TLhs, typename TRhs, typename TTolerance>
		bool IsLogicallyEqual(const SpinfulExpr<TLhs>& lhs,
		                      const SpinfulExpr<TRhs>& rhs,
		                      const TTolerance tolerance) noexcept
		{
			const bool a = DistanceOf{}(lhs.AtIndex(0), rhs.AtIndex(0), 0) <= tolerance;
			if (lhs.IndependentComponentCount() == 1 && rhs.IndependentComponentCount() == 1)
			{
				// Both broadcast slots address Alpha, the second distance would repeat the first.
				return a;
			}
			return a && DistanceOf{}(lhs.AtIndex(1), rhs.AtIndex(1), 0) <= tolerance;
		}
	}  // namespace Detail

	template <typename TOp, typename TOperand>
	class SpinfulUnaryOp : public SpinfulExpr<SpinfulUnaryOp<TOp, TOperand>>
	{
	public:
		explicit SpinfulUnaryOp(const TOperand& operand) noexcept : m_Operand(operand)
		{
			/* NO CODE */
		}

		SpinfulUnaryOp(const TOperand& operand, const TOp& op) noexcept : m_Op(op), m_Operand(operand)
		{
			/* NO CODE */
		}

		/* CRTP OVERRIDE */ int IndependentComponentCount() const noexcept
		{
			return m_Operand.IndependentComponentCount();
		}

		/* CRTP OVERRIDE */ auto AtIndex(const int index) const noexcept
		{
			return m_Op(m_Operand.AtIndex(index));
		}

	private:
		[[no_unique_address]] TOp m_Op;
		Detail::SpinfulOperandStorage<TOperand> m_Operand;
	};

	template <typename TOp, typename TLhs, typename TRhs>
	class SpinfulIntraBinaryOp : public SpinfulExpr<SpinfulIntraBinaryOp<TOp, TLhs, TRhs>>
	{
	public:
		SpinfulIntraBinaryOp(const TLhs& lhs, const TRhs& rhs) noexcept : m_Lhs(lhs), m_Rhs(rhs)
		{
			/* NO CODE */
		}

		/* CRTP OVERRIDE */ int IndependentComponentCount() const noexcept
		{
			const int lhsCount = m_Lhs.IndependentComponentCount();
			const int rhsCount = m_Rhs.IndependentComponentCount();
			return lhsCount > rhsCount ? lhsCount : rhsCount;
		}

		/* CRTP OVERRIDE */ auto AtIndex(const int index) const noexcept
		{
			return m_Op(m_Lhs.AtIndex(index), m_Rhs.AtIndex(index));
		}

	private:
		[[no_unique_address]] TOp m_Op;
		Detail::SpinfulOperandStorage<TLhs> m_Lhs;
		Detail::SpinfulOperandStorage<TRhs> m_Rhs;
	};

	template <typename TOp, typename TLhs, typename TRhs>
	class SpinfulInterBinaryOpComponentRhs : public SpinfulExpr<SpinfulInterBinaryOpComponentRhs<TOp, TLhs, TRhs>>
	{
	public:
		SpinfulInterBinaryOpComponentRhs(const TLhs& lhs, const TRhs& rhs) noexcept : m_Lhs(lhs), m_Rhs(rhs)
		{
			/* NO CODE */
		}

		/* CRTP OVERRIDE */ int IndependentComponentCount() const noexcept
		{
			return m_Lhs.IndependentComponentCount();
		}

		/* CRTP OVERRIDE */ auto AtIndex(const int index) const noexcept
		{
			return m_Op(m_Lhs.AtIndex(index), m_Rhs);
		}

	private:
		[[no_unique_address]] TOp m_Op;
		Detail::SpinfulOperandStorage<TLhs> m_Lhs;
		Detail::ComponentOperandStorage<TRhs> m_Rhs;
	};

	template <typename TOp, typename TLhs, typename TRhs>
	class SpinfulInterBinaryOpComponentLhs : public SpinfulExpr<SpinfulInterBinaryOpComponentLhs<TOp, TLhs, TRhs>>
	{
	public:
		SpinfulInterBinaryOpComponentLhs(const TLhs& lhs, const TRhs& rhs) noexcept : m_Lhs(lhs), m_Rhs(rhs)
		{
			/* NO CODE */
		}

		/* CRTP OVERRIDE */ int IndependentComponentCount() const noexcept
		{
			return m_Rhs.IndependentComponentCount();
		}

		/* CRTP OVERRIDE */ auto AtIndex(const int index) const noexcept
		{
			return m_Op(m_Lhs, m_Rhs.AtIndex(index));
		}

	private:
		[[no_unique_address]] TOp m_Op;
		Detail::ComponentOperandStorage<TLhs> m_Lhs;
		Detail::SpinfulOperandStorage<TRhs> m_Rhs;
	};
}  // namespace SecChem

template <typename TOp, typename TOperand>
struct SecUtility::Traits<SecChem::SpinfulUnaryOp<TOp, TOperand>>
{
	using SpinComponent = decltype(std::declval<const TOp&>()(
	        std::declval<const SecChem::Detail::SpinfulOperandStorage<TOperand>&>().AtIndex(0)));

	using Scalar = typename SecChem::Detail::ComponentTraits<SpinComponent>::Scalar;

	using Evaluated = SecChem::Spinful<typename SecChem::Detail::ComponentTraits<SpinComponent>::Plain>;

	static constexpr bool NestByValue = true;
};

template <typename TOp, typename TLhs, typename TRhs>
struct SecUtility::Traits<SecChem::SpinfulIntraBinaryOp<TOp, TLhs, TRhs>>
{
	using SpinComponent = decltype(std::declval<const TOp&>()(
	        std::declval<const SecChem::Detail::SpinfulOperandStorage<TLhs>&>().AtIndex(0),
	        std::declval<const SecChem::Detail::SpinfulOperandStorage<TRhs>&>().AtIndex(0)));

	using Scalar = std::common_type_t<typename Traits<TLhs>::Scalar, typename Traits<TRhs>::Scalar>;

	using Evaluated = SecChem::Spinful<typename SecChem::Detail::ComponentTraits<SpinComponent>::Plain>;

	static constexpr bool NestByValue = true;
};

template <typename TOp, typename TLhs, typename TRhs>
struct SecUtility::Traits<SecChem::SpinfulInterBinaryOpComponentRhs<TOp, TLhs, TRhs>>
{
	using SpinComponent = decltype(std::declval<const TOp&>()(
	        std::declval<const SecChem::Detail::SpinfulOperandStorage<TLhs>&>().AtIndex(0),
	        std::declval<const SecChem::Detail::ComponentOperandStorage<TRhs>&>()));

	using Scalar = typename Traits<TLhs>::Scalar;

	using Evaluated = SecChem::Spinful<typename SecChem::Detail::ComponentTraits<SpinComponent>::Plain>;

	static constexpr bool NestByValue = true;
};

template <typename TOp, typename TLhs, typename TRhs>
struct SecUtility::Traits<SecChem::SpinfulInterBinaryOpComponentLhs<TOp, TLhs, TRhs>>
{
	using SpinComponent = decltype(std::declval<const TOp&>()(
	        std::declval<const SecChem::Detail::ComponentOperandStorage<TLhs>&>(),
	        std::declval<const SecChem::Detail::SpinfulOperandStorage<TRhs>&>().AtIndex(0)));

	using Scalar = typename Traits<TRhs>::Scalar;

	using Evaluated = SecChem::Spinful<typename SecChem::Detail::ComponentTraits<SpinComponent>::Plain>;

	static constexpr bool NestByValue = true;
};

namespace SecChem
{
	template <typename TLhs, typename TRhs>
	SpinfulIntraBinaryOp<std::plus<>, TLhs, TRhs> operator+(const SpinfulExpr<TLhs>& lhs, const SpinfulExpr<TRhs>& rhs)
	{
		return SpinfulIntraBinaryOp<std::plus<>, TLhs, TRhs>{lhs.AsDerived(), rhs.AsDerived()};
	}

	template <typename TLhs, typename TRhs>
	SpinfulIntraBinaryOp<std::minus<>, TLhs, TRhs> operator-(const SpinfulExpr<TLhs>& lhs, const SpinfulExpr<TRhs>& rhs)
	{
		return SpinfulIntraBinaryOp<std::minus<>, TLhs, TRhs>{lhs.AsDerived(), rhs.AsDerived()};
	}

	template <typename TLhs, typename TRhs>
	SpinfulIntraBinaryOp<std::multiplies<>, TLhs, TRhs> operator*(const SpinfulExpr<TLhs>& lhs,
	                                                              const SpinfulExpr<TRhs>& rhs)
	{
		return SpinfulIntraBinaryOp<std::multiplies<>, TLhs, TRhs>{lhs.AsDerived(), rhs.AsDerived()};
	}

	template <typename TLhs, typename TRhs>
	SpinfulIntraBinaryOp<std::divides<>, TLhs, TRhs> operator/(const SpinfulExpr<TLhs>& lhs,
	                                                           const SpinfulExpr<TRhs>& rhs)
	{
		return SpinfulIntraBinaryOp<std::divides<>, TLhs, TRhs>{lhs.AsDerived(), rhs.AsDerived()};
	}

	template <typename TOperand>
	SpinfulUnaryOp<std::negate<>, TOperand> operator-(const SpinfulExpr<TOperand>& operand)
	{
		return SpinfulUnaryOp<std::negate<>, TOperand>{operand.AsDerived()};
	}

	template <typename TLhs,
	          typename TRhs,
	          typename = std::enable_if_t<!Detail::IsSpinfulType<TRhs> && Detail::is_component<TRhs>::value>>
	SpinfulInterBinaryOpComponentRhs<std::plus<>, TLhs, TRhs> operator+(const SpinfulExpr<TLhs>& lhs, const TRhs& rhs)
	{
		return SpinfulInterBinaryOpComponentRhs<std::plus<>, TLhs, TRhs>{lhs.AsDerived(), rhs};
	}

	template <typename TLhs,
	          typename TRhs,
	          typename = std::enable_if_t<!Detail::IsSpinfulType<TRhs> && Detail::is_component<TRhs>::value>>
	SpinfulInterBinaryOpComponentRhs<std::minus<>, TLhs, TRhs> operator-(const SpinfulExpr<TLhs>& lhs, const TRhs& rhs)
	{
		return SpinfulInterBinaryOpComponentRhs<std::minus<>, TLhs, TRhs>{lhs.AsDerived(), rhs};
	}

	template <typename TLhs,
	          typename TRhs,
	          typename = std::enable_if_t<!Detail::IsSpinfulType<TRhs> && Detail::is_component<TRhs>::value>>
	SpinfulInterBinaryOpComponentRhs<std::multiplies<>, TLhs, TRhs> operator*(const SpinfulExpr<TLhs>& lhs,
	                                                                          const TRhs& rhs)
	{
		return SpinfulInterBinaryOpComponentRhs<std::multiplies<>, TLhs, TRhs>{lhs.AsDerived(), rhs};
	}

	template <typename TLhs,
	          typename TRhs,
	          typename = std::enable_if_t<!Detail::IsSpinfulType<TRhs> && Detail::is_component<TRhs>::value>>
	SpinfulInterBinaryOpComponentRhs<std::divides<>, TLhs, TRhs> operator/(const SpinfulExpr<TLhs>& lhs,
	                                                                       const TRhs& rhs)
	{
		return SpinfulInterBinaryOpComponentRhs<std::divides<>, TLhs, TRhs>{lhs.AsDerived(), rhs};
	}

	template <typename TLhs,
	          typename TRhs,
	          typename = std::enable_if_t<!Detail::IsSpinfulType<TLhs> && Detail::is_component<TLhs>::value>>
	SpinfulInterBinaryOpComponentLhs<std::plus<>, TLhs, TRhs> operator+(const TLhs& lhs, const SpinfulExpr<TRhs>& rhs)
	{
		return SpinfulInterBinaryOpComponentLhs<std::plus<>, TLhs, TRhs>{lhs, rhs.AsDerived()};
	}

	template <typename TLhs,
	          typename TRhs,
	          typename = std::enable_if_t<!Detail::IsSpinfulType<TLhs> && Detail::is_component<TLhs>::value>>
	SpinfulInterBinaryOpComponentLhs<std::minus<>, TLhs, TRhs> operator-(const TLhs& lhs, const SpinfulExpr<TRhs>& rhs)
	{
		return SpinfulInterBinaryOpComponentLhs<std::minus<>, TLhs, TRhs>{lhs, rhs.AsDerived()};
	}

	template <typename TLhs,
	          typename TRhs,
	          typename = std::enable_if_t<!Detail::IsSpinfulType<TLhs> && Detail::is_component<TLhs>::value>>
	SpinfulInterBinaryOpComponentLhs<std::multiplies<>, TLhs, TRhs> operator*(const TLhs& lhs,
	                                                                          const SpinfulExpr<TRhs>& rhs)
	{
		return SpinfulInterBinaryOpComponentLhs<std::multiplies<>, TLhs, TRhs>{lhs, rhs.AsDerived()};
	}

	template <typename TLhs,
	          typename TRhs,
	          typename = std::enable_if_t<!Detail::IsSpinfulType<TLhs> && Detail::is_component<TLhs>::value>>
	SpinfulInterBinaryOpComponentLhs<std::divides<>, TLhs, TRhs> operator/(const TLhs& lhs,
	                                                                       const SpinfulExpr<TRhs>& rhs)
	{
		return SpinfulInterBinaryOpComponentLhs<std::divides<>, TLhs, TRhs>{lhs, rhs.AsDerived()};
	}

	template <typename TLhs, typename TRhs>
	bool operator==(const SpinfulExpr<TLhs>& lhs, const SpinfulExpr<TRhs>& rhs) noexcept
	{
		using TReal = typename Eigen::NumTraits<typename Traits<TLhs>::Scalar>::Real;
		static_assert(std::is_arithmetic_v<TReal>, "Spinful equality requires an arithmetic component scalar");
		return Detail::IsLogicallyEqual(lhs, rhs, TReal{0});
	}

	template <typename TLhs, typename TRhs>
	bool operator!=(const SpinfulExpr<TLhs>& lhs, const SpinfulExpr<TRhs>& rhs) noexcept
	{
		return !(lhs == rhs);
	}

	template <typename TSpinless>
	class Spinful : public SpinfulExpr<Spinful<TSpinless>>
	{
	public:
		using SpinComponent = Detail::ComponentOf<TSpinless>;

		using ToleranceScalar = std::decay_t<decltype(Traits<Spinful>::DefaultEqualityComparisonTolerance)>;

		Spinful() = default;

		explicit Spinful(const TSpinless& component)
		    : m_Data{component, Detail::DefaultValueOf<TSpinless>{}()}, m_Count(1)
		{
			/* NO CODE */
		}

		explicit Spinful(TSpinless&& component)
		    : m_Data{std::move(component), Detail::DefaultValueOf<TSpinless>{}()}, m_Count(1)
		{
			/* NO CODE */
		}

		Spinful(TSpinless alpha, TSpinless beta) : m_Data{std::move(alpha), std::move(beta)}, m_Count(2)
		{
			/* NO CODE */
		}

		template <typename Expr>
		/*IMPLICIT*/ Spinful(const SpinfulExpr<Expr>& expr)
		{
			Assign(expr);
		}

		Spinful(const Spinful&) = default;
		Spinful(Spinful&&) = default;
		Spinful& operator=(const Spinful&) = default;
		Spinful& operator=(Spinful&&) = default;
		~Spinful() noexcept = default;

		template <typename Expr>
		Spinful& operator=(const SpinfulExpr<Expr>& expr)
		{
			Assign(expr);
			return *this;
		}

		template <typename Expr>
		Spinful& operator+=(const SpinfulExpr<Expr>& expr)
		{
			return CompoundAssign(expr, [](auto& block, const auto& value) { block += value; });
		}

		template <typename Expr>
		Spinful& operator-=(const SpinfulExpr<Expr>& expr)
		{
			return CompoundAssign(expr, [](auto& block, const auto& value) { block -= value; });
		}

		template <typename Expr>
		Spinful& operator*=(const SpinfulExpr<Expr>& expr)
		{
			return CompoundAssign(expr, [](auto& block, const auto& value) { block *= value; });
		}

		template <typename Expr>
		Spinful& operator/=(const SpinfulExpr<Expr>& expr)
		{
			return CompoundAssign(expr, [](auto& block, const auto& value) { block /= value; });
		}

		Spinful& operator+=(const TSpinless& component)
		{
			return CompoundAssignValue(component, [](auto& block, const auto& value) { block += value; });
		}

		Spinful& operator-=(const TSpinless& component)
		{
			return CompoundAssignValue(component, [](auto& block, const auto& value) { block -= value; });
		}

		Spinful& operator*=(const TSpinless& component)
		{
			return CompoundAssignValue(component, [](auto& block, const auto& value) { block *= value; });
		}

		Spinful& operator/=(const TSpinless& component)
		{
			return CompoundAssignValue(component, [](auto& block, const auto& value) { block /= value; });
		}

		template <typename TScalar = typename Detail::ComponentTraits<TSpinless>::Scalar,
		          typename = std::enable_if_t<!std::is_same_v<TSpinless, TScalar> && !Detail::IsSpinfulType<TScalar>
		                                      && Detail::is_component<TScalar>::value>>
		Spinful& operator*=(const TScalar scalar)
		{
			return CompoundAssignValue(scalar, [](auto& block, const auto& value) { block *= value; });
		}

		template <typename TScalar = typename Detail::ComponentTraits<TSpinless>::Scalar,
		          typename = std::enable_if_t<!std::is_same_v<TSpinless, TScalar> && !Detail::IsSpinfulType<TScalar>
		                                      && Detail::is_component<TScalar>::value>>
		Spinful& operator/=(const TScalar scalar)
		{
			return CompoundAssignValue(scalar, [](auto& block, const auto& value) { block /= value; });
		}

		/* CRTP OVERRIDE */ int IndependentComponentCount() const noexcept
		{
			return m_Count == 0 ? 1 : m_Count;
		}

		bool IsRestricted() const noexcept
		{
			return m_Count == 1;
		}

		bool IsEmpty() const noexcept
		{
			return m_Count == 0;
		}

		/* CRTP OVERRIDE */ decltype(auto) AtIndex(const int index) noexcept
		{
			assert(index >= 0 && index < 2);
			if (m_Count == 0)
			{
				m_Count = 1;
			}
			return Detail::AsComponent(m_Data[Detail::BroadcastIndex(index, IndependentComponentCount())]);
		}

		/* CRTP OVERRIDE */ decltype(auto) AtIndex(const int index) const noexcept
		{
			assert(index >= 0 && index < 2);
			return Detail::AsComponent(m_Data[Detail::BroadcastIndex(index, IndependentComponentCount())]);
		}

		decltype(auto) Alpha() const noexcept
		{
			return Detail::AsComponent(m_Data[0]);
		}

		decltype(auto) Beta() const noexcept
		{
			return Detail::AsComponent(m_Data[Detail::BroadcastIndex(1, IndependentComponentCount())]);
		}

		decltype(auto) Alpha() noexcept
		{
			if (m_Count == 0)
			{
				m_Count = 1;
			}
			return Detail::AsComponent(m_Data[0]);
		}

		decltype(auto) Beta() noexcept
		{
			if (m_Count == 0)
			{
				m_Count = 1;
			}
			return Detail::AsComponent(m_Data[Detail::BroadcastIndex(1, IndependentComponentCount())]);
		}

		Spinful ToUnrestricted() const
		{
			if (m_Count == 2)
			{
				return *this;
			}
			return Spinful{m_Data[0], m_Data[0]};
		}

		Spinful ToRestricted(
		        const ToleranceScalar tolerance = Traits<Spinful>::DefaultEqualityComparisonTolerance) const
		{
			static_assert(std::is_floating_point_v<ToleranceScalar>,
			              "ToRestricted requires a floating-point component scalar");

			if (const TSpinless &alpha = m_Data[0],
			    beta = m_Data[Detail::BroadcastIndex(1, IndependentComponentCount())];
			    Detail::DistanceOf{}(alpha, beta, 0) > tolerance)
			{
				throw InvalidOperationException{"Cannot restrict a Spinful: the Alpha and Beta components "
				                                "differ by more than the tolerance"};
			}
			return Spinful{m_Data[0]};
		}

		bool IsNearlyRestricted(
		        const ToleranceScalar tolerance = Traits<Spinful>::DefaultEqualityComparisonTolerance) const
		{
			static_assert(std::is_floating_point_v<ToleranceScalar>,
			              "IsNearlyRestricted requires a floating-point component scalar");
			return Detail::DistanceOf{}(Alpha(), Beta(), 0) <= tolerance;
		}

		Spinful ToRestrictedAverage() const
		{
			static_assert(std::is_floating_point_v<ToleranceScalar>,
			              "ToRestrictedAverage requires a floating-point component scalar");
			if (IndependentComponentCount() == 1)
			{
				return *this;
			}

			using TScalar = typename Detail::ComponentTraits<TSpinless>::Scalar;
			return Spinful{(m_Data[0] + m_Data[1]) * TScalar{0.5}};
		}

		Spinful ToSpinFlipped() const
		{
			if (m_Count == 2)
			{
				return Spinful{m_Data[1], m_Data[0]};
			}
			return *this;
		}

		Spinful& FlipSpin()
		{
			if (m_Count == 2)
			{
				std::swap(m_Data[0], m_Data[1]);
			}
			return *this;
		}

		Spinful& SetConstant(const typename Detail::ComponentTraits<TSpinless>::Scalar value)
		{
			for (auto& block : m_Data)
			{
				decltype(auto) component = Detail::AsComponent(block);
				if constexpr (Detail::IsEigenType<decltype(component)>)
				{
					component.setConstant(value);
				}
				else
				{
					component = value;
				}
			}
			return *this;
		}

		Spinful& SetZero()
		{
			for (auto& block : m_Data)
			{
				decltype(auto) component = Detail::AsComponent(block);
				if constexpr (Detail::IsEigenType<decltype(component)>)
				{
					component.setZero();
				}
				else
				{
					component = typename Detail::ComponentTraits<TSpinless>::Scalar{0};
				}
			}
			return *this;
		}

		Spinful& SetOnes()
		{
			for (auto& block : m_Data)
			{
				decltype(auto) component = Detail::AsComponent(block);
				if constexpr (Detail::IsEigenType<decltype(component)>)
				{
					component.setOnes();
				}
				else
				{
					component = typename Detail::ComponentTraits<TSpinless>::Scalar{1};
				}
			}
			return *this;
		}

		Spinful& SetRandom()
		{
			for (auto& block : m_Data)
			{
				decltype(auto) component = Detail::AsComponent(block);
				if constexpr (Detail::IsEigenType<decltype(component)>)
				{
					component.setRandom();
				}
				else
				{
					component = Eigen::internal::random<typename Detail::ComponentTraits<TSpinless>::Scalar>();
				}
			}
			return *this;
		}

		Detail::SpinComponentRange<SpinComponent> LogicalComponents() const noexcept
		{
			return Detail::SpinComponentRange<SpinComponent>{AtIndex(0), AtIndex(1), 2};
		}

		Detail::SpinComponentRange<SpinComponent> IndependentComponents() const noexcept
		{
			return Detail::SpinComponentRange<SpinComponent>{AtIndex(0), AtIndex(1), IndependentComponentCount()};
		}

		template <typename U = TSpinless, typename = std::enable_if_t<Detail::has_set_zero_rows_cols<U>::value>>
		static Spinful Zero(const Eigen::Index rows, const Eigen::Index cols, const int count = 1)
		{
			assert(count == 1 || count == 2);
			TSpinless component;
			component.setZero(rows, cols);
			return count == 2 ? Spinful{component, component} : Spinful{component};
		}

		template <typename U = TSpinless, typename = std::enable_if_t<!Detail::has_set_zero_rows_cols<U>::value>>
		static Spinful Zero(const int count = 1)
		{
			assert(count == 1 || count == 2);
			return count == 2 ? Spinful{TSpinless{0}, TSpinless{0}} : Spinful{TSpinless{0}};
		}

		template <typename U = TSpinless, typename = std::enable_if_t<Detail::has_set_zero_rows_cols<U>::value>>
		static Spinful Ones(const Eigen::Index rows, const Eigen::Index cols, const int count = 1)
		{
			assert(count == 1 || count == 2);
			TSpinless component;
			component.setOnes(rows, cols);
			return count == 2 ? Spinful{component, component} : Spinful{component};
		}

		template <typename U = TSpinless, typename = std::enable_if_t<Detail::has_set_zero_rows_cols<U>::value>>
		static Spinful Identity(const Eigen::Index dimension, const int count = 1)
		{
			assert(count == 1 || count == 2);
			TSpinless component;
			component.setIdentity(dimension, dimension);
			return count == 2 ? Spinful{component, component} : Spinful{component};
		}

		template <typename U = TSpinless, typename = std::enable_if_t<!Detail::has_set_zero_rows_cols<U>::value>>
		static Spinful One(const int count = 1)
		{
			assert(count == 1 || count == 2);
			return count == 2 ? Spinful{TSpinless{1}, TSpinless{1}} : Spinful{TSpinless{1}};
		}

	private:
		template <typename TExpr, typename TCompoundOp>
		Spinful& CompoundAssign(const SpinfulExpr<TExpr>& expr, TCompoundOp&& compoundOp)
		{
			const TExpr& rhs = expr.AsDerived();
			const int sourceCount = rhs.IndependentComponentCount();

			if (m_Count == 0)
			{
				m_Count = sourceCount;
			}
			else if (m_Count == 1 && sourceCount == 2)
			{
				throw InvalidOperationException{"Cannot compound-assign a Spinful expression with independent "
				                                "component count 2 to a Spinful with independent "
				                                "component count 1"};
			}

			if (m_Count == 2)
			{
				compoundOp(Detail::AsComponent(m_Data[0]), rhs.AtIndex(0));
				compoundOp(Detail::AsComponent(m_Data[1]), rhs.AtIndex(1));
				return *this;
			}

			compoundOp(Detail::AsComponent(m_Data[0]), rhs.AtIndex(0));
			return *this;
		}

		template <typename TValue, typename TCompoundOp>
		Spinful& CompoundAssignValue(const TValue& value, TCompoundOp&& compoundOp)
		{
			if (m_Count == 0)
			{
				m_Count = 1;
			}

			compoundOp(Detail::AsComponent(m_Data[0]), value);
			if (m_Count == 2)
			{
				compoundOp(Detail::AsComponent(m_Data[1]), value);
			}
			return *this;
		}

		template <typename Expr>
		void Assign(const SpinfulExpr<Expr>& expr)
		{
			const Expr& rhs = expr.AsDerived();
			const int sourceCount = rhs.IndependentComponentCount();

			if (m_Count == 0)
			{
				m_Count = sourceCount;
			}
			else if (m_Count == 1 && sourceCount == 2)
			{
				throw InvalidOperationException{"Cannot assign a Spinful expression with independent "
				                                "component count 2 to a Spinful with independent "
				                                "component count 1"};
			}

			if (m_Count == 2)
			{
				Detail::AsComponent(m_Data[0]) = rhs.AtIndex(0);
				if (sourceCount == 2)
				{
					Detail::AsComponent(m_Data[1]) = rhs.AtIndex(1);
				}
				else
				{
					// copy beats re-evaluating the component expression.
					Detail::AsComponent(m_Data[1]) = Detail::AsComponent(m_Data[0]);
				}
				return;
			}

			Detail::AsComponent(m_Data[0]) = rhs.AtIndex(0);
		}

		TSpinless m_Data[2]{Detail::DefaultValueOf<TSpinless>{}(), Detail::DefaultValueOf<TSpinless>{}()};
		int m_Count = 0;
	};

	template <typename TExpr>
	std::ostream& operator<<(std::ostream& stream, const SpinfulExpr<TExpr>& expr)
	{
		const TExpr& derived = expr.AsDerived();
		if constexpr (Detail::is_spinful_owner<TExpr>::value)
		{
			if (derived.IsEmpty())
			{
				return stream << "(Empty):\n";
			}
		}
		if (derived.IndependentComponentCount() == 1)
		{
			return stream << "(Alpha|Beta):\n" << derived.AtIndex(0) << '\n';
		}
		return stream << "(Alpha):\n" << derived.AtIndex(0) << "\n(Beta):\n" << derived.AtIndex(1) << '\n';
	}

	template <typename TFunc, typename... TExprs>
	auto ApplyPerIndependentComponent(TFunc&& f, const SpinfulExpr<TExprs>&... exprs)
	{
		static_assert(sizeof...(TExprs) >= 1, "ApplyPerIndependentComponent requires at least one spinful operand");

		const int unifiedCount = Math::Max(1, exprs.IndependentComponentCount()...);

		auto alphaSpinResult = f(exprs.AsDerived().AtIndex(0)...);
		if (unifiedCount == 2)
		{
			auto betaSpinResult = f(exprs.AsDerived().AtIndex(1)...);
			return Detail::AssembleSpinfulResults(std::move(alphaSpinResult), std::move(betaSpinResult));
		}
		return Detail::AssembleSpinfulResults(std::move(alphaSpinResult));
	}
}  // namespace SecChem
