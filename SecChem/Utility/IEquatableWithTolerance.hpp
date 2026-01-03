// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Andy Brown

#pragma once

#include <type_traits>


namespace SecUtility
{
	template <typename>
	struct Traits;

	template <typename Derived>
	class IEquatableWithTolerance
	{
		static constexpr auto DefaultEqualityComparisonTolerance = Traits<Derived>::DefaultEqualityComparisonTolerance;
		using Scalar = std::decay_t<decltype(DefaultEqualityComparisonTolerance)>;

	public:
		constexpr bool EqualsTo(const Derived& other,
							   const Scalar tolerance = DefaultEqualityComparisonTolerance) const noexcept
		{
			return static_cast<const Derived*>(this)->EqualsTo_Impl(other, tolerance);
		}

		constexpr bool NotEqualsTo(const Derived& other,
								  const Scalar tolerance = DefaultEqualityComparisonTolerance) const noexcept
		{
			return !static_cast<const Derived*>(this)->EqualsTo_Impl(other, tolerance);
		}

		constexpr bool operator==(const Derived& other) const noexcept
		{
			return OperatorEquals_Impl(other);
		}

		constexpr bool operator!=(const Derived& other) const noexcept
		{
			return !OperatorEquals_Impl(other);
		}

	private:
		friend Derived;
		constexpr IEquatableWithTolerance() noexcept = default;
		constexpr IEquatableWithTolerance(const IEquatableWithTolerance&) noexcept = default;
		constexpr IEquatableWithTolerance(IEquatableWithTolerance&) noexcept = default;
		constexpr IEquatableWithTolerance& operator=(const IEquatableWithTolerance&) noexcept = default;
		constexpr IEquatableWithTolerance& operator=(IEquatableWithTolerance&) noexcept = default;

#if __cplusplus >= 202002L
		constexpr
#endif
				~IEquatableWithTolerance() noexcept = default;

		bool EqualsTo_Impl(const Derived& other, Scalar tolerance) = delete;

		bool OperatorEquals_Impl(const Derived& other) const noexcept
		{
			return static_cast<const Derived*>(this)->EqualsTo_Impl(other, Scalar{0});
		}
	};
}  // namespace SecUtility
