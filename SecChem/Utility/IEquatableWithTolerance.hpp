//
// Created by Andy on 12/26/2025.
//

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
		constexpr bool EqualTo(const Derived& other,
							   const Scalar tolerance = DefaultEqualityComparisonTolerance) const noexcept
		{
			return static_cast<const Derived*>(this)->EqualTo_Impl(other, tolerance);
		}

		constexpr bool NotEqualTo(const Derived& other,
								  const Scalar tolerance = DefaultEqualityComparisonTolerance) const noexcept
		{
			return !static_cast<const Derived*>(this)->EqualTo_Impl(other, tolerance);
		}

		constexpr bool operator==(const Derived& other) const noexcept
		{
			return static_cast<const Derived*>(this)->EqualTo_Impl(other, static_cast<Scalar>(0));
		}

		constexpr bool operator!=(const Derived& other) const noexcept
		{
			return !static_cast<const Derived*>(this)->EqualTo_Impl(other, static_cast<Scalar>(0));
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

		bool EqualTo_Impl(const Derived& other, Scalar tolerance) = delete;
	};
}  // namespace SecUtility
