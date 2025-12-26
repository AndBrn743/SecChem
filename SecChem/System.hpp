//
// Created by Andy on 12/26/2025.
//

#pragma once


namespace SecChem
{
	using Scalar = double;

	enum class OwnershipSemantics : bool
	{
		Value,
		Reference
	};

	constexpr OwnershipSemantics AlternativeOf(const OwnershipSemantics semantics) noexcept
	{
		return static_cast<OwnershipSemantics>(!static_cast<bool>(semantics));
	}

	template <typename T>
	class Builder;
}  // namespace SecChem
