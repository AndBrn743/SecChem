// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Andy Brown

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
