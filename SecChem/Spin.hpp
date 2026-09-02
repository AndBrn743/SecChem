// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Andy Brown

#pragma once

#include <array>

namespace SecChem
{
	/// <summary>
	/// Spin label of a diagonal block of a spin-resolved object.
	/// </summary>
	/// <remarks>
	/// Underlying type <c>bool</c>: Alpha == 0, Beta == 1 — the value is usable
	/// directly as a logical component index.
	/// </remarks>
	enum class Spin : bool
	{
		Alpha,
		Beta
	};

	/// <summary>
	/// The logical spin components, in iteration order {Alpha, Beta}.
	/// </summary>
	inline constexpr std::array<Spin, 2> AllSpins = {Spin::Alpha, Spin::Beta};
}  // namespace SecChem
