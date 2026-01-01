// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Andy Brown

#pragma once


#if __has_include(<SecUtility/UnitOfMeasurement.hpp>)
#include <SecUtility/UnitOfMeasurement.hpp>
#else
namespace SecUtility::UnitOfMeasurement
{
	constexpr auto BohrRadius2Angstrom = [](const auto& length) -> decltype(auto)
	{
		return length * 5.2917721090380e-1;
	};

	constexpr auto Angstrom2BohrRadius = [](const auto& length) -> decltype(auto)
	{
		return length / 5.2917721090380e-1;
	};

	constexpr auto Degree2Radian = [](const auto& deg) -> decltype(auto)
	{
		static constexpr auto factor = 3.14159265358979323846 / 180;
		return deg * factor;
	};

	constexpr auto Radian2Degree = [](const auto& deg) -> decltype(auto)
	{
		static constexpr auto factor = 180 / 3.14159265358979323846;
		return deg * factor;
	};
}  // namespace SecUtility::UnitOfMeasurement
#endif
