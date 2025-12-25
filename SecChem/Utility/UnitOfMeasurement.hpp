//
// Created by Andy on 12/25/2025.
//

#pragma once


#if __has_include(<SecUtility/UnitOfMeasurement.hpp>)
#include <SecUtility/UnitOfMeasurement.hpp>
#else
namespace SecUtility::UnitOfMeasurement
{
	template <typename T>
	constexpr decltype(auto) BohrRadius2Angstrom(const T length)
	{
		return length * 5.2917721090380e-1;
	}
}  // namespace SecUtility::UnitOfMeasurement
#endif
