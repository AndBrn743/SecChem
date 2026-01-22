// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Andy Brown

#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <memory>
#include <numeric>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

#include <SecChem/Utility/IEquatableWithTolerance.hpp>

#include <SecChem/AzimuthalQuantumNumber.hpp>
#include <SecChem/Element.hpp>
#include <SecChem/System.hpp>

namespace SecChem::BasisSet::Gaussian
{
	class ContractedRadialOrbitalSet;

	class SemiLocalEcpTerm;

	class SemiLocalEcpChannel;

	template <typename>
	class AbstractAzimuthalShell;

	class AzimuthalShell;

	class AzimuthalShellSegmentView;

	namespace Detail
	{
		template <OwnershipSemantics Semantics>
		class BasisSetImpl;

		template <OwnershipSemantics Semantics>
		class BasisSetLibraryImpl;
	}  // namespace Detail

	using BasisSet = Detail::BasisSetImpl<OwnershipSemantics::Value>;

	using SharedBasisSet = Detail::BasisSetImpl<OwnershipSemantics::Reference>;

	using BasisSetLibrary = Detail::BasisSetLibraryImpl<OwnershipSemantics::Value>;

	using SharedBasisSetLibrary = Detail::BasisSetLibraryImpl<OwnershipSemantics::Reference>;
}  // namespace SecChem::BasisSet::Gaussian

template <>
struct SecUtility::Traits<SecChem::BasisSet::Gaussian::ContractedRadialOrbitalSet>
{
	static constexpr auto DefaultEqualityComparisonTolerance = 1e-15;
};

template <>
struct SecUtility::Traits<SecChem::BasisSet::Gaussian::SemiLocalEcpTerm>
{
	static constexpr auto DefaultEqualityComparisonTolerance = 1e-15;
};

template <>
struct SecUtility::Traits<SecChem::BasisSet::Gaussian::SemiLocalEcpChannel>
{
	static constexpr auto DefaultEqualityComparisonTolerance = 1e-15;
};

template <typename Derived>
struct SecUtility::Traits<SecChem::BasisSet::Gaussian::AbstractAzimuthalShell<Derived>>
{
	static constexpr auto DefaultEqualityComparisonTolerance = 1e-15;
};

template <>
struct SecUtility::Traits<SecChem::BasisSet::Gaussian::AzimuthalShell>
{
	static constexpr auto DefaultEqualityComparisonTolerance = 1e-15;
};

#define SECCHEM_GAUSSIAN_BASIS_SET_INTERNAL
// clang-format off
#include "Detail/Gaussian.ContractedRadialOrbitalSet.hpp"
#include "Detail/Gaussian.SemiLocalEcp.hpp"
#include "Detail/Gaussian.AzimuthalShell.hpp"
#include "Detail/Gaussian.BasisSet.hpp"
#include "Detail/Gaussian.MolecularBasisSet.hpp"
// clang-format on
#undef SECCHEM_GAUSSIAN_BASIS_SET_INTERNAL
