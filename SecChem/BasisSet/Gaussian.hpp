// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Andy Brown

#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <filesystem>
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

	class SemiLocalEcp;

	template <typename>
	class AbstractAngularMomentumBlock;

	class AngularMomentumBlock;

	class AngularMomentumBlockSegmentView;

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
struct SecUtility::Traits<SecChem::BasisSet::Gaussian::SemiLocalEcp>
{
	static constexpr auto DefaultEqualityComparisonTolerance = 1e-15;
};

template <typename Derived>
struct SecUtility::Traits<SecChem::BasisSet::Gaussian::AbstractAngularMomentumBlock<Derived>>
{
	static constexpr auto DefaultEqualityComparisonTolerance = 1e-15;
};

template <>
struct SecUtility::Traits<SecChem::BasisSet::Gaussian::AngularMomentumBlock>
{
	static constexpr auto DefaultEqualityComparisonTolerance = 1e-15;
};

#define SECCHEM_GAUSSIAN_BASIS_SET_INTERNAL
#include "Detail/Gaussian.ContractedRadialOrbitalSet.hpp"
#include "Detail/Gaussian.SemiLocalEcp.hpp"
#include "Detail/Gaussian.AngularMomentumBlock.hpp"
#include "Detail/Gaussian.BasisSet.hpp"
#include "Detail/Gaussian.MolecularBasisSet.hpp"
#undef SECCHEM_GAUSSIAN_BASIS_SET_INTERNAL
