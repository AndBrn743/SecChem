//
// Created by Andy on 1/23/2026.
//

#pragma once
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

namespace SecUtility
{
	template <typename>
	struct Traits;
}

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
