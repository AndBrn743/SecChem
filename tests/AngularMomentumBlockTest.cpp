//
// Created by Andy on 12/25/2025.
//

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <Eigen/Dense>
#include <stdexcept>

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinRel;

#include <SecChem/BasisSet/Gaussian.hpp>
using namespace SecChem::BasisSet::Gaussian;
using SecChem::AzimuthalQuantumNumber;
using SecChem::ElectronicSubShell;

static ContractedRadialOrbitalSet MakeSimpleContractedSet(const Eigen::Index primitiveCount,
                                                          const Eigen::Index contractedCount)
{
	Eigen::VectorXd exponents(primitiveCount);
	for (Eigen::Index i = 0; i < primitiveCount; ++i)
	{
		exponents[i] = 0.5 + i;
	}

	Eigen::MatrixXd contractions = Eigen::MatrixXd::Random(primitiveCount, contractedCount);

	return ContractedRadialOrbitalSet{std::move(exponents), std::move(contractions)};
}

static SemiLocalEcp MakeSimpleSemiLocalEcp(const Eigen::Index electronCount, const Eigen::Index termCount)
{
	Eigen::VectorXd coefficients = Eigen::VectorXd::Random(termCount);
	Eigen::VectorXd rExponents = Eigen::VectorXd::Random(termCount);
	Eigen::VectorXd gaussianExponents = Eigen::VectorXd::Random(termCount);
	return {electronCount, coefficients, rExponents, gaussianExponents};
}

//--------------------------------------------------------------------------------------------------------------------//

TEST_CASE("AngularMomentumBlock basic construction without ECP", "[AngularMomentumBlock]")
{
	AzimuthalQuantumNumber l{1};  // p shell
	auto crs = MakeSimpleContractedSet(3, 2);

	AngularMomentumBlock block{l, crs};

	REQUIRE(block.HasOrbital());
	REQUIRE_FALSE(block.HasSemiLocalEcp());
	REQUIRE(block.PrimitiveShellCount() == 3);
	REQUIRE(block.ContractedShellCount() == 2);

	const auto primitiveShells = block.PrimitiveShells();
	REQUIRE(primitiveShells.size() == 3);
	REQUIRE(primitiveShells[0] == ElectronicSubShell{2, 1});
	REQUIRE(primitiveShells[1] == ElectronicSubShell{3, 1});
	REQUIRE(primitiveShells[2] == ElectronicSubShell{4, 1});

	const auto contractedShells = block.ContractedShells();
	REQUIRE(contractedShells.size() == 2);
	REQUIRE(contractedShells[0] == ElectronicSubShell{2, 1});
	REQUIRE(contractedShells[1] == ElectronicSubShell{3, 1});
}

TEST_CASE("AngularMomentumBlock construction with ECP", "[AngularMomentumBlock]")
{
	AzimuthalQuantumNumber l{2};  // d shell
	auto crs = MakeSimpleContractedSet(4, 1);
	SemiLocalEcp ecp = MakeSimpleSemiLocalEcp(14, 2);

	AngularMomentumBlock block{l, crs, ecp};

	REQUIRE(block.HasOrbital());
	REQUIRE(block.HasSemiLocalEcp());
}

TEST_CASE("SemiLocalEcp() throws if no ECP is present", "[AngularMomentumBlock][ECP]")
{
	AzimuthalQuantumNumber l{0};
	auto crs = MakeSimpleContractedSet(1, 1);

	AngularMomentumBlock block{l, crs};

	REQUIRE(block.HasOrbital());
	REQUIRE_FALSE(block.HasSemiLocalEcp());
	REQUIRE_THROWS_AS(block.SemiLocalEcp(), std::logic_error);
}

TEST_CASE("SemiLocalEcp() returns reference when present", "[AngularMomentumBlock][ECP]")
{
	AzimuthalQuantumNumber l{0};
	auto crs = MakeSimpleContractedSet(1, 1);
	SemiLocalEcp ecp = MakeSimpleSemiLocalEcp(14, 4);

	AngularMomentumBlock block{l, crs, ecp};

	REQUIRE(block.HasOrbital());
	REQUIRE(block.HasSemiLocalEcp());
}

TEST_CASE("AddOrOverrideSemiLocalEcp adds ECP when missing", "[AngularMomentumBlock][ECP]")
{
	AzimuthalQuantumNumber l{1};
	auto crs = MakeSimpleContractedSet(2, 1);

	AngularMomentumBlock block{l, crs};
	REQUIRE_FALSE(block.HasSemiLocalEcp());

	SemiLocalEcp ecp = MakeSimpleSemiLocalEcp(8, 1);
	block.AddOrOverrideSemiLocalEcp(ecp);

	REQUIRE(block.HasSemiLocalEcp());
	REQUIRE(block.SemiLocalEcp().Coefficients().size() == 1);

	block.AddOrOverrideSemiLocalEcp(MakeSimpleSemiLocalEcp(8, 3));
	REQUIRE(block.SemiLocalEcp().Coefficients().size() == 3);
}

TEST_CASE("AddOrOverrideSemiLocalEcp replaces existing ECP", "[AngularMomentumBlock][ECP]")
{
	AzimuthalQuantumNumber l{1};
	auto crs = MakeSimpleContractedSet(2, 1);

	SemiLocalEcp ecp1 = MakeSimpleSemiLocalEcp(42, 4);
	SemiLocalEcp ecp2 = MakeSimpleSemiLocalEcp(42, 3);

	AngularMomentumBlock block{l, crs, ecp1};
	block.AddOrOverrideSemiLocalEcp(ecp2);

	REQUIRE(block.HasSemiLocalEcp());
	REQUIRE(block.SemiLocalEcp().RExponents().size() == 3);
	// Optional: if equality is defined
	// REQUIRE(block.SemiLocalEcp() == ecp2);
}

TEST_CASE("Primitive and contracted orbital counts are consistent", "[AngularMomentumBlock][Counts]")
{
	AzimuthalQuantumNumber l{2};  // d shell
	auto crs = MakeSimpleContractedSet(4, 3);

	AngularMomentumBlock block{l, crs};

	const auto nCart = l.CartesianMagneticQuantumNumberCount();
	const auto nSph = l.MagneticQuantumNumberCount();

	REQUIRE(block.PrimitiveCartesianOrbitalCount() == block.PrimitiveShellCount() * nCart);

	REQUIRE(block.PrimitiveSphericalOrbitalCount() == block.PrimitiveShellCount() * nSph);

	REQUIRE(block.ContractedCartesianOrbitalCount() == block.ContractedShellCount() * nCart);

	REQUIRE(block.ContractedSphericalOrbitalCount() == block.ContractedShellCount() * nSph);
}

TEST_CASE("ExponentSet and ContractionSets are forwarded correctly", "[AngularMomentumBlock]")
{
	AzimuthalQuantumNumber l{0};
	auto crs = MakeSimpleContractedSet(6, 3);

	const Eigen::VectorXd& exponents = crs.ExponentSet();
	const Eigen::MatrixXd& contractions = crs.ContractionSets();

	AngularMomentumBlock block{l, crs};

	REQUIRE(block.ExponentSet().isApprox(exponents));
	REQUIRE(block.ContractionSets().isApprox(contractions));
}

TEST_CASE("AngularMomentumBlock::Concat should work")
{
	const auto crs0 = MakeSimpleContractedSet(6, 3);
	const auto crs1 = MakeSimpleContractedSet(3, 2);
	const auto ecp0 = MakeSimpleSemiLocalEcp(32, 3);
	const auto ecp1 = MakeSimpleSemiLocalEcp(32, 1);

	AngularMomentumBlock amb0 = {AzimuthalQuantumNumber::S, crs0, ecp0};
	AngularMomentumBlock amb1 = {AzimuthalQuantumNumber::S, crs1, ecp1};

	AngularMomentumBlock amb3 = {AzimuthalQuantumNumber::S, crs0};
	AngularMomentumBlock amb4 = {AzimuthalQuantumNumber::S, crs1};

	AngularMomentumBlock amb5 = {AzimuthalQuantumNumber::S, std::nullopt, ecp0};
	AngularMomentumBlock amb6 = {AzimuthalQuantumNumber::S, std::nullopt, ecp1};

	{
		auto ambs = {amb0, amb1};
		const auto amb = AngularMomentumBlock::Concat(ambs.begin(), ambs.end());
		REQUIRE(amb.ExponentSet().size() == 6 + 3);
		REQUIRE(amb.ContractedShellCount() == 3 + 2);
		REQUIRE(amb.HasSemiLocalEcp());
		REQUIRE(amb.SemiLocalEcp().Coefficients().size() == 3 + 1);
	}

	{
		auto ambs = {amb0};
		const auto amb = AngularMomentumBlock::Concat(ambs.begin(), ambs.end());
		REQUIRE(amb == amb0);
	}

	{
		auto ambs = {amb0, amb4};
		const auto amb = AngularMomentumBlock::Concat(ambs.begin(), ambs.end());
		REQUIRE(amb.ExponentSet().size() == 6 + 3);
		REQUIRE(amb.ContractedShellCount() == 3 + 2);
		REQUIRE(amb.HasSemiLocalEcp());
		REQUIRE(amb.SemiLocalEcp().Coefficients().size() == 3);

		REQUIRE(amb.SegmentCount() == 2);
		REQUIRE(amb.Segment(0).ExponentSet() == amb0.ExponentSet());
		REQUIRE(amb.Segment(0).ContractionSets() == amb0.ContractionSets());
		REQUIRE(amb.Segment(1).ExponentSet() == amb4.ExponentSet());
		REQUIRE(amb.Segment(1).ContractionSets() == amb4.ContractionSets());
	}

	{
		auto ambs = {amb3, amb4};
		const auto amb = AngularMomentumBlock::Concat(ambs.begin(), ambs.end());
		REQUIRE(amb.ExponentSet().size() == 6 + 3);
		REQUIRE(amb.ContractedShellCount() == 3 + 2);
		REQUIRE(!amb.HasSemiLocalEcp());
	}

	{
		auto ambs = {amb5, amb6};
		const auto amb = AngularMomentumBlock::Concat(ambs.begin(), ambs.end());
		REQUIRE_FALSE(amb.HasOrbital());
		REQUIRE(amb.HasSemiLocalEcp());
	}

	{
		auto ambs = {amb1, amb6};
		const auto amb = AngularMomentumBlock::Concat(ambs.begin(), ambs.end());
		REQUIRE(amb.HasOrbital());
		REQUIRE(amb.HasSemiLocalEcp());
	}

	{
		auto ambs = {amb3, amb6};
		const auto amb = AngularMomentumBlock::Concat(ambs.begin(), ambs.end());
		REQUIRE(amb.HasOrbital());
		REQUIRE(amb.HasSemiLocalEcp());
	}
}

TEST_CASE("AngularMomentumBlock::Concat should throw on empty input range")
{
	std::initializer_list<AngularMomentumBlock> ambs = {};
	REQUIRE_THROWS(AngularMomentumBlock::Concat(ambs.begin(), ambs.end()));
}

TEST_CASE("AngularMomentumBlock::Concat should refuse concat blocks of different angular momentum")
{
	const auto crs0 = MakeSimpleContractedSet(6, 3);
	const auto crs1 = MakeSimpleContractedSet(3, 2);
	const auto ecp0 = MakeSimpleSemiLocalEcp(67, 3);
	const auto ecp1 = MakeSimpleSemiLocalEcp(67, 1);

	AngularMomentumBlock amb0 = {AzimuthalQuantumNumber::S, crs0, ecp0};
	AngularMomentumBlock amb1 = {AzimuthalQuantumNumber::P, crs1, ecp1};

	auto ambs = {amb0, amb1};
	REQUIRE_THROWS(AngularMomentumBlock::Concat(ambs.begin(), ambs.end()));
}

TEST_CASE("Pure ECP AngularMomentumBlock should work", "[AngularMomentumBlock][ECP]")
{
	AzimuthalQuantumNumber l{1};
	SemiLocalEcp ecp = MakeSimpleSemiLocalEcp(8, 1);

	AngularMomentumBlock block{l, std::nullopt, ecp};
	REQUIRE_FALSE(block.HasOrbital());
	REQUIRE(block.HasSemiLocalEcp());

	block.AddOrOverrideSemiLocalEcp(ecp);

	const auto crs0 = MakeSimpleContractedSet(6, 3);
	block.AddOrOverrideContractedRadialOrbitalSet(crs0);
	REQUIRE(block.HasOrbital());
	REQUIRE(block.HasSemiLocalEcp());
}
