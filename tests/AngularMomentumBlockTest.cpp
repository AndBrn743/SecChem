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
using SecChem::ElectronicSubshell;

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

//--------------------------------------------------------------------------------------------------------------------//

TEST_CASE("AngularMomentumBlock basic construction without ECP", "[AngularMomentumBlock]")
{
	AzimuthalQuantumNumber l{1};  // p shell
	auto crs = MakeSimpleContractedSet(3, 2);

	AngularMomentumBlock block{l, crs};

	REQUIRE(block.IsNotEmpty());
	REQUIRE(block.PrimitiveShellCount() == 3);
	REQUIRE(block.ContractedShellCount() == 2);

	const auto primitiveShells = block.PrimitiveShells();
	REQUIRE(primitiveShells.size() == 3);
	REQUIRE(primitiveShells[0] == ElectronicSubshell{2, 1});
	REQUIRE(primitiveShells[1] == ElectronicSubshell{3, 1});
	REQUIRE(primitiveShells[2] == ElectronicSubshell{4, 1});

	const auto contractedShells = block.ContractedShells();
	REQUIRE(contractedShells.size() == 2);
	REQUIRE(contractedShells[0] == ElectronicSubshell{2, 1});
	REQUIRE(contractedShells[1] == ElectronicSubshell{3, 1});
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

TEST_CASE("AngularMomentumBlock::Concat should work", "[AngularMomentumBlock][Concat]")
{
	const auto crs0 = MakeSimpleContractedSet(6, 3);
	const auto crs1 = MakeSimpleContractedSet(3, 2);

	AngularMomentumBlock amb0 = {AzimuthalQuantumNumber::S, crs0};
	AngularMomentumBlock amb1 = {AzimuthalQuantumNumber::S, crs1};

	AngularMomentumBlock amb3 = {AzimuthalQuantumNumber::S, crs0};
	AngularMomentumBlock amb4 = {AzimuthalQuantumNumber::S, crs1};

	AngularMomentumBlock amb5 = {AzimuthalQuantumNumber::S, {}};
	AngularMomentumBlock amb6 = {AzimuthalQuantumNumber::S, {}};

	{
		auto ambs = {amb0, amb1};
		const auto amb = AngularMomentumBlock::Concat(ambs.begin(), ambs.end());
		REQUIRE(amb.ExponentSet().size() == 6 + 3);
		REQUIRE(amb.ContractedShellCount() == 3 + 2);
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
	}

	{
		auto ambs = {amb5, amb6};
		const auto amb = AngularMomentumBlock::Concat(ambs.begin(), ambs.end());
		REQUIRE_FALSE(amb.IsNotEmpty());
	}

	{
		auto ambs = {amb1, amb6};
		const auto amb = AngularMomentumBlock::Concat(ambs.begin(), ambs.end());
		REQUIRE(amb.IsNotEmpty());
	}

	{
		auto ambs = {amb3, amb6};
		const auto amb = AngularMomentumBlock::Concat(ambs.begin(), ambs.end());
		REQUIRE(amb.IsNotEmpty());
	}
}

TEST_CASE("AngularMomentumBlock::Concat should throw on empty input range", "[AngularMomentumBlock][Concat][Exception]")
{
	std::initializer_list<AngularMomentumBlock> ambs = {};
	REQUIRE_THROWS(AngularMomentumBlock::Concat(ambs.begin(), ambs.end()));
}

TEST_CASE("AngularMomentumBlock::Concat should refuse concat blocks of different angular momentum", "[AngularMomentumBlock][Concat][Validation]")
{
	const auto crs0 = MakeSimpleContractedSet(6, 3);
	const auto crs1 = MakeSimpleContractedSet(3, 2);

	AngularMomentumBlock amb0 = {AzimuthalQuantumNumber::S, crs0};
	AngularMomentumBlock amb1 = {AzimuthalQuantumNumber::P, crs1};

	auto ambs = {amb0, amb1};
	REQUIRE_THROWS(AngularMomentumBlock::Concat(ambs.begin(), ambs.end()));
}
