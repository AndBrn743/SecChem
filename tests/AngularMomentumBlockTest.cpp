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

#include "../SecChem/BasisSet/Gaussian.hpp"
using namespace SecChem::BasisSet::Gaussian;
using SecChem::AzimuthalQuantumNumber;

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

static SemiLocalEcp MakeSimpleSemiLocalEcp(const Eigen::Index termCount)
{
	Eigen::VectorXd coefficients = Eigen::VectorXd::Random(termCount);
	Eigen::VectorXd rExponents = Eigen::VectorXd::Random(termCount);
	Eigen::VectorXd gaussianExponents = Eigen::VectorXd::Random(termCount);
	return {coefficients, rExponents, gaussianExponents};
}

//--------------------------------------------------------------------------------------------------------------------//

TEST_CASE("AngularMomentumBlock basic construction without ECP", "[AngularMomentumBlock]")
{
	AzimuthalQuantumNumber l{1}; // p shell
	auto crs = MakeSimpleContractedSet(3, 2);

	AngularMomentumBlock block{l, crs};

	REQUIRE_FALSE(block.HasSemiLocalEcp());
	REQUIRE(block.PrimitiveShellCount() == 3);
	REQUIRE(block.ContractedShellCount() == 2);
}

TEST_CASE("AngularMomentumBlock construction with ECP", "[AngularMomentumBlock]")
{
	AzimuthalQuantumNumber l{2}; // d shell
	auto crs = MakeSimpleContractedSet(4, 1);
	SemiLocalEcp ecp = MakeSimpleSemiLocalEcp(2);

	AngularMomentumBlock block{l, crs, ecp};

	REQUIRE(block.HasSemiLocalEcp());
}

TEST_CASE("SemiLocalEcp() throws if no ECP is present", "[AngularMomentumBlock][ECP]")
{
	AzimuthalQuantumNumber l{0};
	auto crs = MakeSimpleContractedSet(1, 1);

	AngularMomentumBlock block{l, crs};

	REQUIRE_FALSE(block.HasSemiLocalEcp());
	REQUIRE_THROWS_AS(block.SemiLocalEcp(), std::logic_error);
}

TEST_CASE("SemiLocalEcp() returns reference when present", "[AngularMomentumBlock][ECP]")
{
	AzimuthalQuantumNumber l{0};
	auto crs = MakeSimpleContractedSet(1, 1);
	SemiLocalEcp ecp = MakeSimpleSemiLocalEcp(4);

	AngularMomentumBlock block{l, crs, ecp};

	REQUIRE(block.HasSemiLocalEcp());
	REQUIRE(&block.SemiLocalEcp() != nullptr);
}

TEST_CASE("OverrideOrAddSemiLocalEcp adds ECP when missing", "[AngularMomentumBlock][ECP]")
{
	AzimuthalQuantumNumber l{1};
	auto crs = MakeSimpleContractedSet(2, 1);

	AngularMomentumBlock block{l, crs};
	REQUIRE_FALSE(block.HasSemiLocalEcp());

	SemiLocalEcp ecp = MakeSimpleSemiLocalEcp(1);
	block.OverrideOrAddSemiLocalEcp(ecp);

	REQUIRE(block.HasSemiLocalEcp());
	REQUIRE(block.SemiLocalEcp().Coefficients().size() == 1);

	block.OverrideOrAddSemiLocalEcp(MakeSimpleSemiLocalEcp(3));
	REQUIRE(block.SemiLocalEcp().Coefficients().size() == 3);
}

TEST_CASE("OverrideOrAddSemiLocalEcp replaces existing ECP", "[AngularMomentumBlock][ECP]")
{
	AzimuthalQuantumNumber l{1};
	auto crs = MakeSimpleContractedSet(2, 1);

	SemiLocalEcp ecp1 = MakeSimpleSemiLocalEcp(4);
	SemiLocalEcp ecp2 = MakeSimpleSemiLocalEcp(3);

	AngularMomentumBlock block{l, crs, ecp1};
	block.OverrideOrAddSemiLocalEcp(ecp2);

	REQUIRE(block.HasSemiLocalEcp());
	REQUIRE(block.SemiLocalEcp().RExponents().size() == 3);
	// Optional: if equality is defined
	// REQUIRE(block.SemiLocalEcp() == ecp2);
}

TEST_CASE("Primitive and contracted orbital counts are consistent", "[AngularMomentumBlock][Counts]")
{
	AzimuthalQuantumNumber l{2}; // d shell
	auto crs = MakeSimpleContractedSet(4, 3);

	AngularMomentumBlock block{l, crs};

	const auto nCart = l.CartesianMagneticQuantumNumberCount();
	const auto nSph  = l.MagneticQuantumNumberCount();

	REQUIRE(block.PrimitiveCartesianOrbitalCount() ==
			block.PrimitiveShellCount() * nCart);

	REQUIRE(block.PrimitiveSphericalOrbitalCount() ==
			block.PrimitiveShellCount() * nSph);

	REQUIRE(block.ContractedCartesianOrbitalCount() ==
			block.ContractedShellCount() * nCart);

	REQUIRE(block.ContractedSphericalOrbitalCount() ==
			block.ContractedShellCount() * nSph);
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
