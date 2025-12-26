//
// Created by Andy on 12/25/2025.
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <stdexcept>

using Catch::Matchers::WithinRel;
using Catch::Approx;

#include "../SecChem/BasisSet/Gaussian.hpp"
using namespace SecChem::BasisSet::Gaussian;

TEST_CASE("ContractedRadialOrbitalSet: basic construction", "[ContractedRadialOrbitalSet]")
{
    Eigen::VectorXd exponents(3);
    exponents << 1.0, 0.5, 0.25;

    Eigen::MatrixXd contractions(3, 2);
    contractions <<
        0.7,  0.0,
        0.3,  0.8,
        0.0,  0.2;

    ContractedRadialOrbitalSet set(exponents, contractions);

    REQUIRE(set.PrimitiveShellCount() == 3);
    REQUIRE(set.ContractedShellCount() == 2);

    REQUIRE(set.ExponentSet().size() == 3);
    REQUIRE(set.ContractionSets().rows() == 3);
    REQUIRE(set.ContractionSets().cols() == 2);
}

TEST_CASE("ContractionSetView trims leading and trailing zeros", "[ContractedRadialOrbitalSet]")
{
    Eigen::VectorXd exponents(5);
    exponents << 10.0, 5.0, 2.0, 1.0, 0.5;

    Eigen::MatrixXd contractions(5, 1);
    contractions <<
        0.0,
        0.0,
        0.6,
        0.4,
        0.0;

    ContractedRadialOrbitalSet set(exponents, contractions);

    const auto view = set.ContractionSet(0);

    REQUIRE(view.Exponents().size() == 2);
    REQUIRE(view.ContractionCoefficients().size() == 2);

    CHECK_THAT(view.Exponent(0), WithinRel(2.0, 1e-14));
    CHECK_THAT(view.Exponent(1), WithinRel(1.0, 1e-14));

    CHECK_THAT(view.ContractionCoefficient(0), WithinRel(0.6, 1e-14));
    CHECK_THAT(view.ContractionCoefficient(1), WithinRel(0.4, 1e-14));
}

TEST_CASE("Each contracted shell has its own trimming window", "[ContractedRadialOrbitalSet]")
{
    Eigen::VectorXd exponents(4);
    exponents << 4.0, 2.0, 1.0, 0.5;

    Eigen::MatrixXd contractions(4, 2);
    contractions <<
        0.0,  0.2,
        0.5,  0.3,
        0.5,  0.0,
        0.0,  0.0;

    ContractedRadialOrbitalSet set(exponents, contractions);

    const auto view0 = set.ContractionSet(0);
    REQUIRE(view0.Exponents().size() == 2);
    CHECK_THAT(view0.Exponent(0), WithinRel(2.0, 1e-14));
    CHECK_THAT(view0.Exponent(1), WithinRel(1.0, 1e-14));

    const auto view1 = set.ContractionSet(1);
    REQUIRE(view1.Exponents().size() == 2);
    CHECK_THAT(view1.Exponent(0), WithinRel(4.0, 1e-14));
    CHECK_THAT(view1.Exponent(1), WithinRel(2.0, 1e-14));
}

TEST_CASE("ContractionSetView is a true view (no copying)", "[ContractedRadialOrbitalSet]")
{
    Eigen::VectorXd exponents(3);
    exponents << 3.0, 2.0, 1.0;

    Eigen::MatrixXd contractions(3, 1);
    contractions <<
        0.2,
        0.3,
        0.5;

    ContractedRadialOrbitalSet set(exponents, contractions);
    const auto view = set.ContractionSet(0);

    // Compare addresses to ensure Map points into original storage
    REQUIRE(view.Exponents().data() == set.ExponentSet().data());
    REQUIRE(view.ContractionCoefficients().data() == set.ContractionSets().col(0).data());
}

TEST_CASE("AccumulationCoefficientSet returns a row view", "[ContractedRadialOrbitalSet]")
{
    Eigen::VectorXd exponents(3);
    exponents << 3.0, 2.0, 1.0;

    Eigen::MatrixXd contractions(3, 2);
    contractions <<
        0.1, 0.2,
        0.3, 0.4,
        0.5, 0.6;

    ContractedRadialOrbitalSet set(exponents, contractions);

    const auto row1 = set.AccumulationCoefficientSet(1);
    REQUIRE(row1.size() == 2);

    CHECK_THAT(row1(0), WithinRel(0.3, 1e-14));
    CHECK_THAT(row1(1), WithinRel(0.4, 1e-14));
}

TEST_CASE("All-zero contraction set throws", "[ContractedRadialOrbitalSet][exception]")
{
    Eigen::VectorXd exponents(3);
    exponents << 1.0, 0.5, 0.25;

    Eigen::MatrixXd contractions(3, 1);
    contractions <<
        0.0,
        0.0,
        0.0;

    REQUIRE_THROWS_AS(
        ContractedRadialOrbitalSet(exponents, contractions),
        std::runtime_error
    );
}

TEST_CASE("Near-zero coefficients below threshold are ignored", "[ContractedRadialOrbitalSet]")
{
    Eigen::VectorXd exponents(3);
    exponents << 3.0, 2.0, 1.0;

    Eigen::MatrixXd contractions(3, 1);
    contractions <<
        1e-16,
        0.8,
        1e-16;

    ContractedRadialOrbitalSet set(exponents, contractions);
    const auto view = set.ContractionSet(0);

    REQUIRE(view.Exponents().size() == 1);
    CHECK_THAT(view.Exponent(0), WithinRel(2.0, 1e-14));
    CHECK_THAT(view.ContractionCoefficient(0), WithinRel(0.8, 1e-14));
}

TEST_CASE("Single primitive single contraction works", "[ContractedRadialOrbitalSet][edge]")
{
	Eigen::VectorXd exponents(1);
	exponents << 1.234;

	Eigen::MatrixXd contractions(1, 1);
	contractions << 0.9;

	ContractedRadialOrbitalSet set(exponents, contractions);
	const auto view = set.ContractionSet(0);

	REQUIRE(view.Exponents().size() == 1);
	CHECK_THAT(view.Exponent(0), WithinRel(1.234, 1e-14));
	CHECK_THAT(view.ContractionCoefficient(0), WithinRel(0.9, 1e-14));
}

TEST_CASE("Non-zero only at first primitive", "[ContractedRadialOrbitalSet][edge]")
{
    Eigen::VectorXd exponents(4);
    exponents << 4.0, 3.0, 2.0, 1.0;

    Eigen::MatrixXd contractions(4, 1);
    contractions <<
        0.75,
        0.0,
        0.0,
        0.0;

    ContractedRadialOrbitalSet set(exponents, contractions);
    const auto view = set.ContractionSet(0);

    REQUIRE(view.Exponents().size() == 1);
    CHECK_THAT(view.Exponent(0), WithinRel(4.0, 1e-14));
    CHECK_THAT(view.ContractionCoefficient(0), WithinRel(0.75, 1e-14));
}

TEST_CASE("Non-zero only at last primitive", "[ContractedRadialOrbitalSet][edge]")
{
	Eigen::VectorXd exponents(4);
	exponents << 4.0, 3.0, 2.0, 1.0;

	Eigen::MatrixXd contractions(4, 1);
	contractions <<
		0.0,
		0.0,
		0.0,
		0.42;

	ContractedRadialOrbitalSet set(exponents, contractions);
	const auto view = set.ContractionSet(0);

	REQUIRE(view.Exponents().size() == 1);
	CHECK_THAT(view.Exponent(0), WithinRel(1.0, 1e-14));
	CHECK_THAT(view.ContractionCoefficient(0), WithinRel(0.42, 1e-14));
}

TEST_CASE("Negative and positive coefficients are both accepted", "[ContractedRadialOrbitalSet][edge]")
{
	Eigen::VectorXd exponents(5);
	exponents << 5.0, 4.0, 3.0, 2.0, 1.0;

	Eigen::MatrixXd contractions(5, 1);
	contractions <<
		0.0,
		-0.2,
		0.3,
		-0.1,
		0.0;

	ContractedRadialOrbitalSet set(exponents, contractions);
	const auto view = set.ContractionSet(0);

	REQUIRE(view.Exponents().size() == 3);

	CHECK(view.Exponent(0) == Approx(4.0));
	CHECK(view.Exponent(2) == Approx(2.0));

	CHECK(view.ContractionCoefficient(0) == Approx(-0.2));
	CHECK(view.ContractionCoefficient(2) == Approx(-0.1));
}

TEST_CASE("Multiple contraction sets with distinct trimming windows", "[ContractedRadialOrbitalSet][edge]")
{
	Eigen::VectorXd exponents(5);
	exponents << 5, 4, 3, 2, 1;

	Eigen::MatrixXd contractions(5, 3);
	contractions <<
		0.0,  0.9,  0.0,
		0.5,  0.1,  0.0,
		0.5,  0.0,  0.3,
		0.0,  0.0,  0.4,
		0.0,  0.0,  0.0;

	ContractedRadialOrbitalSet set(exponents, contractions);

	const auto v0 = set.ContractionSet(0);
	REQUIRE(v0.Exponents().size() == 2); // 4,3

	const auto v1 = set.ContractionSet(1);
	REQUIRE(v1.Exponents().size() == 2); // 5,4

	const auto v2 = set.ContractionSet(2);
	REQUIRE(v2.Exponents().size() == 2); // 3,2
}

static ContractedRadialOrbitalSet MakeSimpleSet(const std::vector<double>& exponents,
												const Eigen::MatrixXd& contractions)
{
	Eigen::VectorXd e(exponents.size());
	for (Eigen::Index i = 0; i < e.size(); ++i)
	{
		e[i] = exponents[i];
	}

	return ContractedRadialOrbitalSet(e, contractions);
}

TEST_CASE("ContractedRadialOrbitalSet::operator== exact equality")
{
	Eigen::MatrixXd c(3, 2);
	c << 1.0, 0.0,
		 0.5, 0.3,
		 0.0, 0.7;

	auto a = MakeSimpleSet({1.0, 2.0, 3.0}, c);
	auto b = MakeSimpleSet({1.0, 2.0, 3.0}, c);
	auto a2 = a;

	REQUIRE(a == b);
	REQUIRE_FALSE(a != b);
	REQUIRE(a == a2);
	REQUIRE_FALSE(a != a2);
}

TEST_CASE("ContractedRadialOrbitalSet::operator== detects exponent difference")
{
	Eigen::MatrixXd c = Eigen::MatrixXd::Identity(2, 1);

	auto a = MakeSimpleSet({1.0, 2.0}, c);
	auto b = MakeSimpleSet({1.0, 2.0000001}, c);

	REQUIRE(a != b);
	REQUIRE_FALSE(a == b);
}

TEST_CASE("ContractedRadialOrbitalSet::operator== detects contraction difference")
{
	Eigen::MatrixXd c1(2, 1);
	c1 << 1.0, 0.5;

	Eigen::MatrixXd c2(2, 1);
	c2 << 1.0, 0.5000001;

	auto a = MakeSimpleSet({1.0, 2.0}, c1);
	auto b = MakeSimpleSet({1.0, 2.0}, c2);

	REQUIRE(a != b);
}

TEST_CASE("ContractedRadialOrbitalSet::EqualTo respects tolerance")
{
	Eigen::MatrixXd c1(3, 1);
	c1 << 1.0, 0.5, 0.25;

	Eigen::MatrixXd c2 = c1;
	c2(1, 0) += 1e-10;

	auto a = MakeSimpleSet({1.0, 2.0, 3.0}, c1);
	auto b = MakeSimpleSet({1.0, 2.0 + 5e-11, 3.0}, c2);

	REQUIRE(a.EqualTo(b, 1e-9));
	REQUIRE_FALSE(a.EqualTo(b, 1e-12));
}

TEST_CASE("ContractedRadialOrbitalSet::EqualTo detects shape mismatch")
{
	Eigen::MatrixXd c1(2, 1);
	c1 << 1.0, 0.5;

	Eigen::MatrixXd c2(2, 2);
	c2 << 1.0, 0.0,
		  0.5, 1.0;

	auto a = MakeSimpleSet({1.0, 2.0}, c1);
	auto b = MakeSimpleSet({1.0, 2.0}, c2);

	REQUIRE_FALSE(a.EqualTo(b, 1e-6));
	REQUIRE(a.NotEqualTo(b, 1e-6));
}

TEST_CASE("ContractedRadialOrbitalSet::EqualTo is strict with tolerance boundary")
{
	static constexpr double Delta = 0.00390625;  // 2^(-8), it can be expressed exactly with binary float

	Eigen::MatrixXd c1(1, 1);
	c1 << 1.0;

	Eigen::MatrixXd c2(1, 1);
	c2 << 1.0 + Delta;

	auto a = MakeSimpleSet({1.0}, c1);
	auto b = MakeSimpleSet({1.0 - Delta}, c2);

	REQUIRE_FALSE(a.EqualTo(b, Delta - 2e-16));
	REQUIRE(a.EqualTo(b, Delta));
	REQUIRE(a.EqualTo(b, Delta + 2e-16));
}

TEST_CASE("ContractedRadialOrbitalSet::EqualTo is order-sensitive")
{
	Eigen::MatrixXd c(2, 1);
	c << 1.0, 0.5;

	auto a = MakeSimpleSet({1.0, 2.0}, c);
	auto b = MakeSimpleSet({2.0, 1.0}, c);

	REQUIRE_FALSE(a.EqualTo(b, 1e-6));
}
