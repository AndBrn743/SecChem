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

TEST_CASE("SemiLocalEcp constructs correctly from column vectors", "[SemiLocalEcp]")
{
	Eigen::VectorXd coefficients(3);
	coefficients << 1.0, -0.5, 0.25;

	Eigen::VectorXd rExponents(3);
	rExponents << 0.0, 1.0, 2.0;

	Eigen::VectorXd gaussianExponents(3);
	gaussianExponents << 3.0, 2.0, 1.0;

	SemiLocalEcpTerm ecp(coefficients, rExponents, gaussianExponents);

	REQUIRE(ecp.Coefficients().size() == 3);
	REQUIRE(ecp.RExponents().size() == 3);
	REQUIRE(ecp.GaussianExponents().size() == 3);

	CHECK_THAT(ecp.Coefficient(0), WithinRel(1.0, 1e-14));
	CHECK_THAT(ecp.Coefficient(1), WithinRel(-0.5, 1e-14));
	CHECK_THAT(ecp.Coefficient(2), WithinRel(0.25, 1e-14));

	CHECK_THAT(ecp.RExponent(2), WithinRel(2.0, 1e-14));
	CHECK_THAT(ecp.GaussianExponent(0), WithinRel(3.0, 1e-14));
}

TEST_CASE("SemiLocalEcp accepts row vectors", "[SemiLocalEcp]")
{
	Eigen::RowVectorXd coefficients(3);
	coefficients << 0.1, 0.2, 0.3;

	Eigen::RowVectorXd rExponents(3);
	rExponents << 1.0, 2.0, 3.0;

	Eigen::RowVectorXd gaussianExponents(3);
	gaussianExponents << 4.0, 5.0, 6.0;

	SemiLocalEcpTerm ecp(coefficients, rExponents, gaussianExponents);

	CHECK_THAT(ecp.Coefficient(1), WithinRel(0.2, 1e-14));
	CHECK_THAT(ecp.RExponent(2), WithinRel(3.0, 1e-14));
	CHECK_THAT(ecp.GaussianExponent(0), WithinRel(4.0, 1e-14));
}

TEST_CASE("SemiLocalEcp accepts mixed row and column inputs", "[SemiLocalEcp]")
{
	Eigen::VectorXd coefficients(2);
	coefficients << -1.0, 1.0;

	Eigen::RowVectorXd rExponents(2);
	rExponents << 0.5, 1.5;

	Eigen::VectorXd gaussianExponents(2);
	gaussianExponents << 2.0, 3.0;

	SemiLocalEcpTerm ecp(coefficients, rExponents, gaussianExponents);

	CHECK_THAT(ecp.Coefficient(0), WithinRel(-1.0, 1e-14));
	CHECK_THAT(ecp.RExponent(1), WithinRel(1.5, 1e-14));
	CHECK_THAT(ecp.GaussianExponent(1), WithinRel(3.0, 1e-14));
}

TEST_CASE("SemiLocalEcp supports single-term ECP", "[SemiLocalEcp][edge]")
{
	Eigen::VectorXd coefficients(1);
	coefficients << -2.0;

	Eigen::VectorXd rExponents(1);
	rExponents << 1.0;

	Eigen::VectorXd gaussianExponents(1);
	gaussianExponents << 0.75;

	SemiLocalEcpTerm ecp(coefficients, rExponents, gaussianExponents);

	REQUIRE(ecp.Coefficients().size() == 1);
	CHECK_THAT(ecp.Coefficient(0), WithinRel(-2.0, 1e-14));
	CHECK_THAT(ecp.RExponent(0), WithinRel(1.0, 1e-14));
	CHECK_THAT(ecp.GaussianExponent(0), WithinRel(0.75, 1e-14));
}

TEST_CASE("SemiLocalEcp column accessors reference internal storage", "[SemiLocalEcp]")
{
	Eigen::VectorXd coefficients(3);
	coefficients << 1.0, 2.0, 3.0;

	Eigen::VectorXd rExponents(3);
	rExponents << 4.0, 5.0, 6.0;

	Eigen::VectorXd gaussianExponents(3);
	gaussianExponents << 7.0, 8.0, 9.0;

	SemiLocalEcpTerm ecp(coefficients, rExponents, gaussianExponents);

	const auto coeffs = ecp.Coefficients();
	const auto rExps = ecp.RExponents();
	const auto gExps = ecp.GaussianExponents();

	REQUIRE(coeffs.data() == &ecp.Coefficient(0));
	REQUIRE(rExps.data() == &ecp.RExponent(0));
	REQUIRE(gExps.data() == &ecp.GaussianExponent(0));
}

TEST_CASE("SemiLocalEcp handles large and negative values", "[SemiLocalEcp][numeric]")
{
	Eigen::VectorXd coefficients(3);
	coefficients << -1e3, 2e3, -3e3;

	Eigen::VectorXd rExponents(3);
	rExponents << 0.0, 5.0, 10.0;

	Eigen::VectorXd gaussianExponents(3);
	gaussianExponents << 1e-6, 1e-3, 1.0;

	SemiLocalEcpTerm ecp(coefficients, rExponents, gaussianExponents);

	CHECK_THAT(ecp.Coefficient(2), WithinRel(-3e3, 1e-14));
	CHECK_THAT(ecp.RExponent(1), WithinRel(5.0, 1e-14));
	CHECK_THAT(ecp.GaussianExponent(0), WithinRel(1e-6, 1e-14));
}

TEST_CASE("SemiLocalEcp accepts empty coefficient set", "[SemiLocalEcp]")
{
	Eigen::VectorXd coefficients(0);
	Eigen::VectorXd rExponents(0);
	Eigen::VectorXd gaussianExponents(0);

	SemiLocalEcpTerm ecp(coefficients, rExponents, gaussianExponents);

	REQUIRE(ecp.IsEmpty());
}

TEST_CASE("null SemiLocalEcp equal to each other", "[SemiLocalEcp]")
{
	Eigen::VectorXd coefficients(0);
	Eigen::VectorXd rExponents(0);
	Eigen::VectorXd gaussianExponents(0);

	SemiLocalEcpTerm ecp(coefficients, rExponents, gaussianExponents);
	SemiLocalEcpTerm ecp2{};

	REQUIRE(ecp.IsEmpty());
	REQUIRE(ecp2.IsEmpty());
	REQUIRE(ecp == ecp2);
}

TEST_CASE("SemiLocalEcp rejects size mismatch", "[SemiLocalEcp][invalid]")
{
	Eigen::VectorXd coefficients(3);
	coefficients << 1.0, 2.0, 3.0;

	Eigen::VectorXd rExponents(2);
	rExponents << 0.5, 1.5;

	Eigen::VectorXd gaussianExponents(3);
	gaussianExponents << 1.0, 2.0, 3.0;

	REQUIRE_THROWS_MATCHES(SemiLocalEcpTerm(coefficients, rExponents, gaussianExponents),
	                       std::invalid_argument,
	                       MessageMatches(ContainsSubstring("size")));
}

TEST_CASE("SemiLocalEcp rejects non-vector coefficient input", "[SemiLocalEcp][invalid]")
{
	Eigen::MatrixXd coefficients(2, 2);
	coefficients << 1.0, 2.0, 3.0, 4.0;

	Eigen::VectorXd rExponents(4);
	rExponents << 0.1, 0.2, 0.3, 0.4;

	Eigen::VectorXd gaussianExponents(4);
	gaussianExponents << 1.0, 2.0, 3.0, 4.0;

	REQUIRE_THROWS_MATCHES(SemiLocalEcpTerm(coefficients, rExponents, gaussianExponents),
	                       std::invalid_argument,
	                       MessageMatches(ContainsSubstring("vector")));
}

TEST_CASE("SemiLocalEcp rejects non-vector r-exponent input", "[SemiLocalEcp][invalid]")
{
	Eigen::VectorXd coefficients(2);
	coefficients << 1.0, -1.0;

	Eigen::MatrixXd rExponents(2, 2);
	rExponents << 0.5, 1.5, 2.5, 3.5;

	Eigen::VectorXd gaussianExponents(4);
	gaussianExponents << 1.0, 2.0, 3.0, 4.0;

	REQUIRE_THROWS_MATCHES(SemiLocalEcpTerm(coefficients, rExponents, gaussianExponents),
	                       std::invalid_argument,
	                       MessageMatches(ContainsSubstring("size mismatch")));
}

TEST_CASE("SemiLocalEcp rejects non-vector gaussian exponent input", "[SemiLocalEcp][invalid]")
{
	Eigen::VectorXd coefficients(2);
	coefficients << 0.5, 0.25;

	Eigen::VectorXd rExponents(2);
	rExponents << 1.0, 2.0;

	Eigen::MatrixXd gaussianExponents(2, 2);
	gaussianExponents << 0.1, 0.2, 0.3, 0.4;

	REQUIRE_THROWS_MATCHES(SemiLocalEcpTerm(coefficients, rExponents, gaussianExponents),
	                       std::invalid_argument,
	                       MessageMatches(ContainsSubstring("size mismatch")));
}

TEST_CASE("SemiLocalEcp rejects multiple simultaneous input errors", "[SemiLocalEcp][invalid]")
{
	Eigen::MatrixXd coefficients(2, 2);
	coefficients << 1.0, 2.0, 3.0, 4.0;

	Eigen::VectorXd rExponents(3);
	rExponents << 0.1, 0.2, 0.3;

	Eigen::RowVectorXd gaussianExponents(1);
	gaussianExponents << 0.5;

	REQUIRE_THROWS_AS(SemiLocalEcpTerm(coefficients, rExponents, gaussianExponents), std::invalid_argument);
}

TEST_CASE("SemiLocalEcp ValidateInput accepts valid mixed input", "[SemiLocalEcp][valid]")
{
	Eigen::RowVectorXd coefficients(3);
	coefficients << 1.0, -0.5, 0.25;

	Eigen::VectorXd rExponents(3);
	rExponents << 0.0, 1.0, 2.0;

	Eigen::RowVectorXd gaussianExponents(3);
	gaussianExponents << 3.0, 2.0, 1.0;

	REQUIRE_NOTHROW(SemiLocalEcpTerm(coefficients, rExponents, gaussianExponents));
}

TEST_CASE("SemiLocalEcp's equality comparison should work", "[SemiLocalEcp]")
{
	Eigen::RowVectorXd coefficients(3);
	coefficients << 1.0, -0.5, 0.25;

	Eigen::VectorXd rExponents(3);
	rExponents << 0.0, 1.0, 2.0;

	Eigen::RowVectorXd gaussianExponents(3);
	gaussianExponents << 3.0, 2.0, 1.0;

	SemiLocalEcpTerm ecp0(coefficients, rExponents, gaussianExponents);

	{
		SemiLocalEcpTerm ecp1 = ecp0;
		REQUIRE(ecp1 == ecp0);
		REQUIRE_FALSE(ecp1 != ecp0);
	}

	{
		Eigen::RowVectorXd gaussianExponentsToo = gaussianExponents;
		gaussianExponentsToo[0] += 1e-6;

		SemiLocalEcpTerm ecp2(coefficients, rExponents, gaussianExponentsToo);
		REQUIRE_FALSE(ecp2 == ecp0);
		REQUIRE(ecp2 != ecp0);
		REQUIRE_FALSE(ecp2.EqualsTo(ecp0, 1e-9));
		REQUIRE(ecp2.EqualsTo(ecp0, 1e-5));
	}

	{
		SemiLocalEcpTerm ecp3(Eigen::VectorXd::Random(4), Eigen::VectorXd::Random(4), Eigen::VectorXd::Random(4));
		REQUIRE_FALSE(ecp3 == ecp0);
		REQUIRE(ecp3 != ecp0);
		REQUIRE_FALSE(ecp3.EqualsTo(ecp0));
		REQUIRE(ecp3.NotEqualsTo(ecp0));
	}
}

static const auto ecpExampleA = []
{
	Eigen::VectorXd coefficients(3);
	coefficients << 1.0, -0.5, 0.25;

	Eigen::VectorXd rExponents(3);
	rExponents << 0.0, 1.0, 2.0;

	Eigen::VectorXd gaussianExponents(3);
	gaussianExponents << 3.0, 2.0, 1.0;

	SemiLocalEcpTerm ecp(coefficients, rExponents, gaussianExponents);
	return ecp;
}();

static const auto ecpExampleB = []
{
	Eigen::VectorXd coefficients(2);
	coefficients << 10.0, -10.5;

	Eigen::VectorXd rExponents(2);
	rExponents << 1.9, 1.8;

	Eigen::VectorXd gaussianExponents(2);
	gaussianExponents << 3.2, 2.2;

	SemiLocalEcpTerm ecp(coefficients, rExponents, gaussianExponents);
	return ecp;
}();

TEST_CASE("SemiLocalEcp::Concat should work", "[SemiLocalEcp]")
{
	const std::array ecps = {ecpExampleA, ecpExampleB};

	SECTION("Concat iterable")
	{
		auto ecp = SemiLocalEcpTerm::Concat(ecps.begin(), ecps.end());

		REQUIRE(ecp.Coefficients().size() == 5);
		REQUIRE(ecp.Coefficients().head(3) == ecpExampleA.Coefficients());
		REQUIRE(ecp.Coefficients().tail(2) == ecpExampleB.Coefficients());
		REQUIRE(ecp.RExponents().head(3) == ecpExampleA.RExponents());
		REQUIRE(ecp.RExponents().tail(2) == ecpExampleB.RExponents());
		REQUIRE(ecp.GaussianExponents().head(3) == ecpExampleA.GaussianExponents());
		REQUIRE(ecp.GaussianExponents().tail(2) == ecpExampleB.GaussianExponents());
	}

	SECTION("Concat items")
	{
		auto ecp = SemiLocalEcpTerm::Concat(ecps[0], ecps[1]);

		REQUIRE(ecp.Coefficients().size() == 5);
		REQUIRE(ecp.Coefficients().head(3) == ecpExampleA.Coefficients());
		REQUIRE(ecp.Coefficients().tail(2) == ecpExampleB.Coefficients());
		REQUIRE(ecp.RExponents().head(3) == ecpExampleA.RExponents());
		REQUIRE(ecp.RExponents().tail(2) == ecpExampleB.RExponents());
		REQUIRE(ecp.GaussianExponents().head(3) == ecpExampleA.GaussianExponents());
		REQUIRE(ecp.GaussianExponents().tail(2) == ecpExampleB.GaussianExponents());
	}
}

TEST_CASE("SemiLocalEcp::Concat should throw with empty input range", "[SemiLocalEcp]")
{
	std::vector<SemiLocalEcpTerm> empty{};
	REQUIRE_THROWS(SemiLocalEcpTerm::Concat(empty.begin(), empty.end()));
}

TEST_CASE("SemiLocalEcp::Concat should work with single item range", "[SemiLocalEcp]")
{
	const std::array ecps = {ecpExampleA};

	SECTION("Concat iterable")
	{
		auto ecp = SemiLocalEcpTerm::Concat(ecps.begin(), ecps.end());
		REQUIRE(ecp == ecpExampleA);
	}

	SECTION("Concat item")
	{
		auto ecp = SemiLocalEcpTerm::Concat(ecps[0]);
		REQUIRE(ecp == ecpExampleA);
	}
}
