//
// Created by Andy on 12/26/2025.
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
using SecChem::Element;

TEST_CASE("BasisSet should compare by value")
{
	BasisSet bss0;
	BasisSet bss1;

	REQUIRE(bss0 == bss1);
	REQUIRE(bss0.EqualsTo(bss1));
}

TEST_CASE("SharedBasisSet should compare by address unless using named comparisons")
{
	SharedBasisSet bss0;
	SharedBasisSet bss1;

	REQUIRE(bss0 != bss1);
	REQUIRE(bss0.EqualsTo(bss1));
}

TEST_CASE("BasisSet should do deep copy")
{
	BasisSet bss0;
	BasisSet bss1 = bss0;

	REQUIRE(bss0 == bss1);
	REQUIRE(bss0.EqualsTo(bss1));

	bss0.AddEntryFor(Element::Br);
	REQUIRE(bss0 != bss1);
	REQUIRE(bss0.NotEqualsTo(bss1));
}

TEST_CASE("SharedBasisSet should do shallow copy unless .Clone() was called")
{
	SharedBasisSet bss0;
	SharedBasisSet bss1 = bss0;

	REQUIRE(bss0 == bss1);
	REQUIRE(bss0.EqualsTo(bss1));

	bss0.AddEntryFor(Element::Br);
	REQUIRE(bss0 == bss1);
	REQUIRE(bss0.EqualsTo(bss1));

	SharedBasisSet bss2 = bss0.Clone();
	REQUIRE(bss0 != bss2);
	REQUIRE(bss0.EqualsTo(bss2));

	bss0.AddEntryFor(Element::Cu);
	REQUIRE(bss0 != bss2);
	REQUIRE(bss0.NotEqualsTo(bss2));
}

