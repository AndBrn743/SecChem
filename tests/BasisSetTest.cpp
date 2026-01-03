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

#include <SecChem/BasisSet/Gaussian.hpp>
using namespace SecChem::BasisSet::Gaussian;
using SecChem::AzimuthalQuantumNumber;
using SecChem::Element;

TEST_CASE("BasisSet should compare by value", "[BasisSet]")
{
	BasisSet bss0;
	BasisSet bss1;

	REQUIRE_FALSE(bss0.Has(Element::Neutron));
	REQUIRE_FALSE(bss1.Has(Element::Neutron));
	REQUIRE_FALSE(bss0.Has(Element::H));
	REQUIRE_FALSE(bss1.Has(Element::H));
	REQUIRE(bss0 == bss1);
	REQUIRE(bss0.EqualsTo(bss1));
	REQUIRE(bss0.EqualsTo(bss0));
	REQUIRE(bss1.EqualsTo(bss1));
}

TEST_CASE("SharedBasisSet should compare by address unless using named comparisons", "[BasisSet]")
{
	SharedBasisSet bss0;
	SharedBasisSet bss1;

	REQUIRE(bss0 != bss1);
	REQUIRE(bss0.EqualsTo(bss1));
}

TEST_CASE("BasisSet should do deep copy", "[BasisSet]")
{
	BasisSet bss0;
	BasisSet bss1 = bss0;

	REQUIRE(bss0 == bss1);
	REQUIRE(bss0.EqualsTo(bss1));

	bss0.AddEntryFor(Element::Br);
	REQUIRE(bss0 != bss1);
	REQUIRE(bss0.NotEqualsTo(bss1));
}

TEST_CASE("SharedBasisSet should do shallow copy unless .Clone() was called", "[BasisSet]")
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

TEST_CASE("AddEntryFor / OverwriteEntryOf behavior", "[BasisSet]")
{
	BasisSet bs;

	auto& h = bs.AddEntryFor(Element::H);
	REQUIRE(h.empty());

	REQUIRE_THROWS_AS(bs.AddEntryFor(Element::H), std::logic_error);

	auto& h2 = bs.OverwriteEntryOf(Element::H);
	REQUIRE(h2.empty());

	auto& h3 = bs.AddOrOverwriteEntryOf(Element::H);
	REQUIRE(h3.empty());
}

TEST_CASE("operator== is representation-based", "[BasisSet]")
{
	BasisSet a;
	BasisSet b;

	REQUIRE(a == b);  // distinct storage

	a.AddEntryFor(Element::H)
	        .emplace_back(AzimuthalQuantumNumber{0},
	                      ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 1.0),
	                                                 Eigen::MatrixXd::Constant(1, 1, 1.0)});

	SharedBasisSet refA{a};

	REQUIRE(refA == refA);
	STATIC_REQUIRE_FALSE(refA == a);  // different semantics → always false
}

TEST_CASE("EqualsTo performs deep comparison", "[BasisSet]")
{
	BasisSet a;
	a.AddEntryFor(Element::H)
	        .emplace_back(
	                AzimuthalQuantumNumber{0},
	                ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(2, 1.0), Eigen::MatrixXd::Identity(2, 1)});

	SharedBasisSet refA{a};
	BasisSet b{refA};

	REQUIRE(a.EqualsTo(refA));
	REQUIRE(refA.EqualsTo(b));
	REQUIRE(a.EqualsTo(b));
}

TEST_CASE("Clone creates deep copy for Reference semantics", "[BasisSet]")
{
	BasisSet value;
	value.AddEntryFor(Element::H)
	        .emplace_back(AzimuthalQuantumNumber::F,
	                      ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 2.0),
	                                                 Eigen::MatrixXd::Constant(1, 1, 3.0)});

	SharedBasisSet ref{value};
	auto cloned = ref.Clone();

	REQUIRE(ref != cloned);         // different shared_ptr
	REQUIRE(ref.EqualsTo(cloned));  // same value
}

TEST_CASE("StandardizeRepresentation sorts AngularMomentumBlocks", "[BasisSet]")
{
	BasisSet bs;
	REQUIRE(bs.IsInStandardRepresentation());

	auto& c = bs.AddEntryFor(Element::Carbon);
	REQUIRE_FALSE(bs.IsInStandardRepresentation());

	c.emplace_back(AzimuthalQuantumNumber{2},
	               ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 2.0), Eigen::MatrixXd::Ones(1, 1)});
	c.emplace_back(AzimuthalQuantumNumber{0},
	               ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 1.0), Eigen::MatrixXd::Ones(1, 1)});
	c.emplace_back(AzimuthalQuantumNumber{1},
	               ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 1.5), Eigen::MatrixXd::Ones(1, 1)});

	REQUIRE_FALSE(bs.IsInStandardRepresentation());

	bs.StandardizeRepresentation();

	REQUIRE(bs.IsInStandardRepresentation());
	REQUIRE(c[0].AngularMomentum() == AzimuthalQuantumNumber{0});
	REQUIRE(c[1].AngularMomentum() == AzimuthalQuantumNumber{1});
	REQUIRE(c[2].AngularMomentum() == AzimuthalQuantumNumber{2});
}

TEST_CASE("StandardizeRepresentation concatenates same-l blocks", "[BasisSet]")
{
	BasisSet bs;
	auto& h = bs.AddEntryFor(Element::H);

	h.emplace_back(AzimuthalQuantumNumber::S,
	               ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 1.0), Eigen::MatrixXd::Constant(1, 1, 1.0)});

	h.emplace_back(AzimuthalQuantumNumber::S,
	               ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 2.0), Eigen::MatrixXd::Constant(1, 1, 2.0)});

	REQUIRE_FALSE(bs.IsInStandardRepresentation());

	bs.StandardizeRepresentation();

	REQUIRE(bs.IsInStandardRepresentation());
	REQUIRE(h.size() == 1);
}

TEST_CASE("ToStandardizedRepresentation does not mutate original", "[BasisSet]")
{
	BasisSet bs;
	bs.AddEntryFor(Element::Neutron);
	auto& h = bs.AddEntryFor(Element::H);

	h.emplace_back(AzimuthalQuantumNumber::S,
	               ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 1.0), Eigen::MatrixXd::Constant(1, 1, 1.1)});

	h.emplace_back(AzimuthalQuantumNumber::S,
	               ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 2.0), Eigen::MatrixXd::Constant(1, 1, 2.1)});

	h.emplace_back(AzimuthalQuantumNumber::P,
	               ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 3.0), Eigen::MatrixXd::Constant(1, 1, 3.1)});

	REQUIRE_FALSE(bs.IsInStandardRepresentation());

	auto copy = bs;
	auto standardized = bs.ToStandardizedRepresentation();

	REQUIRE(bs == copy);
	REQUIRE_FALSE(bs.IsInStandardRepresentation());
	REQUIRE(standardized.IsInStandardRepresentation());
}


TEST_CASE("Misc .EqualsTo(...)", "[BasisSet]")
{
	SharedBasisSet bss0;
	SharedBasisSet bss1;

	REQUIRE(bss0.EqualsTo(bss1));

	SECTION("Value diff")
	{
		bss0.AddEntryFor(Element::H)
		        .emplace_back(AzimuthalQuantumNumber::P,
		                      ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 3.0),
		                                                 Eigen::MatrixXd::Constant(1, 1, 3.1)});
		bss1.AddEntryFor(Element::H)
		        .emplace_back(AzimuthalQuantumNumber::P,
		                      ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 3.0),
		                                                 Eigen::MatrixXd::Constant(1, 1, 3.2)});
		REQUIRE_FALSE(bss0.EqualsTo(bss1));
	}

	SECTION("AzimuthalQuantumNumber range diff")
	{
		auto& h0 = bss0.AddEntryFor(Element::H);
		h0.emplace_back(
		        AzimuthalQuantumNumber::S,
		        ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 3.0), Eigen::MatrixXd::Constant(1, 1, 3.1)});
		h0.emplace_back(
		        AzimuthalQuantumNumber::P,
		        ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(2, 3.0), Eigen::MatrixXd::Constant(2, 1, 3.1)});
		bss1.AddEntryFor(Element::H)
		        .emplace_back(AzimuthalQuantumNumber::S,
		                      ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 3.0),
		                                                 Eigen::MatrixXd::Constant(1, 1, 3.1)});
		REQUIRE_FALSE(bss0.EqualsTo(bss1));
	}

	SECTION("Element diff")
	{
		bss0.AddEntryFor(Element::H);
		bss1.AddEntryFor(Element::C);
		REQUIRE_FALSE(bss0.EqualsTo(bss1));
	}
}

TEST_CASE("Misc IsInStandardRepresentation and StandardizeRepresentation tests", "[BasisSet]")
{
	BasisSet bs;
	REQUIRE(bs.IsInStandardRepresentation());

	bs.AddEntryFor(Element::Carbon);
	REQUIRE_FALSE(bs.IsInStandardRepresentation());
	bs.StandardizeRepresentation();
	REQUIRE(bs.IsInStandardRepresentation());
	REQUIRE_FALSE(bs.Has(Element::Carbon));
}
