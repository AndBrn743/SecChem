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
	REQUIRE(h.AzimuthalShells.empty());
	REQUIRE(h.EcpElectronCount == 0);

	REQUIRE_THROWS_AS(bs.AddEntryFor(Element::H), std::logic_error);

	auto& h2 = bs.OverwriteEntryOf(Element::H);
	REQUIRE(h2.AzimuthalShells.empty());
	REQUIRE(h2.EcpElectronCount == 0);

	auto& h3 = bs.AddOrOverwriteEntryOf(Element::H);
	REQUIRE(h3.AzimuthalShells.empty());
	REQUIRE(h3.EcpElectronCount == 0);
}

TEST_CASE("operator== is representation-based", "[BasisSet]")
{
	BasisSet a;
	BasisSet b;

	REQUIRE(a == b);  // distinct storage

	a.AddEntryFor(Element::H).AzimuthalShells
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
	a.AddEntryFor(Element::H).AzimuthalShells
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
	value.AddEntryFor(Element::H).AzimuthalShells
	        .emplace_back(AzimuthalQuantumNumber::F,
	                      ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 2.0),
	                                                 Eigen::MatrixXd::Constant(1, 1, 3.0)});

	SharedBasisSet ref{value};
	auto cloned = ref.Clone();

	REQUIRE(ref != cloned);         // different shared_ptr
	REQUIRE(ref.EqualsTo(cloned));  // same value
}

TEST_CASE("StandardizeRepresentation sorts AzimuthalShells", "[BasisSet]")
{
	BasisSet bs;
	REQUIRE(bs.IsInStandardRepresentation());

	auto& c = bs.AddEntryFor(Element::Carbon).AzimuthalShells;
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
	auto& h = bs.AddEntryFor(Element::H).AzimuthalShells;

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
	auto& h = bs.AddEntryFor(Element::H).AzimuthalShells;

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
		bss0.AddEntryFor(Element::H).AzimuthalShells
		        .emplace_back(AzimuthalQuantumNumber::P,
		                      ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 3.0),
		                                                 Eigen::MatrixXd::Constant(1, 1, 3.1)});
		bss1.AddEntryFor(Element::H).AzimuthalShells
		        .emplace_back(AzimuthalQuantumNumber::P,
		                      ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 3.0),
		                                                 Eigen::MatrixXd::Constant(1, 1, 3.2)});
		REQUIRE_FALSE(bss0.EqualsTo(bss1));
	}

	SECTION("AzimuthalQuantumNumber range diff")
	{
		auto& h0 = bss0.AddEntryFor(Element::H).AzimuthalShells;
		h0.emplace_back(
		        AzimuthalQuantumNumber::S,
		        ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 3.0), Eigen::MatrixXd::Constant(1, 1, 3.1)});
		h0.emplace_back(
		        AzimuthalQuantumNumber::P,
		        ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(2, 3.0), Eigen::MatrixXd::Constant(2, 1, 3.1)});
		bss1.AddEntryFor(Element::H).AzimuthalShells
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

// =============================================================================
// BasisSetLibrary Builder Tests
// =============================================================================

TEST_CASE("Builder<BasisSetLibrary> creates empty library", "[BasisSetLibrary][Builder]")
{
	SecChem::Builder<BasisSetLibrary> builder;
	auto library = builder.Build();

	REQUIRE_FALSE(library.Has("sto-3g"));
	REQUIRE_FALSE(library.Has("6-31g"));
}

TEST_CASE("Builder<BasisSetLibrary> can add single basis set", "[BasisSetLibrary][Builder]")
{
	SecChem::Builder<BasisSetLibrary> builder;

	BasisSet bs;
	bs.AddEntryFor(Element::H).AzimuthalShells.emplace_back(
	        AzimuthalQuantumNumber::S,
	        ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 1.0), Eigen::MatrixXd::Constant(1, 1, 1.0)});

	builder.AddBasisSet("minimal-h", std::move(bs));
	auto library = builder.Build();

	REQUIRE(library.Has("minimal-h"));
	REQUIRE_FALSE(library.Has("6-31g"));
	REQUIRE(library["minimal-h"].Has(Element::H));
}

TEST_CASE("Builder<BasisSetLibrary> can add multiple basis sets", "[BasisSetLibrary][Builder]")
{
	SecChem::Builder<BasisSetLibrary> builder;

	BasisSet sto3g;
	sto3g.AddEntryFor(Element::H).AzimuthalShells.emplace_back(
	        AzimuthalQuantumNumber::S,
	        ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 1.0), Eigen::MatrixXd::Constant(1, 1, 1.0)});
	sto3g.AddEntryFor(Element::C).AzimuthalShells.emplace_back(
	        AzimuthalQuantumNumber::S,
	        ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 2.0), Eigen::MatrixXd::Constant(1, 1, 2.0)});

	BasisSet cc_pvdz;
	cc_pvdz.AddEntryFor(Element::O).AzimuthalShells.emplace_back(
	        AzimuthalQuantumNumber::P,
	        ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 3.0), Eigen::MatrixXd::Constant(1, 1, 3.0)});

	builder.AddBasisSet("sto-3g", std::move(sto3g));
	builder.AddBasisSet("cc-pvdz", std::move(cc_pvdz));

	auto library = builder.Build();

	REQUIRE(library.Has("sto-3g"));
	REQUIRE(library.Has("cc-pvdz"));
	REQUIRE(library["sto-3g"].Has(Element::H));
	REQUIRE(library["sto-3g"].Has(Element::C));
	REQUIRE(library["cc-pvdz"].Has(Element::O));
}

TEST_CASE("Builder<BasisSetLibrary> rejects duplicate basis set names", "[BasisSetLibrary][Builder]")
{
	SecChem::Builder<BasisSetLibrary> builder;

	BasisSet bs1;
	bs1.AddEntryFor(Element::H);

	BasisSet bs2;
	bs2.AddEntryFor(Element::C);

	builder.AddBasisSet("basis1", std::move(bs1));

	REQUIRE_THROWS_MATCHES(
	        builder.AddBasisSet("basis1", std::move(bs2)),
	        std::logic_error,
	        MessageMatches(ContainsSubstring("same name")));
}

TEST_CASE("Builder<BasisSetLibrary> throws when already consumed", "[BasisSetLibrary][Builder]")
{
	SecChem::Builder<BasisSetLibrary> builder;

	BasisSet bs;
	bs.AddEntryFor(Element::H);

	builder.AddBasisSet("test", std::move(bs));
	auto library = builder.Build();

	SECTION("Cannot Build twice")
	{
		REQUIRE_THROWS_MATCHES(builder.Build(), std::logic_error, MessageMatches(ContainsSubstring("already consumed")));
	}

	SECTION("Cannot AddBasisSet after Build")
	{
		BasisSet bs2;
		bs2.AddEntryFor(Element::C);

		REQUIRE_THROWS_MATCHES(
		        builder.AddBasisSet("test2", std::move(bs2)),
		        std::logic_error,
		        MessageMatches(ContainsSubstring("already consumed")));
	}
}

TEST_CASE("Builder<BasisSetLibrary> supports chaining", "[BasisSetLibrary][Builder]")
{
	SecChem::Builder<BasisSetLibrary> builder;

	BasisSet bs1;
	bs1.AddEntryFor(Element::H);

	BasisSet bs2;
	bs2.AddEntryFor(Element::C);

	BasisSet bs3;
	bs3.AddEntryFor(Element::O);

	auto library = builder.AddBasisSet("bs1", std::move(bs1))
	                      .AddBasisSet("bs2", std::move(bs2))
	                      .AddBasisSet("bs3", std::move(bs3))
	                      .Build();

	REQUIRE(library.Has("bs1"));
	REQUIRE(library.Has("bs2"));
	REQUIRE(library.Has("bs3"));
}

TEST_CASE("Builder<SharedBasisSetLibrary> creates shared library", "[BasisSetLibrary][Builder]")
{
	SecChem::Builder<SharedBasisSetLibrary> builder;

	SharedBasisSet bs;
	bs.AddEntryFor(Element::H).AzimuthalShells.emplace_back(
	        AzimuthalQuantumNumber::S,
	        ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 1.0), Eigen::MatrixXd::Constant(1, 1, 1.0)});

	builder.AddBasisSet("test", std::move(bs));
	auto library = builder.Build();

	REQUIRE(library.Has("test"));
	REQUIRE(library["test"].Has(Element::H));
}

TEST_CASE("Builder<SharedBasisSetLibrary> throws when already consumed", "[BasisSetLibrary][Builder]")
{
	SecChem::Builder<SharedBasisSetLibrary> builder;

	SharedBasisSet bs;
	bs.AddEntryFor(Element::H);

	builder.AddBasisSet("test", std::move(bs));
	builder.Build();

	REQUIRE_THROWS_MATCHES(builder.Build(), std::logic_error, MessageMatches(ContainsSubstring("already consumed")));
}

TEST_CASE("BasisSetLibrary operator[] throws for missing name", "[BasisSetLibrary]")
{
	SecChem::Builder<BasisSetLibrary> builder;
	auto library = builder.Build();

	REQUIRE_THROWS_AS(library["nonexistent"], std::out_of_range);
}

TEST_CASE("Builder<BasisSetLibrary> with empty basis set", "[BasisSetLibrary][Builder]")
{
	SecChem::Builder<BasisSetLibrary> builder;

	BasisSet emptyBs;
	builder.AddBasisSet("empty", std::move(emptyBs));

	auto library = builder.Build();

	REQUIRE(library.Has("empty"));
	REQUIRE_FALSE(library["empty"].Has(Element::H));
	REQUIRE_FALSE(library["empty"].Has(Element::C));
}

// =============================================================================
// CreatePrincipalQuantumNumberOffsetTable Tests
// =============================================================================
//  S   P   D   F   G
//  2,
//  4, 10,
// 12, 18, 28,
// 30, 36, 46, 60,
// 62, 68, 78, 92, 110.

TEST_CASE("CreatePrincipalQuantumNumberOffsetTable returns all zeros when EcpElectronCount is 0", "[ElementaryBasisSet]")
{
	ElementaryBasisSet ebs;
	ebs.AzimuthalShells.emplace_back(
	        AzimuthalQuantumNumber::G,
	        ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 1.0), Eigen::MatrixXd::Constant(1, 1, 1.0)});

	auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();

	REQUIRE(offsets.size() == 5);  // l=0,1,2,3,4
	REQUIRE(offsets[0] == 0);
	REQUIRE(offsets[1] == 0);
	REQUIRE(offsets[2] == 0);
	REQUIRE(offsets[3] == 0);
	REQUIRE(offsets[4] == 0);
}

TEST_CASE("CreatePrincipalQuantumNumberOffsetTable for l=0 (s orbitals)", "[ElementaryBasisSet]")
{
	ElementaryBasisSet ebs;
	ebs.AzimuthalShells.emplace_back(
	        AzimuthalQuantumNumber::S,
	        ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 1.0), Eigen::MatrixXd::Constant(1, 1, 1.0)});

	SECTION("EcpElectronCount < 2 gives offset 0")
	{
		ebs.EcpElectronCount = 0;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[0] == 0);
	}

	SECTION("EcpElectronCount = 2 gives offset 1")
	{
		ebs.EcpElectronCount = 2;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[0] == 1);
	}

	SECTION("EcpElectronCount = 4 gives offset 2")
	{
		ebs.EcpElectronCount = 4;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[0] == 2);
	}

	SECTION("EcpElectronCount = 12 gives offset 3")
	{
		ebs.EcpElectronCount = 12;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[0] == 3);
	}

	SECTION("EcpElectronCount = 30 gives offset 4")
	{
		ebs.EcpElectronCount = 30;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[0] == 4);
	}

	SECTION("EcpElectronCount = 62 gives offset 5")
	{
		ebs.EcpElectronCount = 62;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[0] == 5);
	}

	SECTION("EcpElectronCount between thresholds uses lower threshold")
	{
		ebs.EcpElectronCount = 20;  // between 12 and 30
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[0] == 3);  // should use 12 -> 3

		ebs.EcpElectronCount = 50;  // between 30 and 62
		offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[0] == 4);  // should use 30 -> 4
	}
}

TEST_CASE("CreatePrincipalQuantumNumberOffsetTable for l=1 (p orbitals)", "[ElementaryBasisSet]")
{
	ElementaryBasisSet ebs;
	ebs.AzimuthalShells.emplace_back(
	        AzimuthalQuantumNumber::P,
	        ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 1.0), Eigen::MatrixXd::Constant(1, 1, 1.0)});

	SECTION("EcpElectronCount < 10 gives offset 0")
	{
		ebs.EcpElectronCount = 4;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[1] == 0);
	}

	SECTION("EcpElectronCount = 10 gives offset 1")
	{
		ebs.EcpElectronCount = 10;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[1] == 1);
	}

	SECTION("EcpElectronCount = 18 gives offset 2")
	{
		ebs.EcpElectronCount = 18;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[1] == 2);
	}

	SECTION("EcpElectronCount = 36 gives offset 3")
	{
		ebs.EcpElectronCount = 36;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[1] == 3);
	}

	SECTION("EcpElectronCount = 68 gives offset 4")
	{
		ebs.EcpElectronCount = 68;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[1] == 4);
	}
}

TEST_CASE("CreatePrincipalQuantumNumberOffsetTable for l=2 (d orbitals)", "[ElementaryBasisSet]")
{
	ElementaryBasisSet ebs;
	ebs.AzimuthalShells.emplace_back(
	        AzimuthalQuantumNumber::D,
	        ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 1.0), Eigen::MatrixXd::Constant(1, 1, 1.0)});

	SECTION("EcpElectronCount < 28 gives offset 0")
	{
		ebs.EcpElectronCount = 18;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[2] == 0);
	}

	SECTION("EcpElectronCount = 28 gives offset 1")
	{
		ebs.EcpElectronCount = 28;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[2] == 1);
	}

	SECTION("EcpElectronCount = 60 gives offset 2")
	{
		ebs.EcpElectronCount = 60;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[2] == 2);
	}

	SECTION("EcpElectronCount = 78 gives offset 3")
	{
		ebs.EcpElectronCount = 78;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[2] == 3);
	}
}

TEST_CASE("CreatePrincipalQuantumNumberOffsetTable for l=3 (f orbitals)", "[ElementaryBasisSet]")
{
	ElementaryBasisSet ebs;
	ebs.AzimuthalShells.emplace_back(
	        AzimuthalQuantumNumber::F,
	        ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 1.0), Eigen::MatrixXd::Constant(1, 1, 1.0)});

	SECTION("EcpElectronCount < 60 gives offset 0")
	{
		ebs.EcpElectronCount = 50;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[3] == 0);
	}

	SECTION("EcpElectronCount = 60 gives offset 1")
	{
		ebs.EcpElectronCount = 60;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[3] == 1);
	}

	SECTION("EcpElectronCount = 92 gives offset 2")
	{
		ebs.EcpElectronCount = 92;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[3] == 2);
	}
}

TEST_CASE("CreatePrincipalQuantumNumberOffsetTable for l=4 (g orbitals)", "[ElementaryBasisSet]")
{
	ElementaryBasisSet ebs;
	ebs.AzimuthalShells.emplace_back(
	        AzimuthalQuantumNumber::G,
	        ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 1.0), Eigen::MatrixXd::Constant(1, 1, 1.0)});

	SECTION("EcpElectronCount < 110 gives offset 0")
	{
		ebs.EcpElectronCount = 100;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[4] == 0);
	}

	SECTION("EcpElectronCount = 110 gives offset 1")
	{
		ebs.EcpElectronCount = 110;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[4] == 1);
	}
}

TEST_CASE("CreatePrincipalQuantumNumberOffsetTable with multiple angular momenta", "[ElementaryBasisSet]")
{
	ElementaryBasisSet ebs;
	ebs.AzimuthalShells.emplace_back(
	        AzimuthalQuantumNumber::S,
	        ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 1.0), Eigen::MatrixXd::Constant(1, 1, 1.0)});
	ebs.AzimuthalShells.emplace_back(
	        AzimuthalQuantumNumber::P,
	        ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 1.0), Eigen::MatrixXd::Constant(1, 1, 1.0)});
	ebs.AzimuthalShells.emplace_back(
	        AzimuthalQuantumNumber::D,
	        ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 1.0), Eigen::MatrixXd::Constant(1, 1, 1.0)});

	ebs.EcpElectronCount = 46;  // common ECP core size (e.g., for Ag, Pd)

	auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();

	REQUIRE(offsets.size() == 3);
	REQUIRE(offsets[0] == 4);  // l=0: 46 >= 30, so offset is 4
	REQUIRE(offsets[1] == 3);  // l=1: 46 >= 36, so offset is 2
	REQUIRE(offsets[2] == 2);  // l=2: 46 >= 28, so offset is 1
}

TEST_CASE("CreatePrincipalQuantumNumberOffsetTable table-driven implementation correctness", "[ElementaryBasisSet]")
{
	// Test specific ECP electron counts that correspond to known elements
	ElementaryBasisSet ebs;
	ebs.AzimuthalShells.emplace_back(
	        AzimuthalQuantumNumber::G,
	        ContractedRadialOrbitalSet{Eigen::VectorXd::Constant(1, 1.0), Eigen::MatrixXd::Constant(1, 1, 1.0)});

	SECTION("EcpElectronCount = 10")
	{
		ebs.EcpElectronCount = 10;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[0] == 2);  // s: 10 >= 4
		REQUIRE(offsets[1] == 1);  // p: 10 >= 10
		REQUIRE(offsets[2] == 0);  // d: 10 < 28
		REQUIRE(offsets[3] == 0);  // f: 10 < 60
		REQUIRE(offsets[4] == 0);  // g: 10 < 110
	}

	SECTION("EcpElectronCount = 28")
	{
		ebs.EcpElectronCount = 28;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[0] == 3);  // s: 28 >= 12
		REQUIRE(offsets[1] == 2);  // p: 28 >= 18
		REQUIRE(offsets[2] == 1);  // d: 28 >= 28
		REQUIRE(offsets[3] == 0);  // f: 28 < 60
		REQUIRE(offsets[4] == 0);  // g: 28 < 110
	}

	SECTION("EcpElectronCount = 60")
	{
		ebs.EcpElectronCount = 60;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[0] == 4);  // s: 60 >= 30
		REQUIRE(offsets[1] == 3);  // p: 60 >= 36
		REQUIRE(offsets[2] == 2);  // d: 60 >= 60
		REQUIRE(offsets[3] == 1);  // f: 60 >= 60
		REQUIRE(offsets[4] == 0);  // g: 60 < 110
	}

	SECTION("EcpElectronCount = 78")
	{
		ebs.EcpElectronCount = 78;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[0] == 5);  // s: 78 >= 62
		REQUIRE(offsets[1] == 4);  // p: 78 >= 68
		REQUIRE(offsets[2] == 3);  // d: 78 >= 78
		REQUIRE(offsets[3] == 1);  // f: 78 >= 60
		REQUIRE(offsets[4] == 0);  // g: 78 < 110
	}

	SECTION("EcpElectronCount = 110 (max supported)")
	{
		ebs.EcpElectronCount = 110;
		auto offsets = ebs.CreatePrincipalQuantumNumberOffsetTable();
		REQUIRE(offsets[0] == 5);  // s: 110 >= 62
		REQUIRE(offsets[1] == 4);  // p: 110 >= 68
		REQUIRE(offsets[2] == 3);  // d: 110 >= 78
		REQUIRE(offsets[3] == 2);  // f: 110 >= 92
		REQUIRE(offsets[4] == 1);  // g: 110 >= 110
	}
}

// =============================================================================
// IsInStandardStorageOrder_OverloadSet for SemiLocalEcpChannel Tests
// =============================================================================

TEST_CASE("StandardizeRepresentation sorts SemiLocalEcpChannels by angular momentum", "[ElementaryBasisSet]")
{
	ElementaryBasisSet ebs;

	// Add channels in wrong order: P, S, D (should be S, P, D)
	Eigen::VectorXd coeffs(1);
	coeffs << 1.0;
	Eigen::VectorXd rExps(1);
	rExps << 0.0;
	Eigen::VectorXd gaussExps(1);
	gaussExps << 2.0;

	ebs.SemiLocalEcpChannels.emplace_back(AzimuthalQuantumNumber::P, coeffs, rExps, gaussExps);
	ebs.SemiLocalEcpChannels.emplace_back(AzimuthalQuantumNumber::S, coeffs, rExps, gaussExps);
	ebs.SemiLocalEcpChannels.emplace_back(AzimuthalQuantumNumber::D, coeffs, rExps, gaussExps);

	REQUIRE_FALSE(ebs.IsInStandardRepresentation());

	ebs.StandardizeRepresentation();

	REQUIRE(ebs.IsInStandardRepresentation());
	REQUIRE(ebs.SemiLocalEcpChannels[0].AngularMomentum() == AzimuthalQuantumNumber::S);
	REQUIRE(ebs.SemiLocalEcpChannels[1].AngularMomentum() == AzimuthalQuantumNumber::P);
	REQUIRE(ebs.SemiLocalEcpChannels[2].AngularMomentum() == AzimuthalQuantumNumber::D);
}

TEST_CASE("StandardizeRepresentation sorts SemiLocalEcpChannels by term count descending", "[ElementaryBasisSet]")
{
	ElementaryBasisSet ebs;

	// Two P channels: one with 3 terms, one with 2 terms
	// The 3-term channel should come first
	Eigen::VectorXd coeffs3(3);
	coeffs3 << 1.0, 0.5, 0.25;
	Eigen::VectorXd rExps3(3);
	rExps3 << 0.0, 1.0, 2.0;
	Eigen::VectorXd gaussExps3(3);
	gaussExps3 << 5.0, 4.0, 3.0;

	Eigen::VectorXd coeffs2(2);
	coeffs2 << 1.0, 0.5;
	Eigen::VectorXd rExps2(2);
	rExps2 << 0.0, 1.0;
	Eigen::VectorXd gaussExps2(2);
	gaussExps2 << 3.0, 2.0;

	// Add in wrong order: 2-term first, then 3-term
	ebs.SemiLocalEcpChannels.emplace_back(AzimuthalQuantumNumber::P, coeffs2, rExps2, gaussExps2);
	ebs.SemiLocalEcpChannels.emplace_back(AzimuthalQuantumNumber::P, coeffs3, rExps3, gaussExps3);
	REQUIRE(ebs.SemiLocalEcpChannels.size() == 2);

	REQUIRE_FALSE(ebs.IsInStandardRepresentation());

	ebs.StandardizeRepresentation();

	REQUIRE(ebs.IsInStandardRepresentation());
	REQUIRE(ebs.SemiLocalEcpChannels.size() == 1);
	REQUIRE(ebs.SemiLocalEcpChannels[0].GaussianTermCount() == 3 + 2);
}

TEST_CASE("StandardizeRepresentation sorts SemiLocalEcpChannels by gaussian exponent descending", "[ElementaryBasisSet]")
{
	ElementaryBasisSet ebs;

	// Two P channels with same term count but different exponents
	// Larger exponent should come first
	Eigen::VectorXd coeffs1(1);
	coeffs1 << 1.0;
	Eigen::VectorXd rExps1(1);
	rExps1 << 0.0;
	Eigen::VectorXd gaussExps1(1);
	gaussExps1 << 2.0;

	Eigen::VectorXd coeffs2(1);
	coeffs2 << 1.0;
	Eigen::VectorXd rExps2(1);
	rExps2 << 0.0;
	Eigen::VectorXd gaussExps2(1);
	gaussExps2 << 5.0;

	// Add in wrong order: smaller exponent first
	ebs.SemiLocalEcpChannels.emplace_back(AzimuthalQuantumNumber::P, coeffs1, rExps1, gaussExps1);
	ebs.SemiLocalEcpChannels.emplace_back(AzimuthalQuantumNumber::P, coeffs2, rExps2, gaussExps2);

	REQUIRE_FALSE(ebs.IsInStandardRepresentation());

	ebs.StandardizeRepresentation();

	REQUIRE(ebs.IsInStandardRepresentation());
	REQUIRE(ebs.SemiLocalEcpChannels.size() == 1);
	REQUIRE(ebs.SemiLocalEcpChannels[0].GaussianExponent(0) == 5.0);
	REQUIRE(ebs.SemiLocalEcpChannels[0].GaussianExponent(1) == 2.0);
}

TEST_CASE("StandardizeRepresentation sorts SemiLocalEcpChannels by coefficient descending", "[ElementaryBasisSet]")
{
	ElementaryBasisSet ebs;

	// Two P channels with same term count and exponent but different coefficients
	// Larger coefficient should come first
	Eigen::VectorXd coeffs1(1);
	coeffs1 << 0.5;
	Eigen::VectorXd rExps1(1);
	rExps1 << 0.0;
	Eigen::VectorXd gaussExps1(1);
	gaussExps1 << 3.0;

	Eigen::VectorXd coeffs2(1);
	coeffs2 << 1.5;
	Eigen::VectorXd rExps2(1);
	rExps2 << 0.0;
	Eigen::VectorXd gaussExps2(1);
	gaussExps2 << 3.0;

	// Add in wrong order: smaller coefficient first
	ebs.SemiLocalEcpChannels.emplace_back(AzimuthalQuantumNumber::P, coeffs1, rExps1, gaussExps1);
	ebs.SemiLocalEcpChannels.emplace_back(AzimuthalQuantumNumber::P, coeffs2, rExps2, gaussExps2);
	REQUIRE(ebs.SemiLocalEcpChannels.size() == 2);

	REQUIRE_FALSE(ebs.IsInStandardRepresentation());

	ebs.StandardizeRepresentation();

	REQUIRE(ebs.IsInStandardRepresentation());
	REQUIRE(ebs.SemiLocalEcpChannels.size() == 1);
	REQUIRE(ebs.SemiLocalEcpChannels[0].Coefficient(0) == 1.5);
	REQUIRE(ebs.SemiLocalEcpChannels[0].Coefficient(1) == 0.5);
}

TEST_CASE("StandardizeRepresentation sorts SemiLocalEcpChannels by r exponent descending", "[ElementaryBasisSet]")
{
	ElementaryBasisSet ebs;

	// Two P channels with same term count, exponent, and coefficient but different r exponents
	// Larger r exponent should come first
	Eigen::VectorXd coeffs(1);
	coeffs << 1.0;
	Eigen::VectorXd rExps1(1);
	rExps1 << 0.5;
	Eigen::VectorXd gaussExps(1);
	gaussExps << 3.0;

	Eigen::VectorXd rExps2(1);
	rExps2 << 2.0;

	// Add in wrong order: smaller r exponent first
	ebs.SemiLocalEcpChannels.emplace_back(AzimuthalQuantumNumber::P, coeffs, rExps1, gaussExps);
	ebs.SemiLocalEcpChannels.emplace_back(AzimuthalQuantumNumber::P, coeffs, rExps2, gaussExps);

	REQUIRE_FALSE(ebs.IsInStandardRepresentation());

	ebs.StandardizeRepresentation();

	REQUIRE(ebs.IsInStandardRepresentation());
	REQUIRE(ebs.SemiLocalEcpChannels.size() == 1);
	REQUIRE(ebs.SemiLocalEcpChannels[0].RExponent(0) == 2.0);
	REQUIRE(ebs.SemiLocalEcpChannels[0].RExponent(1) == 0.5);
}

TEST_CASE("StandardizeRepresentation concatenates same-l SemiLocalEcpChannels", "[ElementaryBasisSet]")
{
	ElementaryBasisSet ebs;

	// Two P channels that should be concatenated
	Eigen::VectorXd coeffs1(2);
	coeffs1 << 1.0, 0.5;
	Eigen::VectorXd rExps1(2);
	rExps1 << 0.0, 1.0;
	Eigen::VectorXd gaussExps1(2);
	gaussExps1 << 5.0, 4.0;

	Eigen::VectorXd coeffs2(1);
	coeffs2 << 0.25;
	Eigen::VectorXd rExps2(1);
	rExps2 << 2.0;
	Eigen::VectorXd gaussExps2(1);
	gaussExps2 << 3.0;

	ebs.SemiLocalEcpChannels.emplace_back(AzimuthalQuantumNumber::P, coeffs1, rExps1, gaussExps1);
	ebs.SemiLocalEcpChannels.emplace_back(AzimuthalQuantumNumber::P, coeffs2, rExps2, gaussExps2);

	REQUIRE_FALSE(ebs.IsInStandardRepresentation());

	ebs.StandardizeRepresentation();

	REQUIRE(ebs.IsInStandardRepresentation());
	REQUIRE(ebs.SemiLocalEcpChannels.size() == 1);
	REQUIRE(ebs.SemiLocalEcpChannels[0].AngularMomentum() == AzimuthalQuantumNumber::P);
	REQUIRE(ebs.SemiLocalEcpChannels[0].GaussianTermCount() == 3);
	REQUIRE(ebs.SemiLocalEcpChannels[0].GaussianExponent(0) == 5.0);
	REQUIRE(ebs.SemiLocalEcpChannels[0].GaussianExponent(1) == 4.0);
	REQUIRE(ebs.SemiLocalEcpChannels[0].GaussianExponent(2) == 3.0);
}

TEST_CASE("IsInStandardRepresentation with empty SemiLocalEcpChannels", "[ElementaryBasisSet]")
{
	ElementaryBasisSet ebs;

	REQUIRE(ebs.IsInStandardRepresentation());

	// Empty channels should be ignored
	Eigen::VectorXd coeffs(0);
	Eigen::VectorXd rExps(0);
	Eigen::VectorXd gaussExps(0);

	ebs.SemiLocalEcpChannels.emplace_back(AzimuthalQuantumNumber::P, coeffs, rExps, gaussExps);

	// standard representation should have empty channels removed
	REQUIRE_FALSE(ebs.IsInStandardRepresentation());
}

TEST_CASE("IsInStandardRepresentation with mixed empty and non-empty SemiLocalEcpChannels", "[ElementaryBasisSet]")
{
	ElementaryBasisSet ebs;

	Eigen::VectorXd coeffs(1);
	coeffs << 1.0;
	Eigen::VectorXd rExps(1);
	rExps << 0.0;
	Eigen::VectorXd gaussExps(1);
	gaussExps << 2.0;

	// Add non-empty channel
	ebs.SemiLocalEcpChannels.emplace_back(AzimuthalQuantumNumber::P, coeffs, rExps, gaussExps);

	// Now add an empty one
	Eigen::VectorXd empty(0);
	ebs.SemiLocalEcpChannels.emplace_back(AzimuthalQuantumNumber::D, empty, empty, empty);

	REQUIRE_FALSE(ebs.IsInStandardRepresentation());

	ebs.StandardizeRepresentation();

	// Empty channel should be removed
	REQUIRE(ebs.IsInStandardRepresentation());
	REQUIRE(ebs.SemiLocalEcpChannels.size() == 1);
	REQUIRE(ebs.SemiLocalEcpChannels[0].AngularMomentum() == AzimuthalQuantumNumber::P);
}
