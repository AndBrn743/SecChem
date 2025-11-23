//
// Created by Andy on 11/23/2025.
//

#include "../SecChem/Element.hpp"
#include <catch2/catch_test_macros.hpp>

using namespace SecChem;

TEST_CASE("Element convertion and basic properties", "CHEMISTRY")
{
	const Element neutron;  // should default to neutron
	CHECK(neutron.Symbol() == "n");
	CHECK(neutron.Name() == "Neutron");
	CHECK(neutron.AtomicNumber() == 0);
	CHECK(Element::Symbol2Element("n") == Element::Neutron);
	CHECK(Element::Symbol2Element("N") == Element::N);
	CHECK(Element::FuzzySymbol2Element("n") == Element::N);
	CHECK(Element::FuzzySymbol2Element("N") == Element::N);

	const Element iron = Element::Name2Element("iron");
	CHECK(iron.Symbol() == "Fe");
	CHECK(iron.Name() == "Iron");
	CHECK(iron.AtomicNumber() == 26);
	CHECK(iron.Period() == 4);
	CHECK(iron.Group() == 8);
	CHECK(iron.IsFromViiiGroup());
	CHECK(iron.IsFromDBlock());
	CHECK(iron.CharacteristicShell() == 3_Diffuse);
	CHECK(iron.NobleGasFromPeriod() == Element::Kr);
	CHECK(Element(Element::H).CharacteristicShell() == 1_Sharp);
	CHECK(Element(Element::He).CharacteristicShell() == 1_Sharp);
	CHECK(Element(Element::Li).CharacteristicShell() == 2_Sharp);
	CHECK(Element(Element::Ne).CharacteristicShell() == 2_Principal);
	CHECK_FALSE(Element(Element::Neutron).IsFromViiiGroup());
	CHECK_FALSE(Element(Element::Zn).IsFromViiiGroup());
	CHECK_FALSE(Element(Element::Hg).IsFromViiiGroup());
	CHECK_FALSE(Element(Element::Cs).IsFromViiiGroup());
	CHECK_FALSE(Element(Element::Cu).IsFromViiiGroup());
	CHECK_FALSE(Element(Element::Au).IsFromViiiGroup());
	CHECK(Element(Element::La).IsLanthanide());
	CHECK(Element(Element::Ce).IsLanthanide());
	CHECK(Element(Element::Pr).IsLanthanide());
	CHECK(Element(Element::Dy).IsLanthanide());
	CHECK(Element(Element::Yb).IsLanthanide());
	CHECK(Element(Element::Lu).IsLanthanide());
	CHECK_FALSE(Element(Element::Au).IsLanthanide());
}

TEST_CASE("Element's generated properties", "CHEMISTRY")
{
	const Element gold = Element::Au;
	CHECK(gold.ElectronConfiguration()
	      == AtomicElectronConfiguration({{1_Sharp, 2},
	                                      {2_Sharp, 2},
	                                      {2_Principal, 6},
	                                      {3_Sharp, 2},
	                                      {3_Diffuse, 10},
	                                      {3_Principal, 6},
	                                      {4_Sharp, 2},
	                                      {4_Principal, 6},
	                                      {4_Diffuse, 10},
	                                      {4_Fundamental, 14},
	                                      {5_Sharp, 2},
	                                      {5_Principal, 6},
	                                      {5_Diffuse, 10},
	                                      {6_Sharp, 1}}));
	CHECK(gold.MadelungElectronConfiguration()
	      == AtomicElectronConfiguration({{4_Sharp, 2},
	                                      {4_Principal, 6},
	                                      {4_Diffuse, 10},
	                                      {4_Fundamental, 14},
	                                      {5_Sharp, 2},
	                                      {5_Principal, 6},
	                                      {5_Diffuse, 9},
	                                      {6_Sharp, 2},
	                                      {1_Sharp, 2},
	                                      {2_Sharp, 2},
	                                      {2_Principal, 6},
	                                      {3_Sharp, 2},
	                                      {3_Principal, 6},
	                                      {3_Diffuse, 10}}));
	CHECK(gold.MadelungElectronConfiguration() == AtomicElectronConfiguration(gold.AtomicNumber()));
}

TEST_CASE("Element: Basic atomic number behavior")
{
	Element h(1);
	Element og(118);

	REQUIRE(h.AtomicNumber() == 1);
	REQUIRE(og.AtomicNumber() == 118);
}

TEST_CASE("Element: Block classification correct for common elements")
{
	Element h(1);    // 1s1
	Element c(6);    // 2p2
	Element fe(26);  // 3d6
	Element ce(58);  // 4f1

	REQUIRE(h.IsFromSBlock());
	REQUIRE_FALSE(h.IsFromPBlock());
	REQUIRE_FALSE(h.IsFromDBlock());
	REQUIRE_FALSE(h.IsFromFBlock());
	REQUIRE(h.CharacteristicOrbitalAzimuthalQuantumNumber() == 0);

	REQUIRE(c.IsFromPBlock());
	REQUIRE_FALSE(c.IsFromSBlock());
	REQUIRE(c.CharacteristicOrbitalAzimuthalQuantumNumber() == 1);

	REQUIRE(fe.IsFromDBlock());
	REQUIRE(fe.CharacteristicOrbitalAzimuthalQuantumNumber() == 2);

	REQUIRE(ce.IsFromFBlock());
	REQUIRE(ce.CharacteristicOrbitalAzimuthalQuantumNumber() == 3);
}

TEST_CASE("Element: Outermost shell principal quantum number")
{
	REQUIRE(Element(1).OuterMostShell() == 1_Sharp);        // H: 1s1
	REQUIRE(Element(10).OuterMostShell() == 2_Principal);   // Ne: 2p6
	REQUIRE(Element(18).OuterMostShell() == 3_Principal);   // Ar: 3p6
	REQUIRE(Element(36).OuterMostShell() == 4_Principal);   // Kr: 4p6
	REQUIRE(Element(54).OuterMostShell() == 5_Principal);   // Xe: 5p6
	REQUIRE(Element(86).OuterMostShell() == 6_Principal);   // Rn: 6p6
	REQUIRE(Element(118).OuterMostShell() == 7_Principal);  // Og: 7p6 (expected)
}

TEST_CASE("Element: Special-case Lawrencium Z=103")
{
	Element lr(103);

	// Modern, relativistic ground state:
	// Lr = [Rn] 5f14 7s2 7p1 → ℓ = 1 (p-block)
	REQUIRE(lr.AtomicNumber() == 103);
	REQUIRE(lr.CharacteristicOrbitalAzimuthalQuantumNumber() == 1);  // p

	// .IsFrom?Block() use classical blocking
	REQUIRE_FALSE(lr.IsFromPBlock());
	REQUIRE_FALSE(lr.IsFromSBlock());
	REQUIRE(lr.IsFromDBlock());
	REQUIRE_FALSE(lr.IsFromFBlock());

	REQUIRE(lr.OuterMostShell() == 7_Principal);  // 7p1 is the valence
}

TEST_CASE("Element: Light elements cross-check")
{
	struct Case
	{
		int Z;
		int l;
	};
	constexpr Case cases[] = {
	        {1, 0},  // H: 1s1
	        {2, 0},  // He: 1s2
	        {3, 0},  // Li: 2s1
	        {4, 0},  // Be: 2s2
	        {5, 1},  // B: 2p1
	        {6, 1},  // C: 2p2
	        {7, 1},  // N: 2p3
	        {8, 1},  // O: 2p4
	        {9, 1},  // F: 2p5
	        {10, 1}  // Ne: 2p6
	};

	for (auto& c : cases)
	{
		Element elem(c.Z);
		REQUIRE(elem.CharacteristicOrbitalAzimuthalQuantumNumber() == c.l);
	}
}

TEST_CASE("Element: First-row transition metals (d-block)")
{
	// 21–30: Sc → Zn
	for (int Z = 21; Z <= 30; ++Z)
	{
		Element e(Z);
		REQUIRE(e.IsFromDBlock());
		REQUIRE(e.CharacteristicOrbitalAzimuthalQuantumNumber() == 2);
	}
}

TEST_CASE("Element: Lanthanides (f-block)")
{
	// 57–71: La -> Lu

	{
		Element e(57);
		REQUIRE(e.CharacteristicOrbitalAzimuthalQuantumNumber() == 3);
		REQUIRE(e.IsFromFBlock());
	}

	for (int Z = 58; Z < 71; ++Z)
	{
		Element e(Z);
		REQUIRE(e.CharacteristicOrbitalAzimuthalQuantumNumber() == 3);
		REQUIRE(e.IsFromFBlock());
	}

	{
		Element e(71);
		REQUIRE(e.CharacteristicOrbitalAzimuthalQuantumNumber() == 2);
		REQUIRE(e.IsFromDBlock());
	}
}
