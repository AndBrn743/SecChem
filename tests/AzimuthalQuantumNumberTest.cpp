//
// Created by Andy on 11/22/2025.
//

// test_azimuthal_quantum_number.cpp
#define CATCH_CONFIG_MAIN
#include "../SecChem/AzimuthalQuantumNumber.hpp"  // your class header
#include <catch2/catch_test_macros.hpp>
#include <sstream>

using namespace SecChem;

TEST_CASE("AzimuthalQuantumNumber: construction", "[AzimuthalQuantumNumber]")
{
	using AQ = AzimuthalQuantumNumber;

	// Implicit from static constants
	AQ s = AQ::S;
	AQ p = AQ::P;

	REQUIRE(static_cast<int>(s) == 0);
	REQUIRE(static_cast<int>(p) == 1);

	// Explicit from int
	AQ d(2);
	REQUIRE(static_cast<int>(d) == 2);

	// From char labels
	AQ f('f');
	REQUIRE(f == AQ::F);

	// Fuzzy label (uppercase)
	AQ g('G');
	REQUIRE(static_cast<int>(g) == 4);
}

TEST_CASE("AzimuthalQuantumNumber: comparisons", "[AzimuthalQuantumNumber]")
{
	using AQ = AzimuthalQuantumNumber;

	AQ s = AQ::S;
	AQ p = AQ::P;
	AQ d = AQ::D;

	REQUIRE(s < p);
	REQUIRE(p > s);
	REQUIRE(d >= p);
	REQUIRE(d != p);
	REQUIRE(p == AQ::P);
}

TEST_CASE("AzimuthalQuantumNumber: increment/decrement", "[AzimuthalQuantumNumber]")
{
	using AQ = AzimuthalQuantumNumber;

	AQ s = AQ::S;
	++s;
	REQUIRE(s == AQ::P);

	AQ d = AQ::D;
	d--;
	REQUIRE(d == AQ::P);

	AQ f = AQ::F;
	f += AQ(2);  // F + 2 -> H
	REQUIRE(f == AQ::H);

	AQ g = AQ::G;
	g -= AQ(2);                         // G - 2 -> E? but limited by max/min
	REQUIRE(static_cast<int>(g) == 2);  // D
}

TEST_CASE("AzimuthalQuantumNumber: Min/Max Magnetic Quantum Numbers", "[AzimuthalQuantumNumber]")
{
	using AQ = AzimuthalQuantumNumber;

	AQ s = AQ::S;
	REQUIRE(s.MagneticQuantumNumberCount() == 1);
	REQUIRE(s.MinMagneticQuantumNumber() == 0);
	REQUIRE(s.MaxMagneticQuantumNumber() == 0);

	AQ p = AQ::P;
	REQUIRE(p.MagneticQuantumNumberCount() == 3);
	REQUIRE(p.MinMagneticQuantumNumber() == -1);
	REQUIRE(p.MaxMagneticQuantumNumber() == 1);
}

TEST_CASE("AzimuthalQuantumNumber: capacity", "[AzimuthalQuantumNumber]")
{
	using AQ = AzimuthalQuantumNumber;

	AQ s = AQ::S;
	AQ p = AQ::P;
	AQ d = AQ::D;

	REQUIRE(s.Capacity() == 2);   // 2 electrons
	REQUIRE(p.Capacity() == 6);   // 3 magnetic numbers * 2
	REQUIRE(d.Capacity() == 10);  // 5 magnetic numbers * 2
}

TEST_CASE("AzimuthalQuantumNumber: Label and Name", "[AzimuthalQuantumNumber]")
{
	using AQ = AzimuthalQuantumNumber;

	AQ s = AQ::S;
	AQ p = AQ::P;
	AQ f = AQ::F;
	AQ h = AQ::H;

	REQUIRE(s.Label() == 's');
	REQUIRE(p.Label() == 'p');
	REQUIRE(f.Label() == 'f');
	REQUIRE(h.Label() == 'h');

	REQUIRE(s.Name() == "Sharp");
	REQUIRE(p.Name() == "Principal");
	REQUIRE(f.Name() == "Fundamental");
	REQUIRE(h.Name() == "ShellH");
}

TEST_CASE("AzimuthalQuantumNumber: ostream/istream", "[AzimuthalQuantumNumber]")
{
	using AQ = AzimuthalQuantumNumber;

	std::ostringstream oss;
	AQ d = AQ::D;
	oss << d;
	REQUIRE(oss.str() == "d");

	std::istringstream iss("f");
	AQ f;
	iss >> f;
	REQUIRE(f == AQ::F);
}
