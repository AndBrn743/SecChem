//
// Created by Andy on 11/22/2025.
//

#include <SecChem/ElectronicSubShell.hpp>
#include <catch2/catch_test_macros.hpp>
#include <sstream>

using namespace SecChem;

TEST_CASE("AzimuthalQuantumNumber basic properties", "[AzimuthalQuantumNumber]")
{
	AzimuthalQuantumNumber s = AzimuthalQuantumNumber::S;
	AzimuthalQuantumNumber p = AzimuthalQuantumNumber::P;
	AzimuthalQuantumNumber d = AzimuthalQuantumNumber::D;

	REQUIRE(s.Value() == 0);
	REQUIRE(p.Value() == 1);
	REQUIRE(d.Value() == 2);

	REQUIRE(s.MinPrincipalQuantumNumber() == 1);
	REQUIRE(p.MinPrincipalQuantumNumber() == 2);
	REQUIRE(d.MinPrincipalQuantumNumber() == 3);

	REQUIRE(s.MagneticQuantumNumberCount() == 1);
	REQUIRE(p.MagneticQuantumNumberCount() == 3);
	REQUIRE(d.MagneticQuantumNumberCount() == 5);

	AzimuthalQuantumNumber l('f');
	REQUIRE(l.Value() == 3);
	REQUIRE(l.Label() == 'f');

	AzimuthalQuantumNumber l2(4);
	REQUIRE(l2.Value() == 4);
	REQUIRE(l2.Label() == 'g');
}

TEST_CASE("ElectronicSubShell construction", "[ElectronicSubShell]")
{
	ElectronicSubShell s1(1, AzimuthalQuantumNumber::S);
	REQUIRE(s1.PrincipalQuantumNumber() == 1);
	REQUIRE(s1.SubshellLabel() == 's');
	REQUIRE(s1.IsSharp());

	ElectronicSubShell p2(2, 'p');
	REQUIRE(p2.PrincipalQuantumNumber() == 2);
	REQUIRE(p2.SubshellLabel() == 'p');
	REQUIRE(p2.IsPrincipal());

	ElectronicSubShell d3(3, 2);  // ℓ=2
	REQUIRE(d3.PrincipalQuantumNumber() == 3);
	REQUIRE(d3.SubshellLabel() == 'd');
	REQUIRE(d3.IsDiffuse());
}

TEST_CASE("ElectronicSubShell comparisons", "[ElectronicSubShell]")
{
	ElectronicSubShell s1(1, 's');
	ElectronicSubShell s2(2, 's');
	ElectronicSubShell p2(2, 'p');

	REQUIRE(s1 < s2);
	REQUIRE(s2 > s1);
	REQUIRE(s2 != p2);
	REQUIRE((s1 <= s1));
	REQUIRE((p2 >= s2));
	REQUIRE((s1 == s1));
}

TEST_CASE("ElectronicSubShell increment/decrement", "[ElectronicSubShell]")
{
	ElectronicSubShell shell(2, 's');  // 2s
	shell++;
	REQUIRE(shell.PrincipalQuantumNumber() == 2);
	REQUIRE(shell.SubshellLabel() == 'p');

	++shell;
	REQUIRE(shell.PrincipalQuantumNumber() == 3);
	REQUIRE(shell.SubshellLabel() == 's');

	shell--;
	REQUIRE(shell.PrincipalQuantumNumber() == 2);
	REQUIRE(shell.SubshellLabel() == 'p');

	--shell;
	REQUIRE(shell.PrincipalQuantumNumber() == 2);
	REQUIRE(shell.SubshellLabel() == 's');
}

TEST_CASE("ElectronicSubShell string conversion and streams", "[ElectronicSubShell]")
{
	ElectronicSubShell s1(1, 's');
	std::stringstream ss;
	ss << s1;
	REQUIRE(ss.str() == "1s");

	ElectronicSubShell readShell;
	ss >> readShell;
	REQUIRE(readShell == s1);

	REQUIRE(s1.ToString() == "1s");
	REQUIRE(s1.SubshellName() == "Sharp");
}

TEST_CASE("ElectronicSubShell Slater effective principal quantum number", "[ElectronicSubShell]")
{
	ElectronicSubShell s1(1, 's');
	REQUIRE(s1.SlaterEffectivePrincipalQuantumNumber() == 1.0);

	ElectronicSubShell n4f(4, 'f');
	REQUIRE(n4f.SlaterEffectivePrincipalQuantumNumber() == 3.7);

	ElectronicSubShell n7p(7, 'p');
	REQUIRE(n7p.SlaterEffectivePrincipalQuantumNumber() == 4.3);
}

TEST_CASE("ElectronicSubShell validity", "[ElectronicSubShell]")
{
	ElectronicSubShell valid(3, 'd');
	REQUIRE(valid.IsValid());
}
