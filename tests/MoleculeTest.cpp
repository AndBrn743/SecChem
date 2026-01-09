//
// Created by Andy on 1/9/2026.
//

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <SecChem/Molecule.hpp>

using namespace SecChem;

TEST_CASE("Molecule basic construction and AtomCount", "[molecule]")
{
	Molecule mol{Atom{Element::H, Eigen::Vector3d{0, 0, 0}},
	             Atom{Element::O, Eigen::Vector3d{0, 0, 1}},
	             Atom{Element::H, Eigen::Vector3d{1, 0, 0}}};

	REQUIRE(mol.AtomCount() == 3);
}

TEST_CASE("Molecule operator[] matches iteration order", "[molecule]")
{
	Molecule mol{Atom{Element::C, Eigen::Vector3d{0, 0, 0}},
	             Atom{Element::H, Eigen::Vector3d{1, 0, 0}},
	             Atom{Element::H, Eigen::Vector3d{-1, 0, 0}}};

	std::size_t i = 0;
	for (const Atom& atom : mol)
	{
		REQUIRE(&atom == &mol[i]);
		++i;
	}

	REQUIRE(i == mol.AtomCount());
}

TEST_CASE("IndexOf returns correct index for referenced atom", "[molecule]")
{
	Molecule mol{Atom{Element::N, Eigen::Vector3d{0, 0, 0}},
	             Atom{Element::H, Eigen::Vector3d{0, 0, 1}},
	             Atom{Element::H, Eigen::Vector3d{1, 0, 0}},
	             Atom{Element::H, Eigen::Vector3d{0, 1, 0}}};

	const Atom& a = mol[2];
	REQUIRE(mol.IndexOf(a) == 2);
}

TEST_CASE("IndexOf finds value-equal atom not belonging to molecule", "[molecule]")
{
	Molecule mol{Atom{Element::He, Eigen::Vector3d{0, 0, 0}}};

	Atom external{Element::He, Eigen::Vector3d{0, 0, 0}};  // same value, different object

	REQUIRE(mol.IndexOf(external) == 0);
}

TEST_CASE("IndexOf throws when atom is not found", "[molecule]")
{
	Molecule mol{Atom{Element::Li, Eigen::Vector3d{0, 0, 0}}};

	Atom external{Element::Be, Eigen::Vector3d{1, 0, 0}};

	REQUIRE_THROWS_AS(mol.IndexOf(external), std::out_of_range);
}

TEST_CASE("IndexOfUnchecked works for contained atoms", "[molecule]")
{
	Molecule mol{Atom{Element::C, Eigen::Vector3d{0, 0, 0}}, Atom{Element::O, Eigen::Vector3d{0, 0, 1}}};

	const Atom& atom = mol[1];
	REQUIRE(mol.IndexOfUnchecked(atom) == 1);  // We intentionally do not test UB cases — unchecked means unchecked
}

TEST_CASE("Contains detects both identity and value equality", "[molecule]")
{
	Molecule mol{Atom{Element::Ne, Eigen::Vector3d{0, 0, 0}}};

	const Atom& a = mol[0];
	Atom copy{Element::Ne, Eigen::Vector3d{0, 0, 0}};
	Atom other{Element::Ar, Eigen::Vector3d{0, 0, 0}};

	REQUIRE(mol.Contains(a));
	REQUIRE(mol.Contains(copy));
	REQUIRE_FALSE(mol.Contains(other));
}

TEST_CASE("SharedMolecule preserves atom identity", "[shared_molecule]")
{
	Molecule mol{Atom{Element::Fe, Eigen::Vector3d{0, 0, 0}}, Atom{Element::Fe, Eigen::Vector3d{0, 0, 2}}};

	SharedMolecule shared{mol};

	const Atom& a0 = mol[0];
	const Atom& b0 = shared[0];

	REQUIRE(&a0 != &b0);  // different storage
	REQUIRE(mol.IndexOf(a0) == 0);
	REQUIRE(shared.IndexOf(b0) == 0);
}

TEST_CASE("Mass and nuclear charge are additive", "[molecule]")
{
	Molecule mol{Atom{Element::H, Eigen::Vector3d{0, 0, 0}}, Atom{Element::H, Eigen::Vector3d{0, 0, 1}}};

	REQUIRE(mol.TotalNuclearCharge() == 2);
	REQUIRE(mol.Mass() > 0.0);
}

TEST_CASE("NuclearRepulsionEnergy is positive and symmetric", "[molecule]")
{
	Molecule mol{Atom{Element::H, Eigen::Vector3d{0, 0, 0}}, Atom{Element::H, Eigen::Vector3d{0, 0, 1}}};

	const double e = mol.NuclearRepulsionEnergy();

	REQUIRE(e > 0.0);

	Molecule reversed{mol[1], mol[0]};

	REQUIRE(e == Catch::Approx(reversed.NuclearRepulsionEnergy()));
}
