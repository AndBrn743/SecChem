//
// Created by Andy on 12/31/2025.
//

#include <Eigen/Core>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <SecChem/Atom.hpp>

using namespace SecChem;

TEST_CASE("Atom: basic construction", "[Atom]")
{
	const Eigen::Vector3d pos{1.0, 2.0, 3.0};
	Atom atom{Element::H, pos};

	REQUIRE(atom.Element() == SecChem::Element::H);
	REQUIRE(atom.Position().isApprox(pos));
	REQUIRE(atom.Mass() == SecChem::Element{SecChem::Element::H}.Mass());

	REQUIRE_FALSE(atom.IsWithFiniteNuclear());
	REQUIRE_FALSE(atom.IsFrozen());
}

TEST_CASE("Atom: finite nuclear model tags are mutually exclusive", "[Atom][AtomTag]")
{
	Atom atom{Element::C, Eigen::Vector3d::Zero()};

	atom.AddTags(AtomTag::GaussianFiniteNuclear);
	REQUIRE(atom.IsWithGaussianFiniteNuclear());
	REQUIRE_FALSE(atom.IsWithThomasFermiFiniteNuclear());

	atom.AddTags(AtomTag::ThomasFermiFiniteNuclear);
	REQUIRE_FALSE(atom.IsWithGaussianFiniteNuclear());
	REQUIRE(atom.IsWithThomasFermiFiniteNuclear());
}

TEST_CASE("Atom: invalid tag combination throws", "[Atom][AtomTag]")
{
	const auto invalidTags = AtomTag::GaussianFiniteNuclear | AtomTag::ThomasFermiFiniteNuclear;

	REQUIRE_THROWS_AS(SecChem::Atom(SecChem::Element::H, Eigen::Vector3d::Zero(), invalidTags), std::invalid_argument);
}

TEST_CASE("Atom: DistanceTo computes Euclidean distance", "[Atom]")
{
	Atom a{Element::H, Eigen::Vector3d{0.0, 0.0, 0.0}};

	Atom b{Element::H, Eigen::Vector3d{1.0, 2.0, 2.0}};

	REQUIRE(a.DistanceTo(b) == Catch::Approx(3.0));
}

TEST_CASE("Atom: custom mass and nuclear radius factories", "[Atom]")
{
	const Eigen::Vector3d pos{0.1, 0.2, 0.3};

	auto atom = Atom::AtomWithMassAndNuclearRadius(Element::C, pos, 13.5, 0.42);

	REQUIRE(atom.Element() == SecChem::Element::C);
	REQUIRE(atom.Position().isApprox(pos));
	REQUIRE(atom.Mass() == Catch::Approx(13.5));
	REQUIRE(atom.NuclearRadius() == Catch::Approx(0.42));
}

TEST_CASE("Atom: exact equality", "[Atom][Equality]")
{
	const Eigen::Vector3d pos{1.0, 1.0, 1.0};

	Atom a{Element::H, pos};
	Atom b{Element::H, pos};

	REQUIRE(a == b);
	REQUIRE_FALSE(a != b);
}

TEST_CASE("Atom: tolerance-based equality (position)", "[Atom][Equality]")
{
	Atom a{Element::C, Eigen::Vector3d{1.0, 1.0, 1.0}};

	Atom b{Element::C, Eigen::Vector3d{1.0 + 1e-13, 1.0, 1.0}};

	REQUIRE_FALSE(a == b);
	REQUIRE(a.EqualsTo(b, 1e-12));
}

TEST_CASE("Atom: tag mismatch breaks equality", "[Atom][Equality]")
{
	Atom a{Element::H, Eigen::Vector3d::Zero()};

	Atom b = a;
	b.AddTags(AtomTag::Frozen);

	REQUIRE(a != b);
}

TEST_CASE("Atom: different elements are never equal", "[Atom][Equality]")
{
	Atom h{Element::H, Eigen::Vector3d::Zero()};

	Atom c{Element::C, Eigen::Vector3d::Zero()};

	REQUIRE(h != c);
}
