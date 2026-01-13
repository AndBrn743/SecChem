//
// Created by Andy on 1/1/2026.
//

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <SecChem/Atom.hpp>
#include <SecChem/Geometry.Input.hpp>
#include <SecChem/Utility/UnitOfMeasurement.hpp>

using namespace SecChem;
using namespace SecChem::Geometry::Input;
using namespace SecUtility::UnitOfMeasurement;

static constexpr auto Bohr2Bohr = [](const double l) { return l; };


TEST_CASE("ParseBasicZMatrixLine basic cases", "[InternalCoordinate][Parser][ZMatrix]")
{
	std::vector<Atom> atoms;

	SECTION("First atom")
	{
		std::vector<std::string> tokens{"C"};
		auto [elem, pos] = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, Degree2Radian);
		REQUIRE(elem == Element::C);
		REQUIRE(pos.isApprox(Eigen::Vector3d{0, 0, 0}));
		atoms.emplace_back(elem, pos);
	}

	SECTION("Second atom, bond length")
	{
		atoms.push_back(Atom(Element::H, Eigen::Vector3d{0, 0, 0}));
		std::vector<std::string> tokens{"O", "0", "1.5"};
		auto [elem, pos] = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, Degree2Radian);
		REQUIRE(elem == Element::O);
		REQUIRE(pos.isApprox(Eigen::Vector3d{0, 0, 1.5}));
		atoms.emplace_back(elem, pos);
	}

	SECTION("Second atom, bond length with conversion")
	{
		atoms.push_back(Atom(Element::H, Eigen::Vector3d{0, 0, 0}));
		std::vector<std::string> tokens{"O", "0", "1.5"};
		auto [elem, pos] =
		        ParseInternalCoordinateLine(tokens, atoms, [](const double l) { return 2 * l; }, Degree2Radian);
		REQUIRE(elem == Element::O);
		REQUIRE(pos.isApprox(Eigen::Vector3d{0, 0, 3}));
		atoms.emplace_back(elem, pos);
	}

	SECTION("Third atom, bond + angle")
	{
		atoms.push_back(Atom(Element::H, Eigen::Vector3d{0, 0, 0}));
		atoms.push_back(Atom(Element::O, Eigen::Vector3d{0, 0, 1.5}));
		std::vector<std::string> tokens{"H", "0", "1.0", "1", "90"};
		auto [elem, pos] = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, Degree2Radian);
		REQUIRE(elem == Element::H);
		// For 90 deg in XY plane
		REQUIRE(pos.isApprox(Eigen::Vector3d{1.0, 0.0, 0}, 1e-12));
	}

	SECTION("Third atom, bond + angle, with 1-based index")
	{
		atoms.push_back(Atom(Element::H, Eigen::Vector3d{0, 0, 0}));
		atoms.push_back(Atom(Element::O, Eigen::Vector3d{0, 0, 1.5}));
		std::vector<std::string> tokens{"H", "1", "1.0", "2", "90"};
		auto [elem, pos] = ParseInternalCoordinateLine<1>(tokens, atoms, Bohr2Bohr, Degree2Radian);
		REQUIRE(elem == Element::H);
		// For 90 deg in XY plane
		REQUIRE(pos.isApprox(Eigen::Vector3d{1.0, 0.0, 0}, 1e-12));
	}

	SECTION("Third atom, bond + angle, uncommon direction")
	{
		atoms.push_back(Atom(Element::H, Eigen::Vector3d{0, 0, 0}));
		atoms.push_back(Atom(Element::O, Eigen::Vector3d{0, 0, 1.5}));
		std::vector<std::string> tokens{"H", "1", "1.0", "0", "120"};
		auto [elem, pos] = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, Degree2Radian);
		REQUIRE(elem == Element::H);
		// For 90 deg in XY plane
		REQUIRE(pos.isApprox(Eigen::Vector3d{std::sin(Degree2Radian(60)), 0.0, 2}, 1e-12));
	}

	SECTION("Fourth atom, bond + angle + dihedral")
	{
		atoms.push_back(Atom(Element::H, Eigen::Vector3d{0, 0, 0}));
		atoms.push_back(Atom(Element::O, Eigen::Vector3d{0, 0, 1.5}));
		atoms.push_back(Atom(Element::H, Eigen::Vector3d{1.0, 0.0, 1.5}));
		std::vector<std::string> tokens{"C", "1", "1.2", "0", "90", "2", "180"};
		auto [elem, pos] = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, Degree2Radian);
		REQUIRE(elem == Element::C);
		// Just sanity check: distance to atom 1 should be close to 1.2
		REQUIRE((pos - atoms[1].Position()).norm() == Catch::Approx(1.2).epsilon(1e-12));
	}
}

TEST_CASE("ParseBasicZMatrixLine negative tests", "[InternalCoordinate][Parser][ZMatrix]")
{
	SECTION("Too few tokens for first atom")
	{
		std::vector<Atom> atoms;
		std::vector<std::string> tokens{};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, Degree2Radian), std::runtime_error);
	}

	SECTION("Invalid bond length")
	{
		std::vector atoms{Atom(Element::H, Eigen::Vector3d{0, 0, 0})};
		std::vector<std::string> tokens{"O", "0", "-1.0"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, Degree2Radian), std::runtime_error);
	}

	SECTION("Invalid atom index")
	{
		std::vector atoms{Atom(Element::H, Eigen::Vector3d{0, 0, 0})};
		std::vector<std::string> tokens{"O", "1", "1.0"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, Degree2Radian), std::out_of_range);
	}

	SECTION("Too few tokens for third atom")
	{
		std::vector atoms{Atom(Element::H, Eigen::Vector3d{0, 0, 0}), Atom(Element::O, Eigen::Vector3d{0, 0, 1.5})};
		std::vector<std::string> tokens{"H", "0", "1.0", "1"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, Degree2Radian), std::runtime_error);
	}
}

TEST_CASE("Parse full Z-matrix input", "[InternalCoordinate][Parser][ZMatrix][integration]")
{
	std::vector<std::string> zMatrixLines = {
	        "C",                       // first atom
	        "H 0 1.089",               // second atom, bond length
	        "H 0 1.089 1 109.5",       // third atom, bond + angle
	        "O 1 1.430 0 104.5 2 180"  // fourth atom, bond + angle + dihedral
	};

	std::vector<Atom> atoms;

	for (const auto& line : zMatrixLines)
	{
		std::vector<std::string> tokens;
		std::istringstream iss(line);
		std::string tok;
		while (iss >> tok)
		{
			tokens.push_back(tok);
		}

		auto [element, pos] = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, Degree2Radian);
		atoms.emplace_back(element, pos);
	}

	REQUIRE(atoms.size() == zMatrixLines.size());

	// Check first atom is at origin
	REQUIRE(atoms[0].Position().isApprox(Eigen::Vector3d{0, 0, 0}));

	// Check second atom is on Z-axis
	REQUIRE(atoms[1].Position().isApprox(Eigen::Vector3d{0, 0, 1.089}, 1e-12));

	// Check third atom distance from first atom ~ bond length
	REQUIRE(atoms[2].DistanceTo(atoms[0]) == Catch::Approx(1.089).epsilon(1e-12));

	// Check fourth atom distance from second atom ~ bond length
	REQUIRE(atoms[3].DistanceTo(atoms[1]) == Catch::Approx(1.430).epsilon(1e-12));

	REQUIRE(atoms[3].DistanceTo(atoms[0]) == Catch::Approx(2.002651181).margin(1e-8));
	REQUIRE(atoms[3].DistanceTo(atoms[2]) == Catch::Approx(3.01512618).margin(1e-8));

	// Optional: check approximate angles using vector dot products
	auto angle = [](const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c)
	{
		const Eigen::Vector3d ab = (a - b).normalized();
		const Eigen::Vector3d cb = (c - b).normalized();
		return std::acos(ab.dot(cb));
	};

	REQUIRE(angle(atoms[0].Position(), atoms[2].Position(), atoms[1].Position())
	        == Catch::Approx(Degree2Radian((180 - 109.5) / 2)).epsilon(1e-3));
	REQUIRE(angle(atoms[0].Position(), atoms[1].Position(), atoms[3].Position())
	        == Catch::Approx(Degree2Radian(104.5)).epsilon(1e-3));
}

TEST_CASE("Parse full Z-matrix input 2", "[InternalCoordinate][Parser][ZMatrix][integration]")
{
	std::vector<std::string> zMatrixLines = {
	        "C",                      // first atom
	        "H 0 1.089",              // second atom, bond length
	        "H 0 1.089 1 109.5",      // third atom, bond + angle
	        "O 1 1.430 0 104.5 2 90"  // fourth atom, bond + angle + dihedral
	};

	std::vector<Atom> atoms;

	for (const auto& line : zMatrixLines)
	{
		std::vector<std::string> tokens;
		std::istringstream iss(line);
		std::string tok;
		while (iss >> tok)
		{
			tokens.push_back(tok);
		}

		auto [element, pos] = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, Degree2Radian);
		atoms.emplace_back(element, pos);
	}

	REQUIRE(atoms.size() == zMatrixLines.size());

	// Check first atom is at origin
	REQUIRE(atoms[0].Position().isApprox(Eigen::Vector3d{0, 0, 0}));

	// Check second atom is on Z-axis
	REQUIRE(atoms[1].Position().isApprox(Eigen::Vector3d{0, 0, 1.089}, 1e-12));

	// Check third atom distance from first atom ~ bond length
	REQUIRE(atoms[2].DistanceTo(atoms[0]) == Catch::Approx(1.089).epsilon(1e-12));

	// Check fourth atom distance from second atom ~ bond length
	REQUIRE(atoms[3].DistanceTo(atoms[1]) == Catch::Approx(1.430).epsilon(1e-12));

	REQUIRE(atoms[3].DistanceTo(atoms[0]) == Catch::Approx(2.002651181).margin(1e-8));
	REQUIRE(atoms[3].DistanceTo(atoms[2]) == Catch::Approx(2.499721270).margin(1e-8));

	// Optional: check approximate angles using vector dot products
	auto angle = [](const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c)
	{
		const Eigen::Vector3d ab = (a - b).normalized();
		const Eigen::Vector3d cb = (c - b).normalized();
		return std::acos(ab.dot(cb));
	};

	REQUIRE(angle(atoms[0].Position(), atoms[2].Position(), atoms[1].Position())
	        == Catch::Approx(Degree2Radian((180 - 109.5) / 2)).epsilon(1e-3));
	REQUIRE(angle(atoms[0].Position(), atoms[1].Position(), atoms[3].Position())
	        == Catch::Approx(Degree2Radian(104.5)).epsilon(1e-3));
}

TEST_CASE("Parse full Z-matrix input 3", "[InternalCoordinate][Parser][ZMatrix][integration]")
{
	std::vector<std::string> zMatrixLines = {
	        "C",                       // first atom
	        "H 0 1.089",               // second atom, bond length
	        "H 0 1.089 1 109.5",       // third atom, bond + angle
	        "O 1 1.430 0 104.5 2 120"  // fourth atom, bond + angle + dihedral
	};

	std::vector<Atom> atoms;

	for (const auto& line : zMatrixLines)
	{
		std::vector<std::string> tokens;
		std::istringstream iss(line);
		std::string tok;
		while (iss >> tok)
		{
			tokens.push_back(tok);
		}

		auto [element, pos] = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, Degree2Radian);
		atoms.emplace_back(element, pos);
	}

	REQUIRE(atoms.size() == zMatrixLines.size());

	// Check first atom is at origin
	REQUIRE(atoms[0].Position().isApprox(Eigen::Vector3d{0, 0, 0}));

	// Check second atom is on Z-axis
	REQUIRE(atoms[1].Position().isApprox(Eigen::Vector3d{0, 0, 1.089}, 1e-12));

	// Check third atom distance from first atom ~ bond length
	REQUIRE(atoms[2].DistanceTo(atoms[0]) == Catch::Approx(1.089).epsilon(1e-12));

	// Check fourth atom distance from second atom ~ bond length
	REQUIRE(atoms[3].DistanceTo(atoms[1]) == Catch::Approx(1.430).epsilon(1e-12));

	REQUIRE(atoms[3].DistanceTo(atoms[0]) == Catch::Approx(2.002651181).margin(1e-8));
	REQUIRE(atoms[3].DistanceTo(atoms[2]) == Catch::Approx(2.76943969).margin(1e-8));

	// Optional: check approximate angles using vector dot products
	auto angle = [](const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c)
	{
		const Eigen::Vector3d ab = (a - b).normalized();
		const Eigen::Vector3d cb = (c - b).normalized();
		return std::acos(ab.dot(cb));
	};

	REQUIRE(angle(atoms[0].Position(), atoms[2].Position(), atoms[1].Position())
	        == Catch::Approx(Degree2Radian((180 - 109.5) / 2)).epsilon(1e-3));
	REQUIRE(angle(atoms[0].Position(), atoms[1].Position(), atoms[3].Position())
	        == Catch::Approx(Degree2Radian(104.5)).epsilon(1e-3));
}
