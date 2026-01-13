//
// Created by Claude on 1/10/2026.
//

// ReSharper disable CppUseStructuredBinding
#include <SecChem/Geometry.Input.hpp>
#include <SecChem/Utility/UnitOfMeasurement.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <sstream>
#include <string>
#include <vector>

using Catch::Approx;
using SecChem::Atom;
using SecChem::Element;
using namespace SecChem::Geometry::Input;
namespace UnitOfMeasurement = SecUtility::UnitOfMeasurement;

TEST_CASE("ParseBasicCartesianCoordinateLine: Valid single atom", "[Geometry][Parser]")
{
	const std::vector<std::string> tokens = {"H", "0.0", "0.0", "0.0"};
	const auto result = ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius);

	CHECK(result.Element == Element::H);
	CHECK(result.Position.x() == Approx(0.0));
	CHECK(result.Position.y() == Approx(0.0));
	CHECK(result.Position.z() == Approx(0.0));
}

TEST_CASE("ParseBasicCartesianCoordinateLine: Valid atom with non-zero position", "[Geometry][Parser]")
{
	const std::vector<std::string> tokens = {"C", "1.0", "2.0", "3.0"};
	const auto result = ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius);

	CHECK(result.Element == Element::C);
	// 1 Angstrom = 1.889726125... Bohr
	constexpr double angstrom2bohr = 1.8897261254578281;
	CHECK(result.Position.x() == Approx(1.0 * angstrom2bohr));
	CHECK(result.Position.y() == Approx(2.0 * angstrom2bohr));
	CHECK(result.Position.z() == Approx(3.0 * angstrom2bohr));
}

TEST_CASE("ParseBasicCartesianCoordinateLine: Various elements", "[Geometry][Parser]")
{
	SECTION("Oxygen")
	{
		const std::vector<std::string> tokens = {"O", "0.5", "1.5", "2.5"};
		const auto result = ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius);
		CHECK(result.Element == Element::O);
	}

	SECTION("Nitrogen")
	{
		const std::vector<std::string> tokens = {"N", "-0.5", "-1.5", "-2.5"};
		const auto result = ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius);
		CHECK(result.Element == Element::N);
	}

	SECTION("Iron")
	{
		const std::vector<std::string> tokens = {"Fe", "0.0", "0.0", "0.0"};
		const auto result = ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius);
		CHECK(result.Element == Element::Fe);
	}

	SECTION("Gold")
	{
		const std::vector<std::string> tokens = {"Au", "1.0", "0.0", "0.0"};
		const auto result = ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius);
		CHECK(result.Element == Element::Au);
	}
}

TEST_CASE("ParseBasicCartesianCoordinateLine: Negative coordinates", "[Geometry][Parser]")
{
	const std::vector<std::string> tokens = {"H", "-1.0", "-2.0", "-3.0"};
	const auto result = ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius);

	CHECK(result.Element == Element::H);
	constexpr double angstrom2bohr = 1.8897261254578281;
	CHECK(result.Position.x() == Approx(-1.0 * angstrom2bohr));
	CHECK(result.Position.y() == Approx(-2.0 * angstrom2bohr));
	CHECK(result.Position.z() == Approx(-3.0 * angstrom2bohr));
}

TEST_CASE("ParseBasicCartesianCoordinateLine: Scientific notation", "[Geometry][Parser]")
{
	const std::vector<std::string> tokens = {"H", "1.0e-10", "2.0e-5", "3.0e5"};
	const auto result = ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius);

	CHECK(result.Element == Element::H);
	CHECK(result.Position.x() == Approx(1.0e-10 * 1.8897261254578281));
	CHECK(result.Position.y() == Approx(2.0e-5 * 1.8897261254578281));
	CHECK(result.Position.z() == Approx(3.0e5 * 1.8897261254578281));
}

TEST_CASE("ParseBasicCartesianCoordinateLine: Token iterator overload", "[Geometry][Parser]")
{
	const std::vector<std::string> tokens = {"C", "1.0", "2.0", "3.0"};
	const auto result = ParseBasicCartesianCoordinateLine(tokens.begin(), UnitOfMeasurement::Angstrom2BohrRadius);

	CHECK(result.Element == Element::C);
	constexpr double angstrom2bohr = 1.8897261254578281;
	CHECK(result.Position.x() == Approx(1.0 * angstrom2bohr));
	CHECK(result.Position.y() == Approx(2.0 * angstrom2bohr));
	CHECK(result.Position.z() == Approx(3.0 * angstrom2bohr));
}

TEST_CASE("ParseBasicCartesianCoordinateLine: Different length units", "[Geometry][Parser]")
{
	const std::vector<std::string> tokens = {"H", "1.0", "2.0", "3.0"};

	SECTION("Angstrom to Bohr")
	{
		const auto result = ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius);
		constexpr double angstrom2bohr = 1.8897261254578281;
		CHECK(result.Position.x() == Approx(1.0 * angstrom2bohr));
	}

	SECTION("No conversion (identity)")
	{
		const auto result = ParseBasicCartesianCoordinateLine(tokens, [](double x) { return x; });
		CHECK(result.Position.x() == Approx(1.0));
		CHECK(result.Position.y() == Approx(2.0));
		CHECK(result.Position.z() == Approx(3.0));
	}
}

//--------------------------------------------------------------------------------------------------------------------//
// Negative tests
//--------------------------------------------------------------------------------------------------------------------//

TEST_CASE("ParseBasicCartesianCoordinateLine: Wrong number of tokens - too few", "[Geometry][Parser][negative]")
{
	SECTION("Empty token list")
	{
		const std::vector<std::string> tokens = {};
		REQUIRE_THROWS_AS(ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius),
		                  std::runtime_error);
	}

	SECTION("Only element symbol")
	{
		const std::vector<std::string> tokens = {"H"};
		REQUIRE_THROWS_AS(ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius),
		                  std::runtime_error);
	}

	SECTION("Element + one coordinate")
	{
		const std::vector<std::string> tokens = {"H", "1.0"};
		REQUIRE_THROWS_AS(ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius),
		                  std::runtime_error);
	}

	SECTION("Element + two coordinates")
	{
		const std::vector<std::string> tokens = {"H", "1.0", "2.0"};
		REQUIRE_THROWS_AS(ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius),
		                  std::runtime_error);
	}
}

TEST_CASE("ParseBasicCartesianCoordinateLine: Wrong number of tokens - too many", "[Geometry][Parser][negative]")
{
	const std::vector<std::string> tokens = {"H", "0.0", "0.0", "0.0", "extra"};
	REQUIRE_THROWS_AS(ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius),
	                  std::runtime_error);
}

TEST_CASE("ParseBasicCartesianCoordinateLine: Invalid element symbol", "[Geometry][Parser][negative]")
{
	SECTION("Completely invalid symbol")
	{
		const std::vector<std::string> tokens = {"Xx", "0.0", "0.0", "0.0"};
		REQUIRE_THROWS_AS(ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius),
		                  std::runtime_error);
	}

	SECTION("Lowercase element symbol")
	{
		const std::vector<std::string> tokens = {"c", "0.0", "0.0", "0.0"};
		REQUIRE_THROWS_AS(ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius),
		                  std::runtime_error);
	}

	SECTION("Number instead of element")
	{
		const std::vector<std::string> tokens = {"123", "0.0", "0.0", "0.0"};
		REQUIRE_THROWS_AS(ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius),
		                  std::runtime_error);
	}
}

TEST_CASE("ParseBasicCartesianCoordinateLine: Invalid coordinate values", "[Geometry][Parser][negative]")
{
	SECTION("Non-numeric coordinate")
	{
		const std::vector<std::string> tokens = {"H", "abc", "0.0", "0.0"};
		REQUIRE_THROWS_AS(ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius),
		                  std::exception);  // Could be various exception types
	}

	SECTION("Empty string as coordinate")
	{
		const std::vector<std::string> tokens = {"H", "", "0.0", "0.0"};
		REQUIRE_THROWS_AS(ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius),
		                  std::exception);
	}
}

TEST_CASE("ParseBasicCartesianCoordinateLine: Error message contains token count", "[Geometry][Parser][negative]")
{
	const std::vector<std::string> tokens = {"H", "1.0", "2.0"};  // Only 2 coordinates
	try
	{
		ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius);
		FAIL("Should have thrown exception");
	}
	catch (const std::runtime_error& e)
	{
		// Error message should mention the expected and actual token counts
		std::string message = e.what();
		CHECK(message.find('3') != std::string::npos);  // Expected 4 tokens (element + 3 coords)
		CHECK(message.find("Unexpected number of tokens") != std::string::npos);
	}
}

TEST_CASE("ParseBasicCartesianCoordinateLine: Case sensitivity", "[Geometry][Parser]")
{
	// Valid capitalization
	const std::vector<std::string> tokens1 = {"C", "0.0", "0.0", "0.0"};
	const auto result1 = ParseBasicCartesianCoordinateLine(tokens1, UnitOfMeasurement::Angstrom2BohrRadius);
	CHECK(result1.Element == Element::C);

	// Invalid capitalization
	const std::vector<std::string> tokens2 = {"c", "0.0", "0.0", "0.0"};
	REQUIRE_THROWS_AS(ParseBasicCartesianCoordinateLine(tokens2, UnitOfMeasurement::Angstrom2BohrRadius),
	                  std::runtime_error);
}

TEST_CASE("ParseBasicCartesianCoordinateLine: Special element symbols", "[Geometry][Parser]")
{
	SECTION("Two-letter symbols - first capital, second lowercase")
	{
		const std::vector<std::string> tokens = {"Fe", "0.0", "0.0", "0.0"};
		const auto result = ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius);
		CHECK(result.Element == Element::Fe);
	}

	SECTION("Two-letter symbol - wrong case")
	{
		const std::vector<std::string> tokens = {"FE", "0.0", "0.0", "0.0"};
		REQUIRE_THROWS_AS(ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius),
		                  std::runtime_error);
	}

	SECTION("Two-letter symbol - both lowercase")
	{
		const std::vector<std::string> tokens = {"fe", "0.0", "0.0", "0.0"};
		REQUIRE_THROWS_AS(ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius),
		                  std::runtime_error);
	}
}

TEST_CASE("ParseBasicCartesianCoordinateLine: Whitespace in tokens", "[Geometry][Parser]")
{
	// This test verifies that tokens with spaces are handled correctly
	// In practice, the tokenization should happen before calling this function
	const std::vector<std::string> tokens = {"H", " 1.0", " 2.0 ", " 3.0 "};

	// Parse<double> should handle leading/trailing whitespace
	const auto result = ParseBasicCartesianCoordinateLine(tokens, UnitOfMeasurement::Angstrom2BohrRadius);

	CHECK(result.Element == Element::H);
	constexpr double angstrom2bohr = 1.8897261254578281;
	CHECK(result.Position.x() == Approx(1.0 * angstrom2bohr));
	CHECK(result.Position.y() == Approx(2.0 * angstrom2bohr));
	CHECK(result.Position.z() == Approx(3.0 * angstrom2bohr));
}

//--------------------------------------------------------------------------------------------------------------------//
// ParseInternalCoordinateLine Tests
//--------------------------------------------------------------------------------------------------------------------//

static constexpr auto Bohr2Bohr = [](const double l) { return l; };

TEST_CASE("ParseInternalCoordinateLine: First atom at origin", "[geometry][z-matrix][parser]")
{
	const std::vector<Atom> atoms;
	const std::vector<std::string> tokens = {"C"};

	const auto result = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian);

	CHECK(result.Element == Element::C);
	CHECK(result.Position.x() == Approx(0.0));
	CHECK(result.Position.y() == Approx(0.0));
	CHECK(result.Position.z() == Approx(0.0));
}

TEST_CASE("ParseInternalCoordinateLine: Second atom with bond length", "[geometry][z-matrix][parser]")
{
	const std::vector atoms = {Atom(Element::H, Eigen::Vector3d{0, 0, 0})};
	const std::vector<std::string> tokens = {"O", "0", "1.5"};

	const auto result = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian);

	CHECK(result.Element == Element::O);
	CHECK(result.Position.x() == Approx(0.0));
	CHECK(result.Position.y() == Approx(0.0));
	CHECK(result.Position.z() == Approx(1.5));
}

TEST_CASE("ParseInternalCoordinateLine: Third atom with bond and angle", "[geometry][z-matrix][parser]")
{
	const std::vector atoms = {Atom(Element::H, Eigen::Vector3d{0, 0, 0}), Atom(Element::O, Eigen::Vector3d{0, 0, 1.5})};
	const std::vector<std::string> tokens = {"H", "0", "1.0", "1", "90"};

	const auto result = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian);

	CHECK(result.Element == Element::H);
	// For 90 deg angle, should be in XY plane
	CHECK(result.Position.x() == Approx(1.0).margin(1e-9));
	CHECK(result.Position.y() == Approx(0.0).margin(1e-9));
	CHECK(result.Position.z() == Approx(0.0).margin(1e-9));
}

TEST_CASE("ParseInternalCoordinateLine: Fourth atom with bond, angle, and dihedral", "[geometry][z-matrix][parser]")
{
	const std::vector atoms = {Atom(Element::H, Eigen::Vector3d{0, 0, 0}),
	                     Atom(Element::O, Eigen::Vector3d{0, 0, 1.5}),
	                     Atom(Element::H, Eigen::Vector3d{1.0, 0.0, 1.5})};
	const std::vector<std::string> tokens = {"C", "1", "1.2", "0", "90", "2", "180"};

	const auto result = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian);

	CHECK(result.Element == Element::C);
	// Distance to atom 1 should be close to 1.2
	const double distance = (result.Position - atoms[1].Position()).norm();
	CHECK(distance == Approx(1.2).margin(1e-9));
}

TEST_CASE("ParseInternalCoordinateLine: 1-based indexing convention", "[geometry][z-matrix][parser]")
{
	const std::vector atoms = {Atom(Element::H, Eigen::Vector3d{0, 0, 0}), Atom(Element::O, Eigen::Vector3d{0, 0, 1.5})};
	const std::vector<std::string> tokens = {"N", "1", "1.0", "2", "109.5"};

	const auto result = ParseInternalCoordinateLine<1>(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian);

	CHECK(result.Element == Element::N);
	// Check distance to first atom
	const double distance = (result.Position - atoms[0].Position()).norm();
	CHECK(distance == Approx(1.0).margin(1e-9));
}

TEST_CASE("ParseInternalCoordinateLine: Various elements", "[geometry][z-matrix][parser]")
{
	std::vector<Atom> atoms;

	SECTION("First atom - Carbon")
	{
		const std::vector<std::string> tokens = {"C"};
		const auto result = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian);
		CHECK(result.Element == Element::C);
	}

	SECTION("Second atom - Iron")
	{
		atoms.emplace_back(Element::H, Eigen::Vector3d{0, 0, 0});
		const std::vector<std::string> tokens = {"Fe", "0", "2.0"};
		const auto result = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian);
		CHECK(result.Element == Element::Fe);
	}

	SECTION("Third atom - Gold")
	{
		atoms.emplace_back(Element::H, Eigen::Vector3d{0, 0, 0});
		atoms.emplace_back(Element::O, Eigen::Vector3d{0, 0, 1.5});
		const std::vector<std::string> tokens = {"Au", "0", "1.5", "1", "120"};
		const auto result = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian);
		CHECK(result.Element == Element::Au);
	}
}

TEST_CASE("ParseInternalCoordinateLine: Unit conversion for bond length", "[geometry][z-matrix][parser]")
{
	const std::vector atoms = {Atom(Element::H, Eigen::Vector3d{0, 0, 0})};
	const std::vector<std::string> tokens = {"O", "0", "1.5"};

	SECTION("No conversion (Bohr)")
	{
		const auto result = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian);
		CHECK(result.Position.z() == Approx(1.5));
	}

	SECTION("Angstrom to Bohr conversion")
	{
		const auto result = ParseInternalCoordinateLine(
		        tokens, atoms, UnitOfMeasurement::Angstrom2BohrRadius, UnitOfMeasurement::Degree2Radian);
		constexpr double angstrom2bohr = 1.8897261254578281;
		CHECK(result.Position.z() == Approx(1.5 * angstrom2bohr));
	}

	SECTION("Custom conversion (2x)")
	{
		const auto result = ParseInternalCoordinateLine(
		        tokens, atoms, [](double x) { return 2 * x; }, UnitOfMeasurement::Degree2Radian);
		CHECK(result.Position.z() == Approx(3.0));
	}
}

TEST_CASE("ParseInternalCoordinateLine: Unit conversion for angle", "[geometry][z-matrix][parser]")
{
	std::vector atoms = {Atom(Element::H, Eigen::Vector3d{0, 0, 0}), Atom(Element::O, Eigen::Vector3d{0, 0, 1.5})};
	const std::vector<std::string> tokens = {"H", "0", "1.0", "1", "90"};

	SECTION("Degree to radian conversion")
	{
		const auto result = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian);
		// sin(90°) = 1, cos(90°) = 0
		CHECK(result.Position.x() == Approx(1).margin(1e-9));
		CHECK(result.Position.z() == Approx(0).margin(1e-9));
	}

	SECTION("No conversion (radians)")
	{
		const std::vector<std::string> radTokens = {"H", "0", "1.0", "1", "1.57079632679"};  // ~90° in radians
		const auto result = ParseInternalCoordinateLine(radTokens, atoms, Bohr2Bohr, [](double x) { return x; });
		CHECK(result.Position.x() == Approx(1).margin(1e-9));
		CHECK(result.Position.z() == Approx(0).margin(1e-9));
	}
}

TEST_CASE("ParseInternalCoordinateLine: Common bond angles", "[geometry][z-matrix][parser]")
{
	std::vector atoms = {Atom(Element::H, Eigen::Vector3d{0, 0, 0}), Atom(Element::O, Eigen::Vector3d{0, 0, 1.5})};

	SECTION("180 degrees - linear")
	{
		const std::vector<std::string> tokens = {"H", "0", "0.5", "1", "180"};
		const auto result = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian);
		CHECK(result.Position.z() == Approx(-0.5).margin(1e-9));
		CHECK(result.Position.x() == Approx(0.0).margin(1e-9));
	}

	SECTION("180 degrees - linear 2")
	{
		const std::vector<std::string> tokens = {"H", "1", "1.0", "0", "0"};
		const auto result = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian);
		CHECK(result.Position.z() == Approx(0.5).margin(1e-9));
		CHECK(result.Position.x() == Approx(0.0).margin(1e-9));
	}

	SECTION("120 degrees - trigonal planar")
	{
		const std::vector<std::string> tokens = {"H", "0", "1.0", "1", "120"};
		const auto result = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian);
		CHECK(result.Position.x() == Approx(std::cos(3.14159265358979323846 / 6)).margin(1e-9));
		CHECK(result.Position.z() == Approx(-std::sin(3.14159265358979323846 / 6)).margin(1e-9));
	}

	SECTION("109.5 degrees - tetrahedral")
	{
		const std::vector<std::string> tokens = {"H", "0", "1.0", "1", "109.5"};
		const auto result = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian);
		// Just check it's in a reasonable range
		CHECK(result.Position.x() > 0);
		CHECK(result.Position.z() < 0);
	}
}

TEST_CASE("ParseInternalCoordinateLine: Common dihedral angles", "[geometry][z-matrix][parser]")
{
	std::vector atoms = {Atom(Element::H, Eigen::Vector3d{0, 0, 0}),
	                     Atom(Element::O, Eigen::Vector3d{0, 0, 1.5}),
	                     Atom(Element::H, Eigen::Vector3d{1.0, 0.0, 1.5})};

	SECTION("0 degrees - cis configuration")
	{
		const std::vector<std::string> tokens = {"C", "1", "1.2", "0", "90", "2", "0"};
		const auto result = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian);
		// Should be in XZ plane
		CHECK(result.Position.y() == Approx(0.0).margin(1e-9));
	}

	SECTION("180 degrees - trans configuration")
	{
		const std::vector<std::string> tokens = {"C", "1", "1.2", "0", "90", "2", "180"};
		const auto result = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian);
		// Should also be in XZ plane
		CHECK(result.Position.y() == Approx(0.0).margin(1e-9));
	}

	SECTION("90 degrees - perpendicular")
	{
		const std::vector<std::string> tokens = {"C", "1", "1.2", "0", "90", "2", "90"};
		const auto result = ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian);
		// Y coordinate should be positive
		CHECK(result.Position.y() > 0);
	}
}

//--------------------------------------------------------------------------------------------------------------------//
// Negative tests for ParseInternalCoordinateLine
//--------------------------------------------------------------------------------------------------------------------//

TEST_CASE("ParseInternalCoordinateLine: Empty token list", "[geometry][z-matrix][parser][negative]")
{
	const std::vector<Atom> atoms;
	const std::vector<std::string> tokens = {};

	REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
	                  std::runtime_error);
}

TEST_CASE("ParseInternalCoordinateLine: Invalid bond length - negative", "[geometry][z-matrix][parser][negative]")
{
	const std::vector atoms = {Atom(Element::H, Eigen::Vector3d{0, 0, 0})};
	const std::vector<std::string> tokens = {"O", "0", "-1.5"};

	REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
	                  std::runtime_error);
}

TEST_CASE("ParseInternalCoordinateLine: Invalid bond length - zero", "[geometry][z-matrix][parser][negative]")
{
	const std::vector atoms = {Atom(Element::H, Eigen::Vector3d{0, 0, 0})};
	const std::vector<std::string> tokens = {"O", "0", "0.0"};

	REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
	                  std::runtime_error);
}

TEST_CASE("ParseInternalCoordinateLine: Out of range atom index - second atom", "[geometry][z-matrix][parser][negative]")
{
	const std::vector atoms = {Atom(Element::H, Eigen::Vector3d{0, 0, 0})};

	SECTION("Index too large")
	{
		const std::vector<std::string> tokens = {"O", "1", "1.5"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::out_of_range);
	}

	SECTION("Index negative with 0-based convention")
	{
		const std::vector<std::string> tokens = {"O", "-1", "1.5"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::exception);
	}
}

TEST_CASE("ParseInternalCoordinateLine: Out of range atom index - third atom", "[geometry][z-matrix][parser][negative]")
{
	const std::vector atoms = {Atom(Element::H, Eigen::Vector3d{0, 0, 0}), Atom(Element::O, Eigen::Vector3d{0, 0, 1.5})};

	SECTION("First index out of range")
	{
		const std::vector<std::string> tokens = {"H", "2", "1.0", "1", "90"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::out_of_range);
	}

	SECTION("Second index out of range")
	{
		const std::vector<std::string> tokens = {"H", "0", "1.0", "2", "90"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::out_of_range);
	}
}

TEST_CASE("ParseInternalCoordinateLine: Out of range atom index - fourth atom", "[geometry][z-matrix][parser][negative]")
{
	const std::vector atoms = {Atom(Element::H, Eigen::Vector3d{0, 0, 0}),
	                     Atom(Element::O, Eigen::Vector3d{0, 0, 1.5}),
	                     Atom(Element::H, Eigen::Vector3d{1.0, 0.0, 1.5})};

	SECTION("First index out of range")
	{
		const std::vector<std::string> tokens = {"C", "3", "1.2", "0", "90", "2", "180"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::out_of_range);
	}

	SECTION("Second index out of range")
	{
		const std::vector<std::string> tokens = {"C", "1", "1.2", "3", "90", "2", "180"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::out_of_range);
	}

	SECTION("Third index out of range")
	{
		const std::vector<std::string> tokens = {"C", "1", "1.2", "0", "90", "3", "180"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::out_of_range);
	}
}

TEST_CASE("ParseInternalCoordinateLine: Duplicate reference atoms - third atom",
          "[geometry][z-matrix][parser][negative]")
{
	const std::vector atoms = {Atom(Element::H, Eigen::Vector3d{0, 0, 0}), Atom(Element::O, Eigen::Vector3d{0, 0, 1.5})};
	const std::vector<std::string> tokens = {"H", "0", "1.0", "0", "90"};  // Both indices 0

	REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
	                  std::runtime_error);
}

TEST_CASE("ParseInternalCoordinateLine: Duplicate reference atoms - fourth atom",
          "[geometry][z-matrix][parser][negative]")
{
	const std::vector atoms = {Atom(Element::H, Eigen::Vector3d{0, 0, 0}),
	                     Atom(Element::O, Eigen::Vector3d{0, 0, 1.5}),
	                     Atom(Element::H, Eigen::Vector3d{1.0, 0.0, 1.5})};

	SECTION("First and second duplicate")
	{
		const std::vector<std::string> tokens = {"C", "1", "1.2", "1", "90", "2", "180"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::runtime_error);
	}

	SECTION("First and third duplicate")
	{
		const std::vector<std::string> tokens = {"C", "1", "1.2", "0", "90", "1", "180"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::runtime_error);
	}

	SECTION("Second and third duplicate")
	{
		const std::vector<std::string> tokens = {"C", "1", "1.2", "0", "90", "0", "180"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::runtime_error);
	}
}

TEST_CASE("ParseInternalCoordinateLine: Collinear reference atoms", "[geometry][z-matrix][parser][negative]")
{
	// All atoms on Z-axis - this will cause collinearity
	const std::vector atoms = {Atom(Element::H, Eigen::Vector3d{0, 0, 0}),
	                     Atom(Element::O, Eigen::Vector3d{0, 0, 1.5}),
	                     Atom(Element::H, Eigen::Vector3d{0, 0, 2.5})};
	const std::vector<std::string> tokens = {"C", "1", "1.2", "0", "90", "2", "180"};

	REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
	                  std::runtime_error);
}

TEST_CASE("ParseInternalCoordinateLine: Wrong token count - first atom", "[geometry][z-matrix][parser][negative]")
{
	SECTION("Too many tokens")
	{
		const std::vector<Atom> atoms;
		const std::vector<std::string> tokens = {"C", "0", "1.5"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::runtime_error);
	}
}

TEST_CASE("ParseInternalCoordinateLine: Wrong token count - second atom", "[geometry][z-matrix][parser][negative]")
{
	const std::vector atoms = {Atom(Element::H, Eigen::Vector3d{0, 0, 0})};

	SECTION("Only element")
	{
		const std::vector<std::string> tokens = {"O"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::runtime_error);
	}

	SECTION("Too many tokens")
	{
		const std::vector<std::string> tokens = {"O", "0", "1.5", "90"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::runtime_error);
	}
}

TEST_CASE("ParseInternalCoordinateLine: Wrong token count - third atom", "[geometry][z-matrix][parser][negative]")
{
	const std::vector atoms = {Atom(Element::H, Eigen::Vector3d{0, 0, 0}), Atom(Element::O, Eigen::Vector3d{0, 0, 1.5})};

	SECTION("Only element and bond")
	{
		const std::vector<std::string> tokens = {"H", "0", "1.0"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::runtime_error);
	}

	SECTION("Too many tokens")
	{
		const std::vector<std::string> tokens = {"H", "0", "1.0", "1", "90", "180"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::runtime_error);
	}
}

TEST_CASE("ParseInternalCoordinateLine: Wrong token count - fourth atom", "[geometry][z-matrix][parser][negative]")
{
	const std::vector atoms = {Atom(Element::H, Eigen::Vector3d{0, 0, 0}),
	                     Atom(Element::O, Eigen::Vector3d{0, 0, 1.5}),
	                     Atom(Element::H, Eigen::Vector3d{1.0, 0.0, 1.5})};

	SECTION("Missing dihedral")
	{
		const std::vector<std::string> tokens = {"C", "1", "1.2", "0", "90", "2"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::runtime_error);
	}

	SECTION("Too many tokens")
	{
		const std::vector<std::string> tokens = {"C", "1", "1.2", "0", "90", "2", "180", "extra"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::runtime_error);
	}
}

TEST_CASE("ParseInternalCoordinateLine: Invalid element symbol", "[geometry][z-matrix][parser][negative]")
{
	std::vector<Atom> atoms;

	SECTION("Invalid symbol for first atom")
	{
		const std::vector<std::string> tokens = {"Xx"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::runtime_error);
	}

	SECTION("Invalid symbol for second atom")
	{
		atoms.emplace_back(Element::H, Eigen::Vector3d{0, 0, 0});
		const std::vector<std::string> tokens = {"Yy", "0", "1.5"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::runtime_error);
	}

	SECTION("Lowercase element symbol")
	{
		atoms.emplace_back(Element::H, Eigen::Vector3d{0, 0, 0});
		const std::vector<std::string> tokens = {"c", "0", "1.5"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::runtime_error);
	}
}

TEST_CASE("ParseInternalCoordinateLine: Invalid numeric values", "[geometry][z-matrix][parser][negative]")
{
	std::vector atoms = {Atom(Element::H, Eigen::Vector3d{0, 0, 0})};

	SECTION("Non-numeric bond length")
	{
		const std::vector<std::string> tokens = {"O", "0", "abc"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::exception);
	}

	SECTION("Non-numeric atom index")
	{
		const std::vector<std::string> tokens = {"O", "x", "1.5"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::exception);
	}

	SECTION("Non-numeric angle")
	{
		atoms.emplace_back(Element::O, Eigen::Vector3d{0, 0, 1.5});
		const std::vector<std::string> tokens = {"H", "0", "1.0", "1", "abc"};
		REQUIRE_THROWS_AS(ParseInternalCoordinateLine(tokens, atoms, Bohr2Bohr, UnitOfMeasurement::Degree2Radian),
		                  std::exception);
	}
}

TEST_CASE("ParseInternalCoordinateLine: Token iterator overload", "[geometry][z-matrix][parser]")
{
	std::vector atoms = {Atom(Element::H, Eigen::Vector3d{0, 0, 0}), Atom(Element::O, Eigen::Vector3d{0, 0, 1.5})};
	const std::vector<std::string> tokens = {"N", "0", "1.0", "1", "120"};

	const auto result = ParseInternalCoordinateLine(
	        tokens.begin(), tokens.end(), atoms.begin(), atoms.end(), Bohr2Bohr, UnitOfMeasurement::Degree2Radian);

	CHECK(result.Element == Element::N);
	// Verify distance to first atom
	const double distance = (result.Position - atoms[0].Position()).norm();
	CHECK(distance == Approx(1.0).margin(1e-9));
}

TEST_CASE("ParseInternalCoordinateLine: Full Z-matrix - water molecule", "[geometry][z-matrix][parser][integration]")
{
	std::vector<std::string> zMatrixLines = {
	        "O",                // first atom
	        "H 0 0.96",         // second atom, bond length
	        "H 0 0.96 1 104.5"  // third atom, bond + angle
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

		const auto [element, pos] = ParseInternalCoordinateLine(
		        tokens, atoms, UnitOfMeasurement::Angstrom2BohrRadius, UnitOfMeasurement::Degree2Radian);
		atoms.emplace_back(element, pos);
	}

	REQUIRE(atoms.size() == 3);
	CHECK(atoms[0].Element() == Element::O);
	CHECK(atoms[1].Element() == Element::H);
	CHECK(atoms[2].Element() == Element::H);

	// Check O-H bond length
	constexpr double angstrom2bohr = 1.8897261254578281;
	CHECK(atoms[1].DistanceTo(atoms[0]) == Approx(0.96 * angstrom2bohr).margin(1e-9));
	CHECK(atoms[2].DistanceTo(atoms[0]) == Approx(0.96 * angstrom2bohr).margin(1e-9));

	// Check H-O-H angle
	const Eigen::Vector3d v1 = (atoms[1].Position() - atoms[0].Position()).normalized();
	const Eigen::Vector3d v2 = (atoms[2].Position() - atoms[0].Position()).normalized();
	const double angle = std::acos(v1.dot(v2));
	CHECK(angle == Approx(UnitOfMeasurement::Degree2Radian(104.5)).margin(1e-9));
}

TEST_CASE("ParseInternalCoordinateLine: Full Z-matrix - methane", "[geometry][z-matrix][parser][integration]")
{
	std::vector<std::string> zMatrixLines = {
	        "C",                      // first atom
	        "H 0 1.09",               // second atom
	        "H 0 1.09 1 109.5",       // third atom
	        "H 0 1.09 1 109.5 2 120"  // fourth atom
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

		const auto [element, pos] = ParseInternalCoordinateLine(
		        tokens, atoms, UnitOfMeasurement::Angstrom2BohrRadius, UnitOfMeasurement::Degree2Radian);
		atoms.emplace_back(element, pos);
	}

	REQUIRE(atoms.size() == 4);
	CHECK(atoms[0].Element() == Element::C);
	CHECK(atoms[1].Element() == Element::H);
	CHECK(atoms[2].Element() == Element::H);
	CHECK(atoms[3].Element() == Element::H);

	// Check all C-H bond lengths
	for (int i = 1; i < 4; ++i)
	{
		constexpr double angstrom2bohr = 1.8897261254578281;
		CHECK(atoms[i].DistanceTo(atoms[0]) == Approx(1.09 * angstrom2bohr).margin(1e-9));
	}
}
