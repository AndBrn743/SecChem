//
// Created by Andy on 12/31/2025.
//

// #include <catch2/catch_test_macros.hpp>

#include "BasisSetExchangeJsonParserExample.hpp"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <SecChem/Geometry.Input.hpp>

using namespace SecChem::Geometry::Input;

class MolecularInputInterpreter
{
	enum class GeometryFormat
	{
		Unknown,
		CartesianCoordinate,
		InternalCoordinate
	};

	struct AtomicAttributes
	{
		SecChem::AtomTag Tags = SecChem::AtomTag::None;
		double Mass = -1;
		double NuclearRadius = -1;
		std::string Basis{};
	};

	using LineParsingReturnType =
	        std::pair<SecChem::Atom, const std::vector<SecChem::BasisSet::Gaussian::AngularMomentumBlock>*>;

public:
	LineParsingReturnType ParseLine(
	        const std::string& line,
	        const std::vector<SecChem::Atom>& atoms,
	        const SecChem::BasisSet::Gaussian::SharedBasisSetLibrary& library,
	        const std::unordered_map<SecChem::Element,
	                                 const std::vector<SecChem::BasisSet::Gaussian::AngularMomentumBlock>*>&
	                defaultBasisSet)
	{
		const auto tokens = SecUtility::SplitRespectingQuotes(line);

		if (tokens.empty())
		{
			throw std::runtime_error("Not enough tokens in line \"" + line + '"');
		}

		const auto geomTokenCount =
		        std::distance(tokens.begin() + 1,
		                      std::find_if(tokens.begin() + 1,
		                                   tokens.end(),
		                                   [](const std::string& token)
		                                   {
			                                   return !token.empty()
			                                          && (token[0] == '_' || (token[0] >= 'A' && token[0] <= 'Z')
			                                              || (token[0] >= 'a' && token[0] <= 'z'));
		                                   }));

		if (atoms.empty())
		{
			if (geomTokenCount == 3)
			{
				m_GeometryFormat = GeometryFormat::CartesianCoordinate;
			}
			else if (geomTokenCount == 0)
			{
				m_GeometryFormat = GeometryFormat::InternalCoordinate;
			}
			else
			{
				throw std::runtime_error("Unexpected number of geometry tokens for the zeroth atom. In the cause of "
				                         "the zeroth atom, 3 geometry tokens was expected for Cartesian coordinate "
				                         "input, and no geometry token was expected for internal coordinate. Line \""
				                         + line + "\" have " + std::to_string(geomTokenCount) + " geometry token(s).");
			}
		}
		else if (m_GeometryFormat == GeometryFormat::Unknown)
		{
			throw std::runtime_error(std::string{__PRETTY_FUNCTION__} + ":" + std::to_string(__LINE__));
		}

		if (m_GeometryFormat == GeometryFormat::CartesianCoordinate && geomTokenCount != 3)
		{
			throw std::runtime_error("A Cartesian coordinate input line must contain exactly 3 geometry fields");
		}

		const auto [element, position] =
		        m_GeometryFormat == GeometryFormat::CartesianCoordinate
		                ? ParseBasicCartesianCoordinateLine(tokens.begin(),
		                                                    SecUtility::UnitOfMeasurement::Angstrom2BohrRadius)
		                : ParseInternalCoordinateLine<1>(tokens.begin(),
		                                                 tokens.begin() + 1 + geomTokenCount,
		                                                 atoms,
		                                                 SecUtility::UnitOfMeasurement::Angstrom2BohrRadius,
		                                                 SecUtility::UnitOfMeasurement::Degree2Radian);
		const auto [tags, mass, nuclearRadius, basis] = ParseAttributes(tokens.begin() + 4, tokens.end(), element);

		SecChem::Atom atom = SecChem::Atom::AtomWithMassAndNuclearRadius(element, position, mass, nuclearRadius, tags);


		if (basis.empty())
		{
			if (defaultBasisSet.find(element) == defaultBasisSet.end())
			{
				throw std::runtime_error("No default basis set configured for element " + element.ToString());
			}

			return {std::move(atom), defaultBasisSet.at(element)};
		}
		else
		{
			if (!library.Has(basis))
			{
				throw std::runtime_error("No basis set named \"" + basis + "\" inside the library");
			}

			if (!library[basis].Has(element))
			{
				throw std::runtime_error("Basis set named \"" + basis + "\" does not contain " + element.ToString());
			}

			return {std::move(atom), &library[basis][element]};
		}
	}


private:
	static SecChem::AtomTag ParseTags(const std::string& tagString)
	{
		using namespace SecChem;
		auto tags = AtomTag::None;
		std::stringstream ss(tagString);
		std::string token;

		while (std::getline(ss, token, ','))
		{
			if (token == "buffer")
			{
				tags |= AtomTag::Buffer;
			}
			else if (token == "link")
			{
				tags |= AtomTag::Link;
			}
			else if (token == "frozen")
			{
				tags |= AtomTag::Frozen;
			}
			else if (token == "gaussian_finite_nuclear")
			{
				tags |= AtomTag::GaussianFiniteNuclear;
			}
			else if (token == "thomas_fermi_finite_nuclear")
			{
				tags |= AtomTag::ThomasFermiFiniteNuclear;
			}
			else
			{
				throw std::invalid_argument("Unknown atomic tag: " + token);
			}
		}

		if (!IsValid(tags))
		{
			throw std::invalid_argument("Atomic tag \"" + tagString + "\" not valid");
		}

		return tags;
	}

	template <typename TokenIterator>
	static AtomicAttributes ParseAttributes(TokenIterator begin,
	                                        const TokenIterator end,
	                                        const SecChem::Element element)
	{
		AtomicAttributes attr{};
		auto& [tags, mass, nuclearRadius, basis] = attr;

		for (/* NO CODE */; begin != end; ++begin)
		{
			const auto& token = *begin;
			const auto eqPos = token.find('=');

			if (eqPos == std::string::npos)
			{
				throw std::invalid_argument("Invalid attribute: " + token);
			}

			const std::string key = token.substr(0, eqPos);
			std::string value = token.substr(eqPos + 1);

			// Strip quotes
			if (!value.empty() && value.front() == '"' && value.back() == '"')
			{
				value = value.substr(1, value.size() - 2);
			}

			if (key == "tag")
			{
				tags = ParseTags(value);
			}
			else if (key == "mass")
			{
				mass = std::stod(value);
			}
			else if (key == "nuclear_radius")
			{
				nuclearRadius = std::stod(value);
			}
			else if (key == "basis")
			{
				basis = value;
			}
			else
			{
				throw std::invalid_argument("Unknown atom attribute: " + key);
			}
		}

		if (mass == -1)
		{
			mass = element.Mass();
		}
		if (nuclearRadius == -1)
		{
			nuclearRadius = element.NuclearRadius();
		}

		return attr;
	}

	GeometryFormat m_GeometryFormat = GeometryFormat::Unknown;
};

static const SecChem::BasisSet::Gaussian::SharedBasisSetLibrary& Library();

TEST_CASE("BasisSetLibrary sanity", "[basis][library]")
{
	using namespace SecChem;
	const auto& lib = Library();

	REQUIRE_FALSE(lib.Has("3-21G"));

	REQUIRE(lib.Has("ANO-R0"));
	REQUIRE(lib["ANO-R0"].Has(Element::C));

	REQUIRE(lib.Has("def2-SVP"));
	REQUIRE(lib["def2-SVP"].Has(Element::C));
	REQUIRE(lib["def2-SVP"].Has(Element::H));
	REQUIRE(lib["def2-SVP"].Has(Element::O));

	REQUIRE(lib.Has("example-basis"));
	REQUIRE(lib["example-basis"].Has(Element::H));
	REQUIRE(lib["example-basis"].Has(Element::Ne));

	REQUIRE(lib["example-basis"][Element::Ne][0].HasSemiLocalEcp());
}

TEST_CASE("MolecularBasisSet build succeeds", "[basis][builder]")
{
	using namespace SecChem;

	const auto basis = SecChem::Builder<BasisSet::Gaussian::MolecularBasisSet>{Library()}
	                           .SetOrOverwriteGlobalDefaultBasisSetTo("def2-SVP")
	                           .SetOrOverwriteElementaryDefaultBasisSetTo("ANO-R0", Element::C)
	                           .BuildWith(MolecularInputInterpreter{},
	                                      R"(C 0 0 0 tag=frozen,gaussian_finite_nuclear
H 1.1 0 0
O 0 0.1e+1 0
Ne 1 1 1 basis="example-basis"
H 1.1 0 0
)");

	const auto& molecule = basis.Molecule();

	REQUIRE(molecule.AtomCount() == 5);

	REQUIRE(molecule[0].Element() == Element::C);
	REQUIRE(molecule[1].Element() == Element::H);
	REQUIRE(molecule[2].Element() == Element::O);
	REQUIRE(molecule[3].Element() == Element::Ne);
	REQUIRE(molecule[4].Element() == Element::H);

	REQUIRE(basis.UniqueElementaryBasisCount() == 4);
	REQUIRE(basis.UniqueElementaryBasis().size() == 4);

	REQUIRE(basis.AtomIndexFromPrimitiveSphericalOrbital(0) == 0);
	REQUIRE(basis.AtomIndexFromPrimitiveSphericalOrbital(1) == 0);
	REQUIRE(basis.AtomIndexFromPrimitiveSphericalOrbital(20) == 0);
	REQUIRE(basis.AtomIndexFromPrimitiveSphericalOrbital(30) == 0);
	REQUIRE(basis.AtomIndexFromPrimitiveSphericalOrbital(32) == 0);
	REQUIRE(basis.AtomIndexFromPrimitiveSphericalOrbital(33) == 1);

	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(0) == 0);
	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(1) == 0);
	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(2) == 0);
	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(3) == 0);
	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(4) == 0);

	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(5) == 1);
	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(6) == 1);
	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(7) == 1);
	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(8) == 1);
	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(9) == 1);

	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(10) == 2);  // 1s
	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(11) == 2);  // 2s
	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(12) == 2);  // 3s
	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(13) == 2);  // 2p(-1)
	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(14) == 2);  // 2p(0)
	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(15) == 2);  // 2p(+1)
	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(16) == 2);  // 3p(-1)
	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(17) == 2);  // 3p(0)
	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(18) == 2);  // 3p(+1)
	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(19) == 2);  // 3d(-2)
	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(20) == 2);  // 3d(-1)
	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(21) == 2);  // 3d(0)
	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(22) == 2);  // 3d(+1)
	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(23) == 2);  // 3d(+2)

	REQUIRE(basis.AtomIndexFromContractedSphericalOrbital(24) == 3);

	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(0) == std::pair{Eigen::Index{0}, 1_Sharp});
	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(1) == std::pair{Eigen::Index{0}, 2_Sharp});
	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(2) == std::pair{Eigen::Index{0}, 2_Principal});
	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(3) == std::pair{Eigen::Index{0}, 2_Principal});
	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(4) == std::pair{Eigen::Index{0}, 2_Principal});

	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(5) == std::pair{Eigen::Index{1}, 1_Sharp});
	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(6) == std::pair{Eigen::Index{1}, 2_Sharp});
	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(7) == std::pair{Eigen::Index{1}, 2_Principal});
	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(8) == std::pair{Eigen::Index{1}, 2_Principal});
	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(9) == std::pair{Eigen::Index{1}, 2_Principal});

	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(10) == std::pair{Eigen::Index{2}, 1_Sharp});
	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(11) == std::pair{Eigen::Index{2}, 2_Sharp});
	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(12) == std::pair{Eigen::Index{2}, 3_Sharp});
	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(13) == std::pair{Eigen::Index{2}, 2_Principal});
	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(14) == std::pair{Eigen::Index{2}, 2_Principal});
	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(15) == std::pair{Eigen::Index{2}, 2_Principal});
	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(16) == std::pair{Eigen::Index{2}, 3_Principal});
	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(17) == std::pair{Eigen::Index{2}, 3_Principal});
	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(18) == std::pair{Eigen::Index{2}, 3_Principal});
	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(19) == std::pair{Eigen::Index{2}, 3_Diffuse});
	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(20) == std::pair{Eigen::Index{2}, 3_Diffuse});
	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(21) == std::pair{Eigen::Index{2}, 3_Diffuse});
	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(22) == std::pair{Eigen::Index{2}, 3_Diffuse});
	REQUIRE(basis.AtomIndexAndSubShellFromContractedSphericalOrbital(23) == std::pair{Eigen::Index{2}, 3_Diffuse});

	REQUIRE(basis.ContractedSubShellCount() == 3 + 3 + 6 + 1 + 3);
	REQUIRE(basis.PrimitiveSubShellCount() == 19 + 5 + 12 + 2 + 5);
}

TEST_CASE("Atom tags are preserved", "[basis][parser]")
{
	using namespace SecChem;

	const auto mbs = SecChem::Builder<BasisSet::Gaussian::MolecularBasisSet>{Library()}
	                         .SetOrOverwriteGlobalDefaultBasisSetTo("def2-SVP")
	                         .SetOrOverwriteElementaryDefaultBasisSetTo("ANO-R0", Element::C)
	                         .BuildWith(MolecularInputInterpreter{},
	                                    R"(C 0 0 0 tag=frozen,gaussian_finite_nuclear
H 1.1 0 0
O 0 1.0 0
Ne 1 1 1 basis="example-basis")");

	const auto& atomC = mbs.Molecule()[0];

	REQUIRE(atomC.IsFrozen());
	REQUIRE(atomC.IsWithGaussianFiniteNuclear());
}

TEST_CASE("Basis assignment precedence is respected", "[basis][assignment]")
{
	using namespace SecChem;
	const auto mbs = SecChem::Builder<BasisSet::Gaussian::MolecularBasisSet>{Library()}
	                         .SetOrOverwriteGlobalDefaultBasisSetTo("def2-SVP")
	                         .SetOrOverwriteElementaryDefaultBasisSetTo("ANO-R0", Element::C)
	                         .BuildWith(MolecularInputInterpreter{},
	                                    R"(C 0 0 0
H 1.1 0 0
O 0 1.0 0
Ne 1 1 1 basis="example-basis")");

	const auto& mol = mbs.Molecule();

	// C → ANO-R0 (element override)
	REQUIRE(&mbs.ElementaryBasisAt(0) == &Library()["ANO-R0"][Element::C]);
	REQUIRE(&mbs.ElementaryBasisOf(mol[0]) == &Library()["ANO-R0"][Element::C]);

	// H → def2-SVP (global default)
	REQUIRE(&mbs.ElementaryBasisAt(1) == &Library()["def2-SVP"][Element::H]);
	REQUIRE(&mbs.ElementaryBasisOf(mol[1]) == &Library()["def2-SVP"][Element::H]);

	// O → def2-SVP (global default)
	REQUIRE(&mbs.ElementaryBasisAt(2) == &Library()["def2-SVP"][Element::O]);
	REQUIRE(&mbs.ElementaryBasisOf(mol[2]) == &Library()["def2-SVP"][Element::O]);

	// Ne → example-basis (per-atom override)
	REQUIRE(&mbs.ElementaryBasisAt(3) == &Library()["example-basis"][Element::Ne]);
	REQUIRE(&mbs.ElementaryBasisOf(mol[3]) == &Library()["example-basis"][Element::Ne]);

	REQUIRE_THROWS_AS(mbs.ElementaryBasisAt(4), std::out_of_range);
	REQUIRE_THROWS_AS(mbs.ElementaryBasisAt(-1), std::out_of_range);
}

TEST_CASE("ECP-only angular momentum blocks are preserved", "[basis][ecp]")
{
	const auto mbs = SecChem::Builder<SecChem::BasisSet::Gaussian::MolecularBasisSet>{Library()}
	                         .SetOrOverwriteGlobalDefaultBasisSetTo("def2-SVP")
	                         .BuildWith(MolecularInputInterpreter{}, R"(Ne 0 0 0 basis="example-basis")");

	REQUIRE(mbs.ElementaryBasisAt(0).size() == 2);
	REQUIRE(mbs.ElementaryBasisAt(0)[0].HasSemiLocalEcp());
	REQUIRE(mbs.ElementaryBasisAt(0)[1].HasSemiLocalEcp());
	REQUIRE(mbs.ElementaryBasisAt(0)[0].HasOrbital());
	REQUIRE_FALSE(mbs.ElementaryBasisAt(0)[1].HasOrbital());
}

TEST_CASE("Elementary default basis cannot be set before global default", "[basis][builder][negative]")
{
	auto builder = SecChem::Builder<SecChem::BasisSet::Gaussian::MolecularBasisSet>{Library()};
	builder.SetOrOverwriteGlobalDefaultBasisSetTo("def2-SVP");

	REQUIRE_THROWS_AS(builder.SetOrOverwriteElementaryDefaultBasisSetTo("ANO-R0", SecChem::Element::Na),
	                  std::runtime_error);
}

TEST_CASE("Setting non-existent global default basis throws", "[basis][builder][negative]")
{
	auto builder = SecChem::Builder<SecChem::BasisSet::Gaussian::MolecularBasisSet>{Library()};

	REQUIRE_THROWS_AS(builder.SetOrOverwriteGlobalDefaultBasisSetTo("3-21G"), std::runtime_error);
}

TEST_CASE("Setting non-existent elementary default basis throws", "[basis][builder][negative]")
{
	auto builder = SecChem::Builder<SecChem::BasisSet::Gaussian::MolecularBasisSet>{Library()};
	builder.SetOrOverwriteGlobalDefaultBasisSetTo("def2-SVP");

	REQUIRE_THROWS_AS(builder.SetOrOverwriteElementaryDefaultBasisSetTo("3-21G", SecChem::Element::C),
	                  std::runtime_error);
}

TEST_CASE("Elementary default basis lacking element throws", "[basis][builder][negative]")
{
	auto builder = SecChem::Builder<SecChem::BasisSet::Gaussian::MolecularBasisSet>{Library()};

	builder.SetOrOverwriteGlobalDefaultBasisSetTo("def2-SVP");

	// example-basis does NOT contain carbon
	REQUIRE_THROWS_AS(builder.SetOrOverwriteElementaryDefaultBasisSetTo("example-basis", SecChem::Element::C),
	                  std::runtime_error);
}

TEST_CASE("Per-atom basis override with unknown basis throws", "[basis][parser][negative]")
{
	auto builder = SecChem::Builder<SecChem::BasisSet::Gaussian::MolecularBasisSet>{Library()}
	                       .SetOrOverwriteGlobalDefaultBasisSetTo("def2-SVP");

	REQUIRE_THROWS_AS(builder.BuildWith(MolecularInputInterpreter{}, R"(H 0 0 0 basis="3-21G")"), std::runtime_error);
}

TEST_CASE("Per-atom basis override lacking element throws", "[basis][parser][negative]")
{
	auto builder = SecChem::Builder<SecChem::BasisSet::Gaussian::MolecularBasisSet>{Library()}
	                       .SetOrOverwriteGlobalDefaultBasisSetTo("def2-SVP");

	// example-basis does not contain oxygen
	REQUIRE_THROWS_AS(builder.BuildWith(MolecularInputInterpreter{}, R"(O 0 0 0 basis="example-basis")"),
	                  std::runtime_error);
}

TEST_CASE("Missing basis assignment throws when no default is set", "[basis][negative]")
{
	auto builder = SecChem::Builder<SecChem::BasisSet::Gaussian::MolecularBasisSet>{Library()};

	REQUIRE_THROWS_AS(builder.BuildWith(MolecularInputInterpreter{}, R"(H 0 0 0)"), std::runtime_error);
}

TEST_CASE("Invalid atom line format throws", "[parser][negative]")
{
	auto builder = SecChem::Builder<SecChem::BasisSet::Gaussian::MolecularBasisSet>{Library()}
	                       .SetOrOverwriteGlobalDefaultBasisSetTo("def2-SVP");

	REQUIRE_THROWS_AS(builder.BuildWith(MolecularInputInterpreter{},
	                                    R"(C 0 0)"),  // missing coordinate
	                  std::runtime_error);
}

TEST_CASE("Unknown element symbol throws", "[parser][negative]")
{
	auto builder = SecChem::Builder<SecChem::BasisSet::Gaussian::MolecularBasisSet>{Library()}
	                       .SetOrOverwriteGlobalDefaultBasisSetTo("def2-SVP");

	REQUIRE_THROWS_AS(builder.BuildWith(MolecularInputInterpreter{}, R"(Xx 0 0 0)"), std::runtime_error);
}

//--------------------------------------------------------------------------------------------------------------------//

static const SecChem::BasisSet::Gaussian::SharedBasisSetLibrary& Library()
{
	using namespace SecChem;
	const auto example0 = nlohmann::json::parse(R"(
{
    "molssi_bse_schema": {
        "schema_type": "complete",
        "schema_version": "0.1"
    },
    "revision_description": "Data from supporting information of erratum",
    "revision_date": "2021-09-25",
    "elements": {
        "6": {
            "electron_shells": [
                {
                    "function_type": "gto",
                    "region": "",
                    "angular_momentum": [
                        0
                    ],
                    "exponents": [
                        "15566.69000000",
                        "2405.49600000",
                        "525.67290000",
                        "146.49690000",
                        "47.97475000",
                        "17.33014000",
                        "6.58396100",
                        "2.53985700",
                        "0.96973270",
                        "0.35953110",
                        "0.12744200",
                        "0.04517401"
                    ],
                    "coefficients": [
                        [
                            "0.0003189872",
                            "0.0019590813",
                            "0.0102385738",
                            "0.0405238695",
                            "0.1258839645",
                            "0.2933115901",
                            "0.4369257167",
                            "0.2362314237",
                            "0.0079326342",
                            "0.0034392481",
                            "-0.0020754636",
                            "0.0005419770"
                        ],
                        [
                            "-0.0000707180",
                            "-0.0004342985",
                            "-0.0022847848",
                            "-0.0091545899",
                            "-0.0297361058",
                            "-0.0760753928",
                            "-0.1509354847",
                            "-0.1515213903",
                            "0.1868849319",
                            "0.6646021627",
                            "0.3158364657",
                            "-0.0073833364"
                        ]
                    ]
                },
                {
                    "function_type": "gto",
                    "region": "",
                    "angular_momentum": [
                        1
                    ],
                    "exponents": [
                        "34.79377000",
                        "7.96239500",
                        "2.37594400",
                        "0.81293870",
                        "0.28812230",
                        "0.10029750",
                        "0.03491430"
                    ],
                    "coefficients": [
                        [
                            "0.0055910270",
                            "0.0375157776",
                            "0.1480134382",
                            "0.3688832987",
                            "0.4813700216",
                            "0.1940050552",
                            "-0.0141156477"
                        ]
                    ]
                }
            ],
            "references": [
                {
                    "reference_description": "ANO-R0",
                    "reference_keys": [
                        "zobel2020a",
                        "zobel2021a"
                    ]
                }
            ]
        }
    },
    "version": "2",
    "function_types": [
        "gto"
    ],
    "names": [
        "ANO-R0"
    ],
    "tags": [],
    "family": "ano",
    "description": "ANO-R0",
    "role": "orbital",
    "auxiliaries": {},
    "name": "ANO-R0"
}
)");

	const auto example1 = nlohmann::json::parse(R"(
{
    "molssi_bse_schema": {
        "schema_type": "complete",
        "schema_version": "0.1"
    },
    "revision_description": "Data from Turbomole 7.3",
    "revision_date": "2018-09-07",
    "elements": {
        "1": {
            "electron_shells": [
                {
                    "function_type": "gto",
                    "region": "",
                    "angular_momentum": [
                        0
                    ],
                    "exponents": [
                        "13.0107010",
                        "1.9622572",
                        "0.44453796"
                    ],
                    "coefficients": [
                        [
                            "0.19682158E-01",
                            "0.13796524",
                            "0.47831935"
                        ]
                    ]
                },
                {
                    "function_type": "gto",
                    "region": "",
                    "angular_momentum": [
                        0
                    ],
                    "exponents": [
                        "0.12194962"
                    ],
                    "coefficients": [
                        [
                            "1.0000000"
                        ]
                    ]
                },
                {
                    "function_type": "gto",
                    "region": "",
                    "angular_momentum": [
                        1
                    ],
                    "exponents": [
                        "0.8000000"
                    ],
                    "coefficients": [
                        [
                            "1.0000000"
                        ]
                    ]
                }
            ],
            "references": [
                {
                    "reference_description": "Common base shells for def2-SV* basis sets",
                    "reference_keys": [
                        "weigend2005a"
                    ]
                },
                {
                    "reference_description": "Additional functions to form def2-SVP",
                    "reference_keys": [
                        "weigend2005a"
                    ]
                }
            ]
        },
        "6": {
            "electron_shells": [
                {
                    "function_type": "gto",
                    "region": "",
                    "angular_momentum": [
                        0
                    ],
                    "exponents": [
                        "1238.4016938",
                        "186.29004992",
                        "42.251176346",
                        "11.676557932",
                        "3.5930506482"
                    ],
                    "coefficients": [
                        [
                            "0.54568832082E-02",
                            "0.40638409211E-01",
                            "0.18025593888",
                            "0.46315121755",
                            "0.44087173314"
                        ]
                    ]
                },
                {
                    "function_type": "gto",
                    "region": "",
                    "angular_momentum": [
                        0
                    ],
                    "exponents": [
                        "0.40245147363"
                    ],
                    "coefficients": [
                        [
                            "1.0000000"
                        ]
                    ]
                },
                {
                    "function_type": "gto",
                    "region": "",
                    "angular_momentum": [
                        0
                    ],
                    "exponents": [
                        "0.13090182668"
                    ],
                    "coefficients": [
                        [
                            "1.0000000"
                        ]
                    ]
                },
                {
                    "function_type": "gto",
                    "region": "",
                    "angular_momentum": [
                        1
                    ],
                    "exponents": [
                        "9.4680970621",
                        "2.0103545142",
                        "0.54771004707"
                    ],
                    "coefficients": [
                        [
                            "0.38387871728E-01",
                            "0.21117025112",
                            "0.51328172114"
                        ]
                    ]
                },
                {
                    "function_type": "gto",
                    "region": "",
                    "angular_momentum": [
                        1
                    ],
                    "exponents": [
                        "0.15268613795"
                    ],
                    "coefficients": [
                        [
                            "1.0000000"
                        ]
                    ]
                },
                {
                    "function_type": "gto_spherical",
                    "region": "",
                    "angular_momentum": [
                        2
                    ],
                    "exponents": [
                        "0.8000000"
                    ],
                    "coefficients": [
                        [
                            "1.0000000"
                        ]
                    ]
                }
            ],
            "references": [
                {
                    "reference_description": "Common base shells for def2-SV* basis sets",
                    "reference_keys": [
                        "weigend2005a"
                    ]
                }
            ]
        },
        "8": {
            "electron_shells": [
                {
                    "function_type": "gto",
                    "region": "",
                    "angular_momentum": [
                        0
                    ],
                    "exponents": [
                        "2266.1767785",
                        "340.87010191",
                        "77.363135167",
                        "21.479644940",
                        "6.6589433124"
                    ],
                    "coefficients": [
                        [
                            "-0.53431809926E-02",
                            "-0.39890039230E-01",
                            "-0.17853911985",
                            "-0.46427684959",
                            "-0.44309745172"
                        ]
                    ]
                },
                {
                    "function_type": "gto",
                    "region": "",
                    "angular_momentum": [
                        0
                    ],
                    "exponents": [
                        "0.80975975668"
                    ],
                    "coefficients": [
                        [
                            "1.0000000"
                        ]
                    ]
                },
                {
                    "function_type": "gto",
                    "region": "",
                    "angular_momentum": [
                        0
                    ],
                    "exponents": [
                        "0.25530772234"
                    ],
                    "coefficients": [
                        [
                            "1.0000000"
                        ]
                    ]
                },
                {
                    "function_type": "gto",
                    "region": "",
                    "angular_momentum": [
                        1
                    ],
                    "exponents": [
                        "17.721504317",
                        "3.8635505440",
                        "1.0480920883"
                    ],
                    "coefficients": [
                        [
                            "0.43394573193E-01",
                            "0.23094120765",
                            "0.51375311064"
                        ]
                    ]
                },
                {
                    "function_type": "gto",
                    "region": "",
                    "angular_momentum": [
                        1
                    ],
                    "exponents": [
                        "0.27641544411"
                    ],
                    "coefficients": [
                        [
                            "1.0000000"
                        ]
                    ]
                },
                {
                    "function_type": "gto_spherical",
                    "region": "",
                    "angular_momentum": [
                        2
                    ],
                    "exponents": [
                        "1.2000000"
                    ],
                    "coefficients": [
                        [
                            "1.0000000"
                        ]
                    ]
                }
            ],
            "references": [
                {
                    "reference_description": "Common base shells for def2-SV* basis sets",
                    "reference_keys": [
                        "weigend2005a"
                    ]
                }
            ]
        }
    },
    "version": "1",
    "function_types": [
        "gto",
        "gto_spherical"
    ],
    "names": [
        "def2-SVP"
    ],
    "tags": [],
    "family": "ahlrichs",
    "description": "def2-SVP",
    "role": "orbital",
    "auxiliaries": {
        "jkfit": "def2-universal-jkfit",
        "jfit": "def2-universal-jfit",
        "rifit": "def2-svp-rifit"
    },
    "name": "def2-SVP"
}
)");

	const auto example2 = nlohmann::json::parse(R"(
{
  "schema_version": "0.1",
  "name": "example-basis",
  "elements": {
    "1": {
      "electron_shells": [
        {
          "angular_momentum": [0],
          "exponents": [
            "13.010000E+00",
            "1.962000E+00"
          ],
          "coefficients": [
            ["0.019685E+00", "0.137977E+00"]
          ]
        },
        {
          "angular_momentum": [0],
          "exponents": [
            "0.444600E+00"
          ],
          "coefficients": [
            ["1.000000E+00"]
          ]
        },
        {
          "angular_momentum": [1],
          "exponents": [
            "0.727000E+00",
            "0.141000E+00"
          ],
          "coefficients": [
            ["0.430128E+00", "0.678913E+00"],
            ["0.500000E+00", "0.500000E+00"]
          ]
        }
      ]
    },

    "10": {
      "electron_shells": [
        {
          "angular_momentum": [0],
          "exponents": [
            "38.360000E+00",
            "5.770000E+00"
          ],
          "coefficients": [
            ["0.023809E+00", "0.154891E+00"]
          ]
        }
      ],

      "ecp_potentials": [
        {
          "ecp_type": "scalar_ecp",
          "angular_momentum": [0],
          "r_exponents": ["2"],
          "gaussian_exponents": ["4.000000E+00"],
          "coefficients": ["-6.000000E+00"]
        },
        {
          "ecp_type": "scalar_ecp",
          "angular_momentum": [1],
          "r_exponents": ["1.6"],
          "gaussian_exponents": ["3.500000E+00"],
          "coefficients": ["-4.000000E+00"]
        }
      ],
      "ecp_electrons": 28
    }
  }
}
)");

	static auto library = Builder<BasisSet::Gaussian::SharedBasisSetLibrary>{}
	                              .AddBasisSet(example0.at("name").get<std::string>(),
	                                           BasisSet::Gaussian::SharedBasisSet{ParseBasisSetExchangeJson(example0)})
	                              .AddBasisSet(example1.at("name").get<std::string>(),
	                                           BasisSet::Gaussian::SharedBasisSet{ParseBasisSetExchangeJson(example1)})
	                              .AddBasisSet(example2.at("name").get<std::string>(),
	                                           BasisSet::Gaussian::SharedBasisSet{ParseBasisSetExchangeJson(example2)})
	                              .Build();

	return library;
}
