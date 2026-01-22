//
// Created by Andy on 12/26/2025.
//

#include <SecChem/BasisSet/Gaussian.hpp>
#include <SecChem/Utility/Parser.hpp>
#include <catch2/catch_test_macros.hpp>

#include "BasisSetExchangeJsonParserExample.hpp"

using namespace SecChem;
using namespace SecChem::BasisSet::Gaussian;

TEST_CASE("Sample parser should able to parse sample BSE JSON", "[BasisParser][BSE][JSON]")
{
	std::string sample = R"(
{
  "schema_version": "0.1",
  "basis_set_name": "test-mixed-basis",
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
        },
        {
          "ecp_type": "scalar_ecp",
          "angular_momentum": [1],
          "r_exponents": ["1"],
          "gaussian_exponents": ["4.500000E+00"],
          "coefficients": ["-5.000000E+00"]
        }
      ],
      "ecp_electrons": 28
    }
  }
}
)";

	const auto sampleBasisSet = ParseBasisSetExchangeJson(nlohmann::json::parse(sample));

	REQUIRE(sampleBasisSet.Has(Element::H));
	REQUIRE(sampleBasisSet.Has(Element::Ne));
	REQUIRE(sampleBasisSet[Element::H].IsInStandardRepresentation());
	REQUIRE(sampleBasisSet[Element::Ne].IsInStandardRepresentation());
	REQUIRE(sampleBasisSet.IsInStandardRepresentation());

	REQUIRE(sampleBasisSet[Element::H].AngularMomentumBlocks.size() == 2);
	REQUIRE(sampleBasisSet[Element::H].AngularMomentumBlocks[0].IsNotEmpty());
	REQUIRE(sampleBasisSet[Element::H].SemiLocalEcpProjectors.empty());
	REQUIRE(sampleBasisSet[Element::H].AngularMomentumBlocks[0].EqualsTo(AngularMomentumBlock{
	        AzimuthalQuantumNumber::S,
	        ContractedRadialOrbitalSet{Eigen::Vector3d{13.01, 1.962, 0.4446},
	                                   Eigen::Matrix<double, 3, 2>{{0.019685, 0}, {0.137977, 0}, {0.000000, 1}}}}));

	REQUIRE(sampleBasisSet[Element::H].AngularMomentumBlocks[1].IsNotEmpty());
	REQUIRE(sampleBasisSet[Element::H].AngularMomentumBlocks[1].EqualsTo(AngularMomentumBlock{
	        AzimuthalQuantumNumber::P,
	        ContractedRadialOrbitalSet{Eigen::Vector2d{0.727, 0.141},
	                                   Eigen::Matrix<double, 2, 2>{{0.430128, 0.5}, {0.678913, 0.5}}}}));

	REQUIRE(sampleBasisSet[Element::H].EcpElectronCount == 0);

	REQUIRE(sampleBasisSet[Element::Ne].AngularMomentumBlocks.size() == 1);
	REQUIRE(sampleBasisSet[Element::Ne].AngularMomentumBlocks[0].IsNotEmpty());
	REQUIRE(sampleBasisSet[Element::Ne].SemiLocalEcpProjectors.size() == 2);
	REQUIRE(sampleBasisSet[Element::Ne].AngularMomentumBlocks[0].EqualsTo(
	        AngularMomentumBlock{AzimuthalQuantumNumber::S,
	                             ContractedRadialOrbitalSet{Eigen::Vector2d{38.36, 5.77},
	                                                        Eigen::Matrix<double, 2, 1>{{0.023809}, {0.154891}}}}));

	const SemiLocalEcp ecp0{Eigen::Vector<double, 1>{-6}, Eigen::Vector<double, 1>{2}, Eigen::Vector<double, 1>{4}};
	REQUIRE(sampleBasisSet[Element::Ne].SemiLocalEcpProjectors[0].EqualsTo({AzimuthalQuantumNumber::S, ecp0}));
	const SemiLocalEcp ecp1{Eigen::Vector2d{-5, -4}, Eigen::Vector2d{1, 1.6}, Eigen::Vector2d{4.5, 3.5}};
	REQUIRE(sampleBasisSet[Element::Ne].SemiLocalEcpProjectors[1].EqualsTo({AzimuthalQuantumNumber::P, ecp1}));

	REQUIRE(sampleBasisSet[Element::Ne].EcpElectronCount == 28);
}