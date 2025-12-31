//
// Created by Andy on 12/31/2025.
//

#include "BasisSetExchangeJsonParserExample.hpp"
#include "../SecChem/Utility/Parser.hpp"

namespace Eigen
{
	template <typename Scalar>
	void from_json(const nlohmann::json& j, VectorX<Scalar>& out_vector)
	{
		out_vector.resize(j.size());
		for (Index i = 0; i < out_vector.size(); i++)
		{
			out_vector[i] = SecUtility::Parse<Scalar>(j[i]);
		}
	}
}  // namespace Eigen

using SecChem::AzimuthalQuantumNumber;
using SecChem::Element;
using SecChem::BasisSet::Gaussian::AngularMomentumBlock;
using SecChem::BasisSet::Gaussian::ContractedRadialOrbitalSet;

SecChem::BasisSet::Gaussian::BasisSet ParseBasisSetExchangeJson(const nlohmann::json& json)
{
	SecChem::BasisSet::Gaussian::BasisSet result;

	if (!json.contains("elements"))
	{
		throw std::runtime_error("Missing 'elements' key");
	}

	const auto& elements = json.at("elements");

	for (auto it = elements.begin(); it != elements.end(); ++it)
	{
		const Element element{SecUtility::Parse<int>(it.key())};

		const auto& elementData = it.value();

		if (!elementData.contains("electron_shells"))
		{
			std::clog << "Warning: An entry of element " << element
			          << " contains no 'electron_shells' key. This entry will be ignored." << std::endl;
			continue;
		}

		auto& angularMomentumBlocks = result.AddEntryFor(element);
		for (const auto& shell : elementData.at("electron_shells"))
		{
			const auto& angularMomenta = shell.at("angular_momentum");
			if (angularMomenta.size() != 1)
			{
				throw std::runtime_error("Only single angular momentum shells are supported. As a workaround, please "
				                         "consider split the input for each angular momentum");
			}
			AzimuthalQuantumNumber angularMomentum{angularMomenta[0].get<int>()};


			auto exponents = shell.at("exponents").get<Eigen::VectorXd>();


			const auto& coefficientsOfJson = shell.at("coefficients");
			const auto contractionSetCount = static_cast<Eigen::Index>(coefficientsOfJson.size());
			Eigen::MatrixXd contractionSets = Eigen::MatrixXd::Zero(exponents.size(), contractionSetCount);
			for (Eigen::Index c = 0; c < contractionSetCount; c++)
			{
				contractionSets.col(c) = coefficientsOfJson[c].get<Eigen::VectorXd>();
			}


			angularMomentumBlocks.emplace_back(
			        angularMomentum, ContractedRadialOrbitalSet{std::move(exponents), std::move(contractionSets)});
		}

		if (!elementData.contains("ecp_potentials"))
		{
			continue;
		}

		for (const auto& ecpJson : elementData.at("ecp_potentials"))
		{
			if (ecpJson.at("ecp_type").get<std::string>() != "scalar_ecp")
			{
				throw std::runtime_error("Only scalar ECPs are supported. Sorry ");
			}

			const auto& angularMomentaJson = ecpJson.at("angular_momentum");
			if (angularMomentaJson.size() != 1)
			{
				throw std::runtime_error("Only single angular momentum ECPs are supported. As a workaround, please "
				                         "consider split ECP blocks for each angular momentum");
			}
			const auto angularMomentum = static_cast<AzimuthalQuantumNumber>(angularMomentaJson[0].get<int>());
			const auto attachPointIterator = std::find_if(result[element].begin(),
			                                              result[element].end(),
			                                              [angularMomentum](const AngularMomentumBlock& amb)
			                                              { return amb.AngularMomentum() == angularMomentum; });
			if (attachPointIterator != result[element].end())
			{
				attachPointIterator->AddOrOverrideSemiLocalEcp(
				        {elementData.at("ecp_electrons").get<int>(),
				         ecpJson.at("coefficients").get<Eigen::VectorXd>(),
				         ecpJson.at("r_exponents").get<Eigen::VectorXd>(),
				         ecpJson.at("gaussian_exponents").get<Eigen::VectorXd>()});
			}
			else
			{
				result[element].emplace_back(angularMomentum,
				                             std::nullopt,
				                             SecChem::BasisSet::Gaussian::SemiLocalEcp{
				                                     elementData.at("ecp_electrons").get<int>(),
				                                     ecpJson.at("coefficients").get<Eigen::VectorXd>(),
				                                     ecpJson.at("r_exponents").get<Eigen::VectorXd>(),
				                                     ecpJson.at("gaussian_exponents").get<Eigen::VectorXd>()});
			}
		}
	}

	result.StandardizeRepresentation();
	return result;
}
