//
// Created by Andy on 12/31/2025.
//

#pragma once

#include <nlohmann/json.hpp>
#include <SecChem/BasisSet/Gaussian.hpp>

SecChem::BasisSet::Gaussian::BasisSet ParseBasisSetExchangeJson(const nlohmann::json& json);
