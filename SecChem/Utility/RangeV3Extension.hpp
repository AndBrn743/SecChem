// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Andy Brown

#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <range/v3/range/conversion.hpp>


template <>
struct ranges::detail::to_container::fn<ranges::detail::from_range<Eigen::VectorX>>
{
	template <typename Rng>
	auto operator()(const Rng& rng) const
	{
		eigen_assert(rng.size() >= 0 && rng.size() < std::numeric_limits<Eigen::Index>::max());
		Eigen::VectorX<std::decay_t<decltype(*rng.begin())>> res(static_cast<Eigen::Index>(rng.size()));
		std::transform(rng.begin(), rng.end(), res.begin(), [](const auto& val) { return val; });
		return res;
	}
};
