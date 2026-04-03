// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Andy Brown

#pragma once

#include <Eigen/Dense>
#include <algorithm>
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-literal-operator"
#pragma	clang diagnostic ignored "-Wc++20-extensions"
#endif
#include <range/v3/range/conversion.hpp>
#if defined(__clang__)
#pragma clang diagnostic pop
#endif


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
