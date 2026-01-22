// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Andy Brown

#pragma once
#define SECCHEM_GAUSSIAN_BASIS_SET_INTERNAL
#include "Detail/Gaussian.Forward.hpp"
#undef SECCHEM_GAUSSIAN_BASIS_SET_INTERNAL

#include <Eigen/Dense>
#include <algorithm>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>

#include <SecChem/AzimuthalQuantumNumber.hpp>
#include <SecChem/Utility/IEquatableWithTolerance.hpp>

namespace SecChem::BasisSet::Gaussian
{
	class SemiLocalEcpTerm : public SecUtility::IEquatableWithTolerance<SemiLocalEcpTerm>
	{
		friend IEquatableWithTolerance;

	public:
		static constexpr Scalar ZeroTolerance =
		        SecUtility::Traits<SemiLocalEcpTerm>::DefaultEqualityComparisonTolerance;

		constexpr SemiLocalEcpTerm() noexcept = default;

		template <typename CoefficientSet, typename RExponentSet, typename GaussianExponentSet>
		SemiLocalEcpTerm(const CoefficientSet& coefficientSet,
		                 const RExponentSet& rExponentSet,
		                 const GaussianExponentSet& gaussianExponentSet)
		    : m_Data(CreateDataSet(coefficientSet, rExponentSet, gaussianExponentSet))
		{
			static_assert(std::is_base_of_v<Eigen::EigenBase<CoefficientSet>, CoefficientSet>);
			static_assert(std::is_base_of_v<Eigen::EigenBase<RExponentSet>, RExponentSet>);
			static_assert(std::is_base_of_v<Eigen::EigenBase<GaussianExponentSet>, GaussianExponentSet>);
		}

		bool IsEmpty() const noexcept
		{
			return m_Data.rows() == 0;
		}

		bool IsNotEmpty() const noexcept
		{
			return m_Data.rows() != 0;
		}

		Eigen::Index GaussianTermCount() const noexcept
		{
			return m_Data.rows();
		}

		auto Coefficients() const noexcept
		{
			return m_Data.col(0);
		}

		auto RExponents() const noexcept
		{
			return m_Data.col(1);
		}

		auto GaussianExponents() const noexcept
		{
			return m_Data.col(2);
		}

		const auto& Coefficient(const Eigen::Index index) const noexcept
		{
			return m_Data(index, 0);
		}

		const auto& RExponent(const Eigen::Index index) const noexcept
		{
			return m_Data(index, 1);
		}

		const auto& GaussianExponent(const Eigen::Index index) const noexcept
		{
			return m_Data(index, 2);
		}

		template <typename ForwardIterator, typename Getter>
		static std::enable_if_t<!std::is_same_v<std::decay_t<decltype(*std::declval<ForwardIterator>())>, void>,
		                        SemiLocalEcpTerm>
		Concat(const ForwardIterator begin, const ForwardIterator end, const Getter get)
		{
			if (begin == end)
			{
				throw std::runtime_error("SemiLocalEcp: input range is empty");
			}

			if (std::distance(begin, end) == 1)
			{
				return get(*begin);
			}

			Eigen::Matrix<Scalar, Eigen::Dynamic, 3> data(std::accumulate(begin,
			                                                              end,
			                                                              Eigen::Index{0},
			                                                              [get](const Eigen::Index acc, const auto& ecp)
			                                                              { return acc + get(ecp).m_Data.rows(); }),
			                                              3);

			auto it = begin;
			for (Eigen::Index offset = 0; it != end; ++it)
			{
				data.middleRows(offset, get(*it).m_Data.rows()) = get(*it).m_Data;
				offset += get(*it).m_Data.rows();
			}

			return SemiLocalEcpTerm{std::move(data)};
		}

		template <typename ForwardIterator>
		static std::enable_if_t<!std::is_same_v<std::decay_t<decltype(*std::declval<ForwardIterator>())>, void>,
		                        SemiLocalEcpTerm>
		Concat(const ForwardIterator begin, const ForwardIterator end)
		{
			return Concat(begin, end, [](auto&& item) -> decltype(auto) { return std::forward<decltype(item)>(item); });
		}

		template <typename... Sets>
		static std::enable_if_t<(std::is_same_v<std::decay_t<Sets>, SemiLocalEcpTerm> && ...), SemiLocalEcpTerm> Concat(
		        const Sets&... sets)
		{
			static_assert(sizeof...(Sets) > 0);
			const std::initializer_list<std::reference_wrapper<const SemiLocalEcpTerm>> list = {std::ref(sets)...};

			if constexpr (sizeof...(Sets) == 1)
			{
				return list.begin()->get();
			}
			else
			{
				return Concat(std::begin(list),
				              std::end(list),
				              [](const std::reference_wrapper<const SemiLocalEcpTerm> ref) -> const SemiLocalEcpTerm&
				              { return ref.get(); });
			}
		}

	private:
		explicit SemiLocalEcpTerm(Eigen::Matrix<Scalar, Eigen::Dynamic, 3> data) : m_Data(std::move(data))
		{
			/* NO CODE */
		}

		template <typename CoefficientSet, typename RExponentSet, typename GaussianExponentSet>
		static Eigen::Matrix<Scalar, Eigen::Dynamic, 3> CreateDataSet(const CoefficientSet& coefficientSet,
		                                                              const RExponentSet& rExponentSet,
		                                                              const GaussianExponentSet& gaussianExponentSet)
		{
			ValidateInput(coefficientSet, rExponentSet, gaussianExponentSet);
			Eigen::Matrix<Scalar, Eigen::Dynamic, 3> data(coefficientSet.size(), 3);
			data.col(0) = coefficientSet;
			data.col(1) = rExponentSet;
			data.col(2) = gaussianExponentSet;
			return data;
		}

		template <typename C, typename R, typename G>
		static void ValidateInput(const C& c, const R& r, const G& g)
		{
			if (c.size() != r.size() || c.size() != g.size())
			{
				throw std::invalid_argument("SemiLocalEcp: size mismatch");
			}

			if (!(c.rows() == 1 || c.cols() == 1))
			{
				throw std::invalid_argument("SemiLocalEcp: coefficients must be a vector");
			}

			if (!(r.rows() == 1 || r.cols() == 1))
			{
				throw std::invalid_argument("SemiLocalEcp: r exponents must be a vector");
			}

			if (!(g.rows() == 1 || g.cols() == 1))
			{
				throw std::invalid_argument("SemiLocalEcp: gaussian exponents must be a vector");
			}
		}

		bool EqualsTo_Impl(const SemiLocalEcpTerm& other, const Scalar tolerance = ZeroTolerance) const noexcept
		{
			return m_Data.rows() == other.m_Data.rows()
			       && (m_Data.size() == 0 || (m_Data - other.m_Data).cwiseAbs().maxCoeff() <= tolerance);
		}

		Eigen::Matrix<Scalar, Eigen::Dynamic, 3> m_Data;
	};


	class SemiLocalEcpChannel : public SecUtility::IEquatableWithTolerance<SemiLocalEcpChannel>
	{
		friend IEquatableWithTolerance;

	public:
		static constexpr Scalar ZeroTolerance =
		        SecUtility::Traits<SemiLocalEcpTerm>::DefaultEqualityComparisonTolerance;

		explicit SemiLocalEcpChannel(const AzimuthalQuantumNumber angularMomentum)
		    : m_AzimuthalQuantumNumber(angularMomentum)
		{
			/* NO CODE */
		}

		template <typename... Args>
		// ReSharper disable once CppNonExplicitConvertingConstructor
		SemiLocalEcpChannel(const AzimuthalQuantumNumber angularMomentum, Args&&... args)
		    : m_AzimuthalQuantumNumber(angularMomentum), m_SemiLocalEcpTerm(std::forward<Args>(args)...)
		{
			/* NO CODE */
		}

		AzimuthalQuantumNumber AngularMomentum() const noexcept
		{
			return m_AzimuthalQuantumNumber;
		}

		bool IsEmpty() const noexcept
		{
			return m_SemiLocalEcpTerm.IsEmpty();
		}

		bool IsNotEmpty() const noexcept
		{
			return m_SemiLocalEcpTerm.IsNotEmpty();
		}

		Eigen::Index GaussianTermCount() const noexcept
		{
			return m_SemiLocalEcpTerm.GaussianTermCount();
		}

		auto Coefficients() const noexcept
		{
			return m_SemiLocalEcpTerm.Coefficients();
		}

		auto RExponents() const noexcept
		{
			return m_SemiLocalEcpTerm.RExponents();
		}

		auto GaussianExponents() const noexcept
		{
			return m_SemiLocalEcpTerm.GaussianExponents();
		}

		const auto& Coefficient(const Eigen::Index index) const noexcept
		{
			return m_SemiLocalEcpTerm.Coefficient(index);
		}

		const auto& RExponent(const Eigen::Index index) const noexcept
		{
			return m_SemiLocalEcpTerm.RExponent(index);
		}

		const auto& GaussianExponent(const Eigen::Index index) const noexcept
		{
			return m_SemiLocalEcpTerm.GaussianExponent(index);
		}

		template <typename ForwardIterator, typename Getter>
		static std::enable_if_t<!std::is_same_v<std::decay_t<decltype(*std::declval<ForwardIterator>())>, void>,
		                        SemiLocalEcpChannel>
		Concat(const ForwardIterator begin, const ForwardIterator end, Getter&& get)
		{
			if (begin == end)
			{
				throw std::invalid_argument("Cannot concat empty range");
			}

			const auto l = begin->m_AzimuthalQuantumNumber;
			if (std::find_if(
			            begin, end, [&get, l](const auto& item) { return get(item).m_AzimuthalQuantumNumber != l; })
			    != end)
			{
				throw std::invalid_argument(
				        "Cannot concatenate SemiLocalEcpProjector objects of different angular momenta");
			}

			return SemiLocalEcpChannel{l,
			                           SemiLocalEcpTerm::Concat(begin,
			                                                    end,
			                                                    [&get](const auto& item) -> const SemiLocalEcpTerm&
			                                                    { return get(item).m_SemiLocalEcpTerm; })};
		}

		template <typename ForwardIterator>
		static std::enable_if_t<!std::is_same_v<std::decay_t<decltype(*std::declval<ForwardIterator>())>, void>,
		                        SemiLocalEcpChannel>
		Concat(const ForwardIterator begin, const ForwardIterator end)
		{
			return Concat(begin, end, [](const SemiLocalEcpChannel& p) -> const SemiLocalEcpChannel& { return p; });
		}

		template <typename Arg0, typename... Args>
		static std::enable_if_t<std::is_same_v<std::decay_t<Arg0>, SemiLocalEcpChannel>
		                                && (std::is_same_v<std::decay_t<Args>, SemiLocalEcpChannel> && ...),
		                        SemiLocalEcpChannel>
		Concat(Arg0&& arg0, Args&&... args)
		{
			const auto l = arg0.m_AzimuthalQuantumNumber;
			if (((args.m_AzimuthalQuantumNumber != l) || ...))
			{
				throw std::invalid_argument(
				        "Cannot concatenate SemiLocalEcpProjector objects of different angular momenta");
			}

			return SemiLocalEcpChannel{arg0.m_AzimuthalQuantumNumber,
			                           SemiLocalEcpTerm::Concat(arg0.m_SemiLocalEcpTerm, args.m_SemiLocalEcpTerm...)};
		}


	private:
		bool EqualsTo_Impl(const SemiLocalEcpChannel& other, const Scalar tolerance = ZeroTolerance) const noexcept
		{
			return m_AzimuthalQuantumNumber == other.m_AzimuthalQuantumNumber
			       && m_SemiLocalEcpTerm.EqualsTo(other.m_SemiLocalEcpTerm, tolerance);
		}


	private:
		AzimuthalQuantumNumber m_AzimuthalQuantumNumber;
		SemiLocalEcpTerm m_SemiLocalEcpTerm;
	};
}  // namespace SecChem::BasisSet::Gaussian
