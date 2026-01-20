// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Andy Brown

#pragma once
#if !defined(SECCHEM_GAUSSIAN_BASIS_SET_INTERNAL)
#error Do not include internal header files directly
#endif

namespace SecChem::BasisSet::Gaussian
{
	class SemiLocalEcp : public SecUtility::IEquatableWithTolerance<SemiLocalEcp>
	{
		friend IEquatableWithTolerance;

	public:
		static constexpr Scalar ZeroTolerance = 1e-15;

		constexpr SemiLocalEcp() noexcept = default;

		template <typename CoefficientSet, typename RExponentSet, typename GaussianExponentSet>
		SemiLocalEcp(const CoefficientSet& coefficientSet,
		             const RExponentSet& rExponentSet,
		             const GaussianExponentSet& gaussianExponentSet)
		    : m_Data(CreateDataSet(coefficientSet, rExponentSet, gaussianExponentSet))
		{
			static_assert(std::is_base_of_v<Eigen::EigenBase<CoefficientSet>, CoefficientSet>);
			static_assert(std::is_base_of_v<Eigen::EigenBase<RExponentSet>, RExponentSet>);
			static_assert(std::is_base_of_v<Eigen::EigenBase<GaussianExponentSet>, GaussianExponentSet>);
		}

		Eigen::Index TermCount() const noexcept
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
		static SemiLocalEcp Concat(const ForwardIterator begin, const ForwardIterator end, const Getter get)
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

			return SemiLocalEcp{std::move(data)};
		}

		template <typename ForwardIterator>
		static SemiLocalEcp Concat(const ForwardIterator begin, const ForwardIterator end)
		{
			return Concat(begin, end, [](auto&& item) -> decltype(auto) { return std::forward<decltype(item)>(item); });
		}


	private:
		explicit SemiLocalEcp(Eigen::Matrix<Scalar, Eigen::Dynamic, 3> data) : m_Data(std::move(data))
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

		bool EqualsTo_Impl(const SemiLocalEcp& other, const Scalar tolerance = ZeroTolerance) const noexcept
		{
			return m_Data.rows() == other.m_Data.rows()
			       && (m_Data.size() == 0 || (m_Data - other.m_Data).cwiseAbs().maxCoeff() <= tolerance);
		}

		Eigen::Matrix<Scalar, Eigen::Dynamic, 3> m_Data;
	};
}  // namespace SecChem::BasisSet::Gaussian
