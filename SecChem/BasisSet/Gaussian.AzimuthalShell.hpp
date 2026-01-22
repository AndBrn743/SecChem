// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Andy Brown

#pragma once
#define SECCHEM_GAUSSIAN_BASIS_SET_INTERNAL
#include "Detail/Gaussian.Forward.hpp"
#undef SECCHEM_GAUSSIAN_BASIS_SET_INTERNAL

#include <Eigen/Dense>
#include <algorithm>
#include <numeric>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>
#include <type_traits>
#include <utility>
#include <vector>

#include <SecChem/AzimuthalQuantumNumber.hpp>
#include <SecChem/Element.hpp>
#include <SecChem/System.hpp>
#include <SecChem/Utility/IEquatableWithTolerance.hpp>

#include "Gaussian.ContractedRadialOrbitalSet.hpp"

namespace SecChem::BasisSet::Gaussian
{
	template <typename Derived>
	class AbstractAzimuthalShell : public SecUtility::IEquatableWithTolerance<AbstractAzimuthalShell<Derived>>
	{
		friend SecUtility::IEquatableWithTolerance<AbstractAzimuthalShell>;
		friend Derived;

	public:
		const auto& ExponentSet() const noexcept
		{
			return static_cast<const Derived&>(*this).ExponentSet_Impl();
		}

		const auto& ContractionSets() const noexcept
		{
			return static_cast<const Derived&>(*this).ContractionSets_Impl();
		}

		Eigen::Index PrimitiveShellCount() const noexcept
		{
			return static_cast<const Derived&>(*this).PrimitiveShellCount_Impl();
		}

		auto PrimitiveShells(const int principalQuantumNumberOffset = 0) const noexcept
		{
			const auto l = AngularMomentum();
			const auto n0 = l.MinPrincipalQuantumNumber() + principalQuantumNumberOffset;

			return ranges::views::iota(0 + n0, static_cast<int>(PrimitiveShellCount()) + n0)
			       | ranges::views::transform([l](const int n) { return ElectronicSubshell{n, l}; });
		}

		Eigen::Index ContractedShellCount() const noexcept
		{
			return static_cast<const Derived&>(*this).ContractedShellCount_Impl();
		}

		auto ContractedShells(const int principalQuantumNumberOffset = 0) const noexcept
		{
			const auto l = AngularMomentum();
			const auto n0 = l.MinPrincipalQuantumNumber() + principalQuantumNumberOffset;

			return ranges::views::iota(0 + n0, static_cast<int>(ContractedShellCount()) + n0)
			       | ranges::views::transform([l](const int n) { return ElectronicSubshell{n, l}; });
		}

		Eigen::Index PrimitiveCartesianOrbitalCount() const noexcept
		{
			return PrimitiveShellCount() * m_AzimuthalQuantumNumber.CartesianMagneticQuantumNumberCount();
		}

		Eigen::Index PrimitiveSphericalOrbitalCount() const noexcept
		{
			return PrimitiveShellCount() * m_AzimuthalQuantumNumber.MagneticQuantumNumberCount();
		}

		Eigen::Index ContractedCartesianOrbitalCount() const noexcept
		{
			return ContractedShellCount() * m_AzimuthalQuantumNumber.CartesianMagneticQuantumNumberCount();
		}

		Eigen::Index ContractedSphericalOrbitalCount() const noexcept
		{
			return ContractedShellCount() * m_AzimuthalQuantumNumber.MagneticQuantumNumberCount();
		}

		AzimuthalQuantumNumber AngularMomentum() const noexcept
		{
			return m_AzimuthalQuantumNumber;
		}

	private:
		/* CRTP PURE VIRTUAL */ const auto& ExponentSet_Impl() const noexcept = delete;

		/* CRTP PURE VIRTUAL */ const auto& ContractionSets_Impl() const noexcept = delete;

		bool EqualsTo_Impl(const AbstractAzimuthalShell& other, const Scalar tolerance) const noexcept
		{
			return PrimitiveShellCount() == other.PrimitiveShellCount()
			       && ContractedShellCount() == other.ContractedShellCount()
			       && (ExponentSet() - other.ExponentSet()).cwiseAbs().maxCoeff() <= tolerance
			       && (ContractionSets() - other.ContractionSets()).cwiseAbs().maxCoeff() <= tolerance;
		}

		/* CRTP VIRTUAL */ Eigen::Index PrimitiveShellCount_Impl() const noexcept
		{
			return ExponentSet().size();
		}

		/* CRTP VIRTUAL */ Eigen::Index ContractedShellCount_Impl() const noexcept
		{
			return ContractionSets().cols();
		}

		explicit constexpr AbstractAzimuthalShell(const AzimuthalQuantumNumber l) noexcept : m_AzimuthalQuantumNumber(l)
		{
			/* NO CODE */
		}

		constexpr AbstractAzimuthalShell(const AbstractAzimuthalShell&) noexcept = default;
		constexpr AbstractAzimuthalShell(AbstractAzimuthalShell&&) noexcept = default;
		constexpr AbstractAzimuthalShell& operator=(const AbstractAzimuthalShell&) noexcept = default;
		constexpr AbstractAzimuthalShell& operator=(AbstractAzimuthalShell&&) noexcept = default;

#if __cplusplus >= 202002L
		constexpr
#endif
		        ~AbstractAzimuthalShell() noexcept = default;

		AzimuthalQuantumNumber m_AzimuthalQuantumNumber;
	};

	class AzimuthalShellSegmentView : public AbstractAzimuthalShell<AzimuthalShellSegmentView>
	{
		using Base = AbstractAzimuthalShell;
		friend Base;
		friend AzimuthalShell;
		using Scalar = SecChem::Scalar;

		AzimuthalShellSegmentView(const AzimuthalQuantumNumber azimuthalQuantumNumber,
		                          const Eigen::Map<const Eigen::VectorX<Scalar>>& exponentSet,
		                          const Eigen::Block<const Eigen::MatrixX<Scalar>>& contractionSets) noexcept
		    : Base(azimuthalQuantumNumber), m_ExponentSetView(exponentSet), m_ContractionSetsView(contractionSets)
		{
			/* NO CODE */
		}

		const auto& ExponentSet_Impl() const noexcept
		{
			return m_ExponentSetView;
		}

		const auto& ContractionSets_Impl() const noexcept
		{
			return m_ContractionSetsView;
		}

		Eigen::Map<const Eigen::VectorX<SecChem::Scalar>> m_ExponentSetView;
		Eigen::Block<const Eigen::MatrixX<Scalar>> m_ContractionSetsView;
	};

	class AzimuthalShell : public SecUtility::IEquatableWithTolerance<AzimuthalShell>,
	                       public AbstractAzimuthalShell<AzimuthalShell>
	{
		friend IEquatableWithTolerance<AzimuthalShell>;
		using Base = AbstractAzimuthalShell;
		friend Base;
		using SegmentationTable = std::vector<std::pair<Eigen::Index, Eigen::Index>>;

	public:
		AzimuthalShell(const AzimuthalQuantumNumber angularMomentum,
		               ContractedRadialOrbitalSet contractedRadialOrbitalSet)
		    : Base(angularMomentum), m_ContractedRadialOrbitalSet(std::move(contractedRadialOrbitalSet)),
		      m_SegmentationTable{{0, 0},
		                          {m_ContractedRadialOrbitalSet.PrimitiveShellCount(),
		                           m_ContractedRadialOrbitalSet.ContractedShellCount()}}
		{
			/* NO CODE */
		}

		AzimuthalShell& AddOrOverrideContractedRadialOrbitalSet(ContractedRadialOrbitalSet contractedRadialOrbitalSet)
		{
			m_ContractedRadialOrbitalSet = std::move(contractedRadialOrbitalSet);
			m_SegmentationTable = {{0, 0},
			                       {m_ContractedRadialOrbitalSet.PrimitiveShellCount(),
			                        m_ContractedRadialOrbitalSet.ContractedShellCount()}};
			return *this;
		}

		bool IsEmpty() const noexcept
		{
			return m_ContractedRadialOrbitalSet.PrimitiveShellCount() == 0;
		}

		bool IsNotEmpty() const noexcept
		{
			return m_ContractedRadialOrbitalSet.PrimitiveShellCount() != 0;
		}

		// GCC require elaborated type specifier
		// ReSharper disable once CppRedundantElaboratedTypeSpecifier
		const class ContractedRadialOrbitalSet& ContractedRadialOrbitalSet() const noexcept
		{
			return m_ContractedRadialOrbitalSet;
		}

		/// <summary>
		/// Concatenates multiple AzimuthalShell instances.
		/// Throws std::runtime_error if the input range is empty or if blocks have different angular momenta.
		/// </summary>
		template <typename ForwardIterator, typename Getter>
		static AzimuthalShell Concat(ForwardIterator begin, const ForwardIterator end, Getter get)
		{
			if (begin == end)
			{
				throw std::runtime_error("AzimuthalShell: input range is empty");
			}

			if (std::distance(begin, end) == 1)
			{
				return get(*begin);
			}

			const auto azimuthalQuantumNumber = get(*begin).m_AzimuthalQuantumNumber;

			if (std::any_of(begin,
			                end,
			                [get, azimuthalQuantumNumber](const auto& item)
			                { return get(item).m_AzimuthalQuantumNumber != azimuthalQuantumNumber; }))
			{
				throw std::runtime_error("AzimuthalShell: cannot concat blocks with different angular momenta");
			}

			// clang-format off
			auto crs = ContractedRadialOrbitalSet::Concat(
					begin, end, [get](const auto& block) -> const auto& { return get(block).m_ContractedRadialOrbitalSet; });
			auto segTable = ConcatSegmentationTable(
					begin, end, [get](const auto& block) -> const auto& { return get(block).m_SegmentationTable; });
			// clang-format on

			return {azimuthalQuantumNumber, std::move(crs), std::move(segTable)};
		}

		/// <summary>
		/// Concatenates multiple AzimuthalShell instances.
		/// Throws std::runtime_error if the input range is empty or if blocks have different angular momenta.
		/// </summary>
		template <typename ForwardIterator>
		static AzimuthalShell Concat(ForwardIterator begin, const ForwardIterator end)
		{
			return Concat(begin, end, [](auto&& item) -> decltype(auto) { return std::forward<decltype(item)>(item); });
		}

#if __cplusplus >= 202002L
		constexpr Eigen::Index SegmentCount() const noexcept
#else
		Eigen::Index SegmentCount() const noexcept
#endif
		{
			assert(IsNotEmpty());
			return static_cast<Eigen::Index>(m_SegmentationTable.size() - 1);
		}

		AzimuthalShellSegmentView Segment(const Eigen::Index index) const noexcept
		{
			assert(IsNotEmpty());
			assert(index >= 0 && index < SegmentCount());

			const auto [p0, c0] = m_SegmentationTable[index];
			const auto [p1, c1] = m_SegmentationTable[index + 1];

			return AzimuthalShellSegmentView{AngularMomentum(),
			                                 {ExponentSet().data() + p0, p1 - p0},
			                                 ContractionSets().block(p0, c0, p1 - p0, c1 - c0)};
		}

		using IEquatableWithTolerance<AzimuthalShell>::EqualsTo;
		using IEquatableWithTolerance<AzimuthalShell>::NotEqualsTo;
		using IEquatableWithTolerance<AzimuthalShell>::operator==;
		using IEquatableWithTolerance<AzimuthalShell>::operator!=;

	private:
		AzimuthalShell(const AzimuthalQuantumNumber angularMomentum,
		               Gaussian::ContractedRadialOrbitalSet nullableContractedRadialOrbitalSet,
		               std::vector<std::pair<Eigen::Index, Eigen::Index>> segmentationTable)
		    : Base(angularMomentum), m_ContractedRadialOrbitalSet(std::move(nullableContractedRadialOrbitalSet)),
		      m_SegmentationTable(std::move(segmentationTable))
		{
			/* NO CODE */
		}

		bool EqualsTo_Impl(const AzimuthalShell& other, const SecChem::Scalar tolerance) const noexcept
		{
			return m_ContractedRadialOrbitalSet.EqualsTo(other.m_ContractedRadialOrbitalSet, tolerance);
			// return m_ContractedRadialOrbitalSet.EqualsTo(other.m_ContractedRadialOrbitalSet, tolerance)
			//        && m_SegmentationTable == other.m_SegmentationTable;
		}

		/* CRTP OVERRIDE */ const Eigen::VectorXd& ExponentSet_Impl() const noexcept
		{
			return ContractedRadialOrbitalSet().ExponentSet();
		}

		/* CRTP OVERRIDE */ const Eigen::MatrixXd& ContractionSets_Impl() const noexcept
		{
			return ContractedRadialOrbitalSet().ContractionSets();
		}

		/* CRTP OVERRIDE */ Eigen::Index PrimitiveShellCount_Impl() const noexcept
		{
			return IsNotEmpty() ? ExponentSet_Impl().size() : 0;
		}

		/* CRTP OVERRIDE */ Eigen::Index ContractedShellCount_Impl() const noexcept
		{
			return IsNotEmpty() ? ContractionSets_Impl().cols() : 0;
		}

		template <typename InputIterator, typename Getter>
		static SegmentationTable ConcatSegmentationTable(const InputIterator begin, const InputIterator end, Getter get)
		{
			static_assert(std::is_lvalue_reference_v<decltype(get(*begin))>,
			              "Getter must return a reference type (lvalue reference)");
			const auto segmentCount = std::accumulate(begin,
			                                          end,
			                                          std::size_t{1},  // start from 1, not 0
			                                          [get](const std::size_t acc, const auto& segmentationTable)
			                                          {
				                                          if (get(segmentationTable).empty())
				                                          {
					                                          return acc;
				                                          }
				                                          return acc + get(segmentationTable).size() - 1;
			                                          });
			if (segmentCount == 1)
			{
				return {};
			}

			SegmentationTable result;
			result.reserve(segmentCount);

			std::pair<Eigen::Index, Eigen::Index> offsetPair = {0, 0};
			result.emplace_back(offsetPair);
			for (auto it = begin; it != end; ++it)
			{
				const SegmentationTable& segmentationTable = get(*it);

				if (segmentationTable.empty())
				{
					continue;
				}

				std::transform(std::next(segmentationTable.begin()),
				               segmentationTable.end(),
				               std::back_inserter(result),
				               [offsetPair](const std::pair<Eigen::Index, Eigen::Index>& seg)
				               { return std::pair{offsetPair.first + seg.first, offsetPair.second + seg.second}; });
				offsetPair = result.back();
			}

			return result;
		}

		Gaussian::ContractedRadialOrbitalSet m_ContractedRadialOrbitalSet = {};
		SegmentationTable m_SegmentationTable = {};
	};
}  // namespace SecChem::BasisSet::Gaussian
