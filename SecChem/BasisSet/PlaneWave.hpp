// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Andy Brown

#pragma once

#include <Eigen/Dense>
#include <SecChem/System.hpp>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>


namespace SecChem::BasisSet::PlaneWave
{
	class BasisSet
	{
		friend Builder<BasisSet>;

	public:
		template <typename TRealSpaceLatticeVectors>
		BasisSet(const Scalar cutoffEnergy, TRealSpaceLatticeVectors&& realSpaceLatticeVectors)
		    : m_CutoffEnergy(cutoffEnergy),
		      m_RealSpaceLatticeVectors(std::forward<TRealSpaceLatticeVectors>(realSpaceLatticeVectors))
		{
			PopulateThis();
		}

		template <typename... TRealSpaceLatticeVectors>
		BasisSet(const Scalar cutoffEnergy, TRealSpaceLatticeVectors&&... realSpaceLatticeVectors)
		    : m_CutoffEnergy(cutoffEnergy),
		      m_RealSpaceLatticeVectors({std::forward<TRealSpaceLatticeVectors>(realSpaceLatticeVectors)...})
		{
			static_assert(sizeof...(TRealSpaceLatticeVectors) == 3);
			PopulateThis();
		}

		BasisSet(const Scalar cutoffEnergy, const Scalar dimension)
		    : m_CutoffEnergy(cutoffEnergy), m_RealSpaceLatticeVectors{Eigen::Vector3<Scalar>{dimension, 0, 0},
		                                                              Eigen::Vector3<Scalar>{0, dimension, 0},
		                                                              Eigen::Vector3<Scalar>{0, 0, dimension}}
		{
			PopulateThis();
		}

		Scalar UnitCellVolume() const noexcept
		{
			return m_UnitCellVolume;
		}

		const Eigen::Vector3<Scalar>& RealSpaceLatticeVector(const std::size_t index) const noexcept
		{
			assert(index < 3);
			return m_RealSpaceLatticeVectors[index];
		}

		const Eigen::Vector3<Scalar>& ReciprocalLatticeVector(const std::size_t index) const noexcept
		{
			assert(index < 3);
			return m_ReciprocalLatticeVectors[index];
		}

		const Eigen::Vector3<Eigen::Index>& ReciprocalIndexUpperBoundTriplet() const noexcept
		{
			return m_ReciprocalIndexUpperBoundTriplet;
		}


	private:
		void PopulateThis()
		{
			m_UnitCellVolume =
			        m_RealSpaceLatticeVectors[0].dot(m_RealSpaceLatticeVectors[1].cross(m_RealSpaceLatticeVectors[2]));

			const auto factor = static_cast<Scalar>(2 * 3.1415926535897932384626433832795) / m_UnitCellVolume;
			m_ReciprocalLatticeVectors[0] = factor * m_RealSpaceLatticeVectors[1].cross(m_RealSpaceLatticeVectors[2]);
			m_ReciprocalLatticeVectors[1] = factor * m_RealSpaceLatticeVectors[2].cross(m_RealSpaceLatticeVectors[0]);
			m_ReciprocalLatticeVectors[2] = factor * m_RealSpaceLatticeVectors[0].cross(m_RealSpaceLatticeVectors[1]);

			const auto sqrtOfTwoTimesTheCutoffEnergy = std::sqrt(2 * m_CutoffEnergy);
			m_ReciprocalIndexUpperBoundTriplet[0] = static_cast<Eigen::Index>(
			        std::ceil(sqrtOfTwoTimesTheCutoffEnergy / m_ReciprocalLatticeVectors[0].norm()));
			m_ReciprocalIndexUpperBoundTriplet[1] = static_cast<Eigen::Index>(
			        std::ceil(sqrtOfTwoTimesTheCutoffEnergy / m_ReciprocalLatticeVectors[1].norm()));
			m_ReciprocalIndexUpperBoundTriplet[2] = static_cast<Eigen::Index>(
			        std::ceil(sqrtOfTwoTimesTheCutoffEnergy / m_ReciprocalLatticeVectors[2].norm()));
		}


	private:
		Scalar m_CutoffEnergy = std::numeric_limits<Scalar>::signaling_NaN();
		std::array<Eigen::Vector3<Scalar>, 3> m_RealSpaceLatticeVectors;
		std::array<Eigen::Vector3<Scalar>, 3> m_ReciprocalLatticeVectors;
		Scalar m_UnitCellVolume = std::numeric_limits<Scalar>::signaling_NaN();
		Eigen::Vector3<Eigen::Index> m_ReciprocalIndexUpperBoundTriplet;
	};
}  // namespace SecChem::BasisSet::PlaneWave
