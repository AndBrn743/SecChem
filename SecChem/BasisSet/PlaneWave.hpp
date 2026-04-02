// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Andy Brown

#pragma once

#include <Eigen/Dense>
#include <SecChem/System.hpp>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <utility>


namespace SecChem::BasisSet::PlaneWave
{
	class BasisSet
	{
		friend Builder<BasisSet>;

		class PlaneWaves
		{
			friend BasisSet;

			using Iterable = Eigen::VectorwiseOp<const Eigen::Matrix<Scalar, 3, Eigen::Dynamic>, Eigen::Vertical>;

		public:
			PlaneWaves(const PlaneWaves& other)
			    : m_KineticEnergies(other.m_KineticEnergies), m_WaveVectors(other.m_WaveVectors),
			      m_RealSpaceGridPoints(other.m_RealSpaceGridPoints)
			{
				CreateOrSyncIterables();
				// NO CODE
			}

			PlaneWaves& operator=(const PlaneWaves& other) &
			{
				m_KineticEnergies = other.m_KineticEnergies;
				m_WaveVectors = other.m_WaveVectors;
				m_RealSpaceGridPoints = other.m_RealSpaceGridPoints;
				CreateOrSyncIterables();
				return *this;
			}

			PlaneWaves(PlaneWaves&& other) noexcept
			    : m_KineticEnergies(std::move(other).m_KineticEnergies), m_WaveVectors(std::move(other).m_WaveVectors),
			      m_RealSpaceGridPoints(std::move(other).m_RealSpaceGridPoints)
			{
				CreateOrSyncIterables();
			}

			PlaneWaves& operator=(PlaneWaves&& other) noexcept
			{
				m_KineticEnergies = std::move(other).m_KineticEnergies;
				m_WaveVectors = std::move(other).m_WaveVectors;
				m_RealSpaceGridPoints = std::move(other).m_RealSpaceGridPoints;
				CreateOrSyncIterables();
				return *this;
			}

			~PlaneWaves() = default;

			constexpr const auto& KineticEnergies() const noexcept
			{
				return m_KineticEnergies;
			}

			const auto& IterableWaveVectors() const noexcept
			{
				return *m_IterableWaveVectorsPtr;  // lvalue is required for to use range-v3's pipe operator
			}

			constexpr const auto& WaveVectors() const noexcept
			{
				return m_WaveVectors;
			}

			decltype(auto) WaveVector(const Eigen::Index index) const noexcept
			{
				return m_WaveVectors.col(index);
			}

			const auto& IterableRealSpaceGridPoints() const noexcept
			{
				return *m_IterableRealSpaceGridPointsPtr;  // lvalue is required for to use range-v3's pipe operator
			}

			constexpr const auto& RealSpaceGridPoints() const noexcept
			{
				return m_RealSpaceGridPoints;
			}

			decltype(auto) RealSpaceGridPoint(const Eigen::Index index) const noexcept
			{
				return m_RealSpaceGridPoints.col(index);
			}

			Eigen::Index Count() const noexcept
			{
				return m_KineticEnergies.size();
			}


		private:
			PlaneWaves() noexcept = default;

			void CreateOrSyncIterables() noexcept
			{
				m_IterableWaveVectorsPtr = std::make_unique<const Iterable>(std::as_const(m_WaveVectors).colwise());
				m_IterableRealSpaceGridPointsPtr =
				        std::make_unique<const Iterable>(std::as_const(m_RealSpaceGridPoints).colwise());
			}


		private:
			Eigen::VectorX<Scalar> m_KineticEnergies;
			Eigen::Matrix<Scalar, 3, Eigen::Dynamic> m_WaveVectors;
			Eigen::Matrix<Scalar, 3, Eigen::Dynamic> m_RealSpaceGridPoints;

			// lvalue is required for to use range-v3's pipe operator, therefore we can't lazy eval them.
			// Eigen::VectorwiseOp nests by reference, which forbids assignment, therefore we need store them w/ ptr.
			// not ideal, but it's the best thing I have for now. if anyone got better idea, I'd love to hear them
			std::unique_ptr<const Iterable> m_IterableWaveVectorsPtr = nullptr;
			std::unique_ptr<const Iterable> m_IterableRealSpaceGridPointsPtr = nullptr;
		};


	public:
		template <typename TRealSpaceLatticeVectors>
		BasisSet(const Scalar cutoffEnergy, TRealSpaceLatticeVectors&& realSpaceLatticeVectors)
		    : m_CutoffEnergy(cutoffEnergy),
		      m_RealSpaceLatticeVectors(std::forward<TRealSpaceLatticeVectors>(realSpaceLatticeVectors))
		{
			PopulateThis();
		}

		template <typename... TRealSpaceLatticeVectors>
		// ReSharper disable once CppNonExplicitConvertingConstructor
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

		const PlaneWaves& Basis() const noexcept
		{
			return m_PlaneWaves;
		}

		PlaneWaves& Basis() noexcept
		{
			return m_PlaneWaves;
		}

		Eigen::Index BasisCount() const noexcept
		{
			return m_PlaneWaves.Count();
		}


	private:
		void PopulateThis()
		{
			m_UnitCellVolume =
			        m_RealSpaceLatticeVectors[0].dot(m_RealSpaceLatticeVectors[1].cross(m_RealSpaceLatticeVectors[2]));

			// ReSharper disable once CppRedundantCastExpression
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

			m_PlaneWaves = PlaneWavesInFftOrdering();
		}

		PlaneWaves PlaneWavesInFftOrdering() const noexcept
		{
			PlaneWaves pws;
			std::vector<Scalar> kineticEnergies;
			std::vector<Eigen::Vector3<Scalar>> waveVectors;
			std::vector<Eigen::Vector3<Scalar>> rs;

			static constexpr auto wrapIndex = [](const Eigen::Index i, const Eigen::Index n)
			{ return i <= n / 2 ? i : i - n; };

			for (Eigen::Index k = 0; k < m_ReciprocalIndexUpperBoundTriplet[2]; k++)
			{
				for (Eigen::Index j = 0; j < m_ReciprocalIndexUpperBoundTriplet[1]; j++)
				{
					for (Eigen::Index i = 0; i < m_ReciprocalIndexUpperBoundTriplet[0]; i++)
					{
						const Eigen::Vector3<Eigen::Index> latticeIndex(
						        wrapIndex(i, m_ReciprocalIndexUpperBoundTriplet[0]),
						        wrapIndex(j, m_ReciprocalIndexUpperBoundTriplet[1]),
						        wrapIndex(k, m_ReciprocalIndexUpperBoundTriplet[2]));

						// `Eigen::Vector3<Scalar>` is not movable, and coping 3 `double`s is cheap anyway
						const Eigen::Vector3<Scalar> waveVector = latticeIndex[0] * m_ReciprocalLatticeVectors[0]
						                                          + latticeIndex[1] * m_ReciprocalLatticeVectors[1]
						                                          + latticeIndex[2] * m_ReciprocalLatticeVectors[2];
						const auto kineticEnergy = waveVector.squaredNorm() / 2;

						// if (kineticEnergy <= m_CutoffEnergy)
						{
							waveVectors.emplace_back(waveVector);
							rs.emplace_back(i * m_RealSpaceLatticeVectors[0] / m_ReciprocalIndexUpperBoundTriplet[0]
							                + j * m_RealSpaceLatticeVectors[1] / m_ReciprocalIndexUpperBoundTriplet[1]
							                + k * m_RealSpaceLatticeVectors[2] / m_ReciprocalIndexUpperBoundTriplet[2]);
							kineticEnergies.emplace_back(kineticEnergy);
						}
					}
				}
			}

			pws.m_KineticEnergies = Eigen::VectorX<Scalar>::Map(kineticEnergies.data(), static_cast<Eigen::Index>(kineticEnergies.size()));
			pws.m_WaveVectors = Eigen::Matrix<Scalar, 3, Eigen::Dynamic>(3, pws.m_KineticEnergies.size());
			pws.m_RealSpaceGridPoints = Eigen::Matrix<Scalar, 3, Eigen::Dynamic>(3, pws.m_KineticEnergies.size());
			for (Eigen::Index i = 0; i < pws.m_WaveVectors.cols(); i++)
			{
				pws.m_WaveVectors.col(i) = waveVectors[i];
				pws.m_RealSpaceGridPoints.col(i) = rs[i];
			}

			pws.CreateOrSyncIterables();

			return pws;
		}


	private:
		Scalar m_CutoffEnergy = std::numeric_limits<Scalar>::signaling_NaN();
		std::array<Eigen::Vector3<Scalar>, 3> m_RealSpaceLatticeVectors;
		std::array<Eigen::Vector3<Scalar>, 3> m_ReciprocalLatticeVectors;
		Scalar m_UnitCellVolume = std::numeric_limits<Scalar>::signaling_NaN();
		Eigen::Vector3<Eigen::Index> m_ReciprocalIndexUpperBoundTriplet;
		PlaneWaves m_PlaneWaves;
	};
}  // namespace SecChem::BasisSet::PlaneWave
