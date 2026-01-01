// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Andy Brown

#pragma once

#include "ElectronicSubShell.hpp"
#include <algorithm>
#include <array>


namespace SecChem
{
	class AtomicElectronConfiguration
	{
	public:
		static constexpr int MaxPrincipalQuantumNumber = 7;  // Periodic table only reach 7-th period for now
		static constexpr ElectronicSubShell MaxElectronShell{MaxPrincipalQuantumNumber, MaxPrincipalQuantumNumber - 1};

		AtomicElectronConfiguration() = default;


		constexpr explicit AtomicElectronConfiguration(int electronCount, const int diagonalPassCount = 18)
		{
			// Diagonal pass of Aufbau ordering (n + l rule)
			// The "diagonal index" corresponds to the sum of principal + azimuthal quantum numbers.
			// This loop walks the diagonal, filling shells in n+l order.
			for (int diagonalIndex = 0; diagonalIndex < diagonalPassCount; diagonalIndex++)
			{
				auto principalMinusOneTimesTwo = diagonalIndex;
				auto azimuthalTimesTwo = diagonalIndex;

				if (principalMinusOneTimesTwo % 2 != 0)
				{
					principalMinusOneTimesTwo++;
					azimuthalTimesTwo--;
				}

				while (azimuthalTimesTwo >= 0)
				{
					m_ElectronCounts[ElectronicSubShell(principalMinusOneTimesTwo / 2 + 1, azimuthalTimesTwo / 2)
					                         .UnderlyingId()] =
					        electronCount < (azimuthalTimesTwo + 1) * 2 ? electronCount : (azimuthalTimesTwo + 1) * 2;

					electronCount -= (azimuthalTimesTwo + 1) * 2;
					if (electronCount <= 0)
					{
						return;
					}
					principalMinusOneTimesTwo += 2;
					azimuthalTimesTwo -= 2;
				}
			}

			throw std::runtime_error("Loop limit reach while generating electron configuration");
		}

		AtomicElectronConfiguration(const std::initializer_list<std::pair<ElectronicSubShell, int>>& config)
		{
			for (const auto& [subShell, electronCount] : config)
			{
				if (electronCount < 0 || electronCount > subShell.Capacity())
				{
					throw std::runtime_error(
					        "Electron count was negative or excess the capacity of electronic subshell");
				}

				m_ElectronCounts[subShell.UnderlyingId()] = electronCount;
			}
		}

		constexpr bool operator==(const AtomicElectronConfiguration& other) const
		{
			for (std::size_t i = 0; i < m_ElectronCounts.size(); i++)
			{
				if (m_ElectronCounts[i] != other.m_ElectronCounts[i])
				{
					return false;
				}
			}

			return true;
		}

		constexpr bool operator!=(const AtomicElectronConfiguration& other) const
		{
			return !(*this == other);
		}

		constexpr int operator[](const ElectronicSubShell shell) const
		{
			return m_ElectronCounts[shell.UnderlyingId()];
		}

		constexpr void MoveElectron(const ElectronicSubShell source, const ElectronicSubShell destination)
		{
			MoveElectrons(source, destination, 1);
		}

		constexpr void MoveElectrons(const ElectronicSubShell source, const ElectronicSubShell destination, int amount)
		{
			amount = std::min(m_ElectronCounts[source.UnderlyingId()], amount);
			amount = std::min(destination.Capacity() - m_ElectronCounts[destination.UnderlyingId()], amount);

			if (amount < 1)
			{
				return;
			}

			m_ElectronCounts[source.UnderlyingId()] -= amount;
			m_ElectronCounts[destination.UnderlyingId()] += amount;
		}

		friend std::ostream& operator<<(std::ostream& os, const AtomicElectronConfiguration& config)
		{
			os << ElectronicSubShell(0) << '(' << config.m_ElectronCounts[0]
			   << ") ";  // NOTE: We always prints 1s shell

			// NOTE: `i` starts from 1 as 1s shell was printed earlier
			for (int i = 1; i <= MaxElectronShell.UnderlyingId(); i++)
			{
				if (const auto n = config.m_ElectronCounts[i]; n > 0)
				{
					os << ElectronicSubShell(i) << '(' << n << ") ";
				}
			}

			return os;
		}

	private:
		std::array<int, MaxElectronShell.UnderlyingId() + 1> m_ElectronCounts{};
	};
}  // namespace SecChem
