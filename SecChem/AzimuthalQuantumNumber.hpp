// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Andy Brown

#pragma once

#include <iostream>
#include <string>

#include "Utility/StringUtility.hpp"


namespace SecChem
{
	class AzimuthalQuantumNumber
	{
		enum class AzimuthalQuantumNumberId : int
		{
			// symbols has in uppercase to conform naming convention
			S,
			P,
			D,
			F,
			G,
			H,
			I,
			J,
			K,
			L,
			M,
			N,
			O,

			Sharp = S,
			Principal = P,
			Diffuse = D,
			Fundamental = F
		};


	public:
		static constexpr auto S = AzimuthalQuantumNumberId::S;
		static constexpr auto P = AzimuthalQuantumNumberId::P;
		static constexpr auto D = AzimuthalQuantumNumberId::D;
		static constexpr auto F = AzimuthalQuantumNumberId::F;
		static constexpr auto G = AzimuthalQuantumNumberId::G;
		static constexpr auto H = AzimuthalQuantumNumberId::H;
		static constexpr auto I = AzimuthalQuantumNumberId::I;
		static constexpr auto J = AzimuthalQuantumNumberId::J;
		static constexpr auto K = AzimuthalQuantumNumberId::K;
		static constexpr auto L = AzimuthalQuantumNumberId::L;
		static constexpr auto M = AzimuthalQuantumNumberId::M;
		static constexpr auto N = AzimuthalQuantumNumberId::N;
		static constexpr auto O = AzimuthalQuantumNumberId::O;
		static constexpr auto Sharp = AzimuthalQuantumNumberId::Sharp;
		static constexpr auto Principal = AzimuthalQuantumNumberId::Principal;
		static constexpr auto Diffuse = AzimuthalQuantumNumberId::Diffuse;
		static constexpr auto Fundamental = AzimuthalQuantumNumberId::Fundamental;
		static constexpr auto SupportedAzimuthalQuantumNumberCount = static_cast<int>(AzimuthalQuantumNumberId::O) + 1;


		constexpr AzimuthalQuantumNumber() = default;

		// ReSharper disable once CppNonExplicitConvertingConstructor
		/* IMPLICIT */ constexpr AzimuthalQuantumNumber(  // NOLINT(*-explicit-constructor)
		        const AzimuthalQuantumNumberId id)
		    : m_AzimuthalQuantumNumberValue(id)
		{
			/* NO CODE */
		}

		explicit constexpr AzimuthalQuantumNumber(int value)
		    : m_AzimuthalQuantumNumberValue(static_cast<AzimuthalQuantumNumberId>(value))
		{
			/* NO CODE */
		}

		explicit constexpr AzimuthalQuantumNumber(const char label)
		    : m_AzimuthalQuantumNumberValue(
		              static_cast<AzimuthalQuantumNumberId>(FuzzySubshellLabel2AzimuthalQuantumNumberValue(label)))
		{
			/* NO CODE */
		}

		explicit constexpr operator int() const noexcept  // NOLINT(*-explicit-constructor)
		{
			return this->Value();
		}

		constexpr int Value() const noexcept
		{
			return static_cast<int>(m_AzimuthalQuantumNumberValue);
		}

		friend std::ostream& operator<<(std::ostream& os, const AzimuthalQuantumNumber azimuthalQuantumNumber)
		{
			return os << azimuthalQuantumNumber.Label();
		}

		friend std::istream& operator>>(std::istream& is, AzimuthalQuantumNumber& out_azimuthalQuantumNumber)
		{
			char l = '\0';
			is >> l;
			out_azimuthalQuantumNumber = AzimuthalQuantumNumber{l};
			return is;
		}

		AzimuthalQuantumNumber operator++()
		{
			auto temp = static_cast<std::underlying_type_t<AzimuthalQuantumNumberId>>(m_AzimuthalQuantumNumberValue);
			++temp;
			m_AzimuthalQuantumNumberValue = static_cast<AzimuthalQuantumNumberId>(temp);
			return *this;
		}

		AzimuthalQuantumNumber operator++(int)
		{
			// ReSharper disable once CppLocalVariableMayBeConst
			auto copy = *this;
			++*this;
			return copy;
		}

		AzimuthalQuantumNumber operator--()
		{
			auto temp = static_cast<std::underlying_type_t<AzimuthalQuantumNumberId>>(m_AzimuthalQuantumNumberValue);
			--temp;
			m_AzimuthalQuantumNumberValue = static_cast<AzimuthalQuantumNumberId>(temp);
			return *this;
		}

		AzimuthalQuantumNumber operator--(int)
		{
			// ReSharper disable once CppLocalVariableMayBeConst
			auto copy = *this;
			--*this;
			return copy;
		}

		void operator+=(const AzimuthalQuantumNumber delta)
		{
			m_AzimuthalQuantumNumberValue = static_cast<AzimuthalQuantumNumberId>(Value() + delta.Value());
		}

		void operator-=(const AzimuthalQuantumNumber delta)
		{
			m_AzimuthalQuantumNumberValue = static_cast<AzimuthalQuantumNumberId>(Value() - delta.Value());
		}

#define SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(LHS, RHS, OP)                                             \
	friend constexpr bool operator OP(LHS lhs, RHS rhs)                                                                \
	{                                                                                                                  \
		return static_cast<int>(lhs) OP static_cast<int>(rhs);                                                         \
	}
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumber, AzimuthalQuantumNumber, ==)
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumber, AzimuthalQuantumNumber, !=)
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumber, AzimuthalQuantumNumber, >)
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumber, AzimuthalQuantumNumber, <)
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumber, AzimuthalQuantumNumber, >=)
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumber, AzimuthalQuantumNumber, <=)
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumber, AzimuthalQuantumNumberId, ==)
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumber, AzimuthalQuantumNumberId, !=)
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumber, AzimuthalQuantumNumberId, >)
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumber, AzimuthalQuantumNumberId, <)
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumber, AzimuthalQuantumNumberId, >=)
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumber, AzimuthalQuantumNumberId, <=)
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumberId, AzimuthalQuantumNumber, ==)
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumberId, AzimuthalQuantumNumber, !=)
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumberId, AzimuthalQuantumNumber, >)
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumberId, AzimuthalQuantumNumber, <)
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumberId, AzimuthalQuantumNumber, >=)
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumberId, AzimuthalQuantumNumber, <=)
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumberId, AzimuthalQuantumNumberId, >)
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumberId, AzimuthalQuantumNumberId, <)
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumberId, AzimuthalQuantumNumberId, >=)
		SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR(AzimuthalQuantumNumberId, AzimuthalQuantumNumberId, <=)
#undef SECCHEM_AZIMUTHAL_QUANTUM_NUMBER_COMPARISON_OPERATOR

		constexpr int MinPrincipalQuantumNumber() const noexcept
		{
			return static_cast<int>(m_AzimuthalQuantumNumberValue) + 1;
		}

		constexpr int MagneticQuantumNumberCount() const noexcept
		{
			return 2 * this->Value() + 1;
		}

		constexpr int MaxMagneticQuantumNumber() const noexcept
		{
			return this->Value();
		}

		constexpr int MinMagneticQuantumNumber() const noexcept
		{
			return -MaxMagneticQuantumNumber();
		}

		constexpr int CartesianMagneticQuantumNumberCount() const noexcept
		{
			return (this->Value() + 1) * (this->Value() + 2) / 2;
		}

		constexpr int Capacity() const noexcept  // as in how much electron can the subshell contain
		{
			return MagneticQuantumNumberCount() * 2;
		}

		std::string Name() const
		{
			switch (Value())
			{
				case 0:
					return "Sharp";
				case 1:
					return "Principal";
				case 2:
					return "Diffuse";
				case 3:
					return "Fundamental";
				default:
					char name[] = "Shell?\0";
					name[5] = static_cast<char>('F' + static_cast<int>(*this) - 3);
					return name;
			}
		}

		constexpr char Label() const noexcept
		{
			switch (Value())
			{
				case 0:
					return 's';
				case 1:
					return 'p';
				case 2:
					return 'd';
				default:
					return static_cast<char>('f' + static_cast<int>(*this) - 3);
			}
		}

		constexpr bool IsSharp() const noexcept
		{
			return m_AzimuthalQuantumNumberValue == S;
		}

		constexpr bool IsPrincipal() const noexcept
		{
			return m_AzimuthalQuantumNumberValue == P;
		}

		constexpr bool IsDiffuse() const noexcept
		{
			return m_AzimuthalQuantumNumberValue == D;
		}

		constexpr bool IsFundamental() const noexcept
		{
			return m_AzimuthalQuantumNumberValue == F;
		}

		static constexpr int SubshellLabel2AzimuthalQuantumNumberValue(const char label)
		{
			switch (label)
			{
				case 's':
					return 0;
				case 'p':
					return 1;
				case 'd':
					return 2;
				default:
				{
					if (label >= 'f' && label <= 'o')
					{
						return label + 3 - static_cast<int>('f');
					}
					throw std::runtime_error("Invalid electronic subshell label");
				}
			}
		}

		static constexpr int FuzzySubshellLabel2AzimuthalQuantumNumberValue(const char label)
		{
			return SubshellLabel2AzimuthalQuantumNumberValue(static_cast<char>(SecUtility::AsciiToLower(label)));
		}

	private:
		AzimuthalQuantumNumberId m_AzimuthalQuantumNumberValue = S;
	};
}  // namespace SecChem
