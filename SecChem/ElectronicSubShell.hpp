//
// Created by Andy on 11/22/2025.
//

#pragma once

#include "AzimuthalQuantumNumber.hpp"
#include <cassert>


namespace SecChem
{
	class ElectronicSubShell
	{
	public:
		ElectronicSubShell() = default;

		constexpr ElectronicSubShell(const int principalQuantumNumber,
		                             const AzimuthalQuantumNumber azimuthalQuantumNumber)
		    : m_Id(GenerateId(principalQuantumNumber, azimuthalQuantumNumber))
		{
			assert(principalQuantumNumber >= 1 && "Principal quantum number shall not less than 1");
			assert(principalQuantumNumber >= azimuthalQuantumNumber.MinPrincipalQuantumNumber()
			       && "Invalid principal-azimuthal quantum number combination");
		}

		constexpr ElectronicSubShell(const int principalQuantumNumber, const int azimuthalQuantumNumber)
		    : m_Id(GenerateId(principalQuantumNumber, SecChem::AzimuthalQuantumNumber(azimuthalQuantumNumber)))
		{
			assert(principalQuantumNumber >= 1 && "Principal quantum number shall not less than 1");
			assert(azimuthalQuantumNumber >= 0 && "Azimuthal quantum number shall not less than 0");
			assert(principalQuantumNumber >= this->AzimuthalQuantumNumber().MinPrincipalQuantumNumber()
			       && "Invalid principal-azimuthal quantum number combination");
		}

		constexpr ElectronicSubShell(const int principalQuantumNumber, const char subshellLabel)
		    : m_Id(GenerateId(principalQuantumNumber, SecChem::AzimuthalQuantumNumber(subshellLabel)))
		{
			assert(principalQuantumNumber >= 1 && "Principal quantum number shall not less than 1");
			assert(principalQuantumNumber >= this->AzimuthalQuantumNumber().MinPrincipalQuantumNumber()
			       && "Invalid principal-azimuthal quantum number combination");
		}

		explicit constexpr ElectronicSubShell(const int id) : m_Id(id)
		{
			/* NO CODE */
		}


		constexpr int UnderlyingId() const
		{
			return m_Id;
		}

		constexpr bool IsValid() const
		{
			// return m_Id >= 0;
			return m_Id >= 0
			       && this->PrincipalQuantumNumber() >= this->AzimuthalQuantumNumber().MinPrincipalQuantumNumber();
		}

		constexpr int PrincipalQuantumNumber() const
		{
			return m_Id / AzimuthalQuantumNumber::SupportedAzimuthalQuantumNumberCount + 1;
		}

		std::string SubshellName() const
		{
			return AzimuthalQuantumNumber().Name();
		}

		std::string ToString() const
		{
			return std::to_string(PrincipalQuantumNumber()) + SubshellLabel();
		}

		friend std::ostream& operator<<(std::ostream& os, const ElectronicSubShell shell)
		{
			return os << shell.ToString();
		}

		friend std::istream& operator>>(std::istream& is, ElectronicSubShell& out_shell)
		{
			int n = 0;
			SecChem::AzimuthalQuantumNumber l;
			is >> n >> l;
			out_shell = ElectronicSubShell(n, l);
			return is;
		}

		friend constexpr bool operator<(const ElectronicSubShell lhs, const ElectronicSubShell rhs)
		{
			return lhs.m_Id < rhs.m_Id;
		}

		friend constexpr bool operator>(const ElectronicSubShell lhs, const ElectronicSubShell rhs)
		{
			return lhs.m_Id > rhs.m_Id;
		}

		friend constexpr bool operator<=(const ElectronicSubShell lhs, const ElectronicSubShell rhs)
		{
			return lhs.m_Id <= rhs.m_Id;
		}

		friend constexpr bool operator>=(const ElectronicSubShell lhs, const ElectronicSubShell rhs)
		{
			return lhs.m_Id >= rhs.m_Id;
		}

		friend constexpr bool operator==(const ElectronicSubShell lhs, const ElectronicSubShell rhs)
		{
			return lhs.m_Id == rhs.m_Id;
		}

		friend constexpr bool operator!=(const ElectronicSubShell lhs, const ElectronicSubShell rhs)
		{
			return lhs.m_Id != rhs.m_Id;
		}

		constexpr ElectronicSubShell& operator++()
		{
			/*
			 1s
			 2s  2p
			 3s  3p  3d
			 4s  4p  4d  4f
			 5s  5p  5d  5f  5g
			 6s  6p  6d  6f  6g  6h
			 7s  7p  7d  7f  7g  7h  7i
			 8s  8p  8d  8f  8g  8h  8i  8j
			 9s  9p  9d  9f  9g  9h  9i  9j  9k
			10s 10p 10d 10f 10g 10h 10i 10j 10k 10l
			11s 11p 11d 11f 11g 11h 11i 11j 11k 11l 11m
			12s 12p 12d 12f 12g 12h 12i 12j 12k 12l 12m 12n
			13s 13p 13d 13f 13g 13h 13i 13j 13k 13l 13m 13n 13o
			14s 14p 14d 14f 14g 14h 14i 14j 14k 14l 14m 14n 14o
			15s 15p 15d 15f 15g 15h 15i 15j 15k 15l 15m 15n 15o
			16s 16p 16d 16f 16g 16h 16i 16j 16k 16l 16m 16n 16o
			17s 17p 17d 17f 17g 17h 17i 17j 17k 17l 17m 17n 17o
			18s 18p 18d 18f 18g 18h 18i 18j 18k 18l 18m 18n 18o
			19s 19p 19d 19f 19g 19h 19i 19j 19k 19l 19m 19n 19o
			 :   :   :   :   :   :   :   :   :   :   :   :   :
			 */
			do
			{
				m_Id++;
			// } while (this->PrincipalQuantumNumber() < this->AzimuthalQuantumNumber().MinPrincipalQuantumNumber());
			} while (!this->IsValid());
			return *this;
		}

		ElectronicSubShell operator++(int)
		{
			// ReSharper disable once CppLocalVariableMayBeConst
			auto tmp = *this;
			++*this;
			return tmp;
		}

		constexpr ElectronicSubShell& operator--()
		{
			do
			{
				m_Id--;
			} while (m_Id >= 0 && !this->IsValid() /*this->PrincipalQuantumNumber() < this->AzimuthalQuantumNumber().MinPrincipalQuantumNumber()
			         && this->IsValid()*/);
			return *this;
		}

		ElectronicSubShell operator--(int)
		{
			// ReSharper disable once CppLocalVariableMayBeConst
			auto tmp = *this;
			--*this;
			return tmp;
		}

		constexpr double SlaterEffectivePrincipalQuantumNumber() const
		{
			switch (PrincipalQuantumNumber())
			{
				case 1:
				case 2:
				case 3:
					return PrincipalQuantumNumber();
				case 4:
					return 3.7;
				case 5:
					return 4.0;
				case 6:
					return 4.2;
				case 7:
					return 4.3;
				default:
					throw std::runtime_error(
					        "Slater effective principal quantum number is only available for first 7 periods");
			}
		}


		constexpr AzimuthalQuantumNumber AzimuthalQuantumNumber() const
		{
			return SecChem::AzimuthalQuantumNumber{m_Id % AzimuthalQuantumNumber::SupportedAzimuthalQuantumNumberCount};
		}

		constexpr int MagneticQuantumNumberCount() const noexcept
		{
			return AzimuthalQuantumNumber().MagneticQuantumNumberCount();
		}

		constexpr int MaxMagneticQuantumNumber() const noexcept
		{
			return AzimuthalQuantumNumber().MaxMagneticQuantumNumber();
		}

		constexpr int MinMagneticQuantumNumber() const noexcept
		{
			return AzimuthalQuantumNumber().MinMagneticQuantumNumber();
		}

		constexpr int CartesianMagneticQuantumNumberCount() const noexcept
		{
			return AzimuthalQuantumNumber().CartesianMagneticQuantumNumberCount();
		}

		constexpr int Capacity() const noexcept
		{
			return AzimuthalQuantumNumber().Capacity();
		}

		constexpr char SubshellLabel() const
		{
			return AzimuthalQuantumNumber().Label();
		}

		constexpr bool IsSharp() const noexcept
		{
			return AzimuthalQuantumNumber().IsSharp();
		}

		constexpr bool IsPrincipal() const noexcept
		{
			return AzimuthalQuantumNumber().IsPrincipal();
		}

		constexpr bool IsDiffuse() const noexcept
		{
			return AzimuthalQuantumNumber().IsDiffuse();
		}

		constexpr bool IsFundamental() const noexcept
		{
			return AzimuthalQuantumNumber().IsFundamental();
		}



	private:
		static constexpr int GenerateId(const int principalQuantumNumber,
		                                const SecChem::AzimuthalQuantumNumber azimuthalQuantumNumber)
		{
			return (principalQuantumNumber - 1) * AzimuthalQuantumNumber::SupportedAzimuthalQuantumNumberCount
			       + static_cast<int>(azimuthalQuantumNumber);
		}


	private:
		int m_Id = 0;
	};


	constexpr ElectronicSubShell operator""_Sharp(const unsigned long long principalQuantumNumber)
	{
		return ElectronicSubShell{static_cast<int>(principalQuantumNumber), 0};
	}

	constexpr ElectronicSubShell operator""_Principal(const unsigned long long principalQuantumNumber)
	{
		return ElectronicSubShell{static_cast<int>(principalQuantumNumber), 1};
	}

	constexpr ElectronicSubShell operator""_Diffuse(const unsigned long long principalQuantumNumber)
	{
		return ElectronicSubShell{static_cast<int>(principalQuantumNumber), 2};
	}

	constexpr ElectronicSubShell operator""_Fundamental(const unsigned long long principalQuantumNumber)
	{
		return ElectronicSubShell{static_cast<int>(principalQuantumNumber), 3};
	}
}  // namespace SecChem
