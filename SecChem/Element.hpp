//
// Created by Andy on 11/22/2025.
//

#pragma once

#include "AtomicElectronConfiguration.hpp"
#include "AzimuthalQuantumNumber.hpp"
#include <cassert>
#include <iostream>
#include <string>

#define CONSTEXPR23

#if __has_include(<SecUtility/IO/Stringtility.hpp>)
#include <SecUtility/IO/Stringtility.hpp>
#else
namespace SecUtility
{
	inline int CaseInsensitiveStringComparison(const char* lhs, const char* rhs)
	{
#if defined(_WIN32)
		return _stricmp(lhs, rhs);
#else
		return strcasecmp(lhs, rhs);
#endif
	}

	inline int CaseInsensitiveStringComparison(const std::string& lhs, const char* rhs)
	{
		return CaseInsensitiveStringComparison(lhs.c_str(), rhs);
	}

	inline int CaseInsensitiveStringComparison(const char* lhs, const std::string& rhs)
	{
		return CaseInsensitiveStringComparison(lhs, rhs.c_str());
	}

	inline int CaseInsensitiveStringComparison(const std::string& lhs, const std::string& rhs)
	{
		return CaseInsensitiveStringComparison(lhs.c_str(), rhs.c_str());
	}

	template <typename LhsString, typename RhsString>
	bool IsCaseInsensitiveStringEqual(LhsString&& lhs, RhsString&& rhs)
	{
		return CaseInsensitiveStringComparison(std::forward<decltype(lhs)>(lhs), std::forward<decltype(rhs)>(rhs)) == 0;
	}
}  // namespace SecUtility
#endif
#if __has_include(<SecUtility/UnitOfMeasurement.hpp>)
#include <SecUtility/UnitOfMeasurement.hpp>
#else
namespace SecUtility::UnitOfMeasurement
{
	template <typename T>
	constexpr decltype(auto) BohrRadius2Angstrom(const T length)
	{
		return length * 5.2917721090380e-1;
	}
}  // namespace SecUtility::UnitOfMeasurement
#endif


namespace SecChem
{
	class Element
	{
		enum class Id : int
		{
			Neutron,
			// clang-format off
				 H,                                                                                                                         He,
				Li, Be,                                                                                                  B,  C,  N,  O,  F, Ne,
				Na, Mg,                                                                                                 Al, Si,  P,  S, Cl, Ar,
				 K, Ca,                                                         Sc, Ti,  V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr,
				Rb, Sr,                                                          Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te,  I, Xe,
				Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta,  W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,
				Fr, Ra, Ac, Th, Pa,  U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og
			// clang-format on
		};

	public:
		static constexpr int MinAtomicNumber = 0;
		static constexpr int MaxAtomicNumber = 118;
		static constexpr int PeriodicTableElementCount = MaxAtomicNumber + 1;

#pragma region
		static constexpr auto Neutron = Id::Neutron;
		static constexpr auto Hydrogen = Id::H;
		static constexpr auto Helium = Id::He;
		static constexpr auto Lithium = Id::Li;
		static constexpr auto Beryllium = Id::Be;
		static constexpr auto Boron = Id::B;
		static constexpr auto Carbon = Id::C;
		static constexpr auto Nitrogen = Id::N;
		static constexpr auto Oxygen = Id::O;
		static constexpr auto Fluorine = Id::F;
		static constexpr auto Neon = Id::Ne;
		static constexpr auto Sodium = Id::Na;
		static constexpr auto Magnesium = Id::Mg;
		static constexpr auto Aluminum = Id::Al;
		static constexpr auto Silicon = Id::Si;
		static constexpr auto Phosphorus = Id::P;
		static constexpr auto Sulfur = Id::S;
		static constexpr auto Chlorine = Id::Cl;
		static constexpr auto Argon = Id::Ar;
		static constexpr auto Potassium = Id::K;
		static constexpr auto Calcium = Id::Ca;
		static constexpr auto Scandium = Id::Sc;
		static constexpr auto Titanium = Id::Ti;
		static constexpr auto Vanadium = Id::V;
		static constexpr auto Chromium = Id::Cr;
		static constexpr auto Manganese = Id::Mn;
		static constexpr auto Iron = Id::Fe;
		static constexpr auto Cobalt = Id::Co;
		static constexpr auto Nickel = Id::Ni;
		static constexpr auto Copper = Id::Cu;
		static constexpr auto Zinc = Id::Zn;
		static constexpr auto Gallium = Id::Ga;
		static constexpr auto Germanium = Id::Ge;
		static constexpr auto Arsenic = Id::As;
		static constexpr auto Selenium = Id::Se;
		static constexpr auto Bromine = Id::Br;
		static constexpr auto Krypton = Id::Kr;
		static constexpr auto Rubidium = Id::Rb;
		static constexpr auto Strontium = Id::Sr;
		static constexpr auto Yttrium = Id::Y;
		static constexpr auto Zirconium = Id::Zr;
		static constexpr auto Niobium = Id::Nb;
		static constexpr auto Molybdenum = Id::Mo;
		static constexpr auto Technetium = Id::Tc;
		static constexpr auto Ruthenium = Id::Ru;
		static constexpr auto Rhodium = Id::Rh;
		static constexpr auto Palladium = Id::Pd;
		static constexpr auto Silver = Id::Ag;
		static constexpr auto Cadmium = Id::Cd;
		static constexpr auto Indium = Id::In;
		static constexpr auto Tin = Id::Sn;
		static constexpr auto Antimony = Id::Sb;
		static constexpr auto Tellurium = Id::Te;
		static constexpr auto Iodine = Id::I;
		static constexpr auto Xenon = Id::Xe;
		static constexpr auto Cesium = Id::Cs;
		static constexpr auto Barium = Id::Ba;
		static constexpr auto Lanthanum = Id::La;
		static constexpr auto Cerium = Id::Ce;
		static constexpr auto Praseodymium = Id::Pr;
		static constexpr auto Neodymium = Id::Nd;
		static constexpr auto Promethium = Id::Pm;
		static constexpr auto Samarium = Id::Sm;
		static constexpr auto Europium = Id::Eu;
		static constexpr auto Gadolinium = Id::Gd;
		static constexpr auto Terbium = Id::Tb;
		static constexpr auto Dysprosium = Id::Dy;
		static constexpr auto Holmium = Id::Ho;
		static constexpr auto Erbium = Id::Er;
		static constexpr auto Thulium = Id::Tm;
		static constexpr auto Ytterbium = Id::Yb;
		static constexpr auto Lutetium = Id::Lu;
		static constexpr auto Hafnium = Id::Hf;
		static constexpr auto Tantalum = Id::Ta;
		static constexpr auto Tungsten = Id::W;
		static constexpr auto Rhenium = Id::Re;
		static constexpr auto Osmium = Id::Os;
		static constexpr auto Iridium = Id::Ir;
		static constexpr auto Platinum = Id::Pt;
		static constexpr auto Gold = Id::Au;
		static constexpr auto Mercury = Id::Hg;
		static constexpr auto Thallium = Id::Tl;
		static constexpr auto Lead = Id::Pb;
		static constexpr auto Bismuth = Id::Bi;
		static constexpr auto Polonium = Id::Po;
		static constexpr auto Astatine = Id::At;
		static constexpr auto Radon = Id::Rn;
		static constexpr auto Francium = Id::Fr;
		static constexpr auto Radium = Id::Ra;
		static constexpr auto Actinium = Id::Ac;
		static constexpr auto Thorium = Id::Th;
		static constexpr auto Protactinium = Id::Pa;
		static constexpr auto Uranium = Id::U;
		static constexpr auto Neptunium = Id::Np;
		static constexpr auto Plutonium = Id::Pu;
		static constexpr auto Americium = Id::Am;
		static constexpr auto Curium = Id::Cm;
		static constexpr auto Berkelium = Id::Bk;
		static constexpr auto Californium = Id::Cf;
		static constexpr auto Einsteinium = Id::Es;
		static constexpr auto Fermium = Id::Fm;
		static constexpr auto Mendelevium = Id::Md;
		static constexpr auto Nobelium = Id::No;
		static constexpr auto Lawrencium = Id::Lr;
		static constexpr auto Rutherfordium = Id::Rf;
		static constexpr auto Dubnium = Id::Db;
		static constexpr auto Seaborgium = Id::Sg;
		static constexpr auto Bohrium = Id::Bh;
		static constexpr auto Hassium = Id::Hs;
		static constexpr auto Meitnerium = Id::Mt;
		static constexpr auto Darmstadtium = Id::Ds;
		static constexpr auto Roentgenium = Id::Rg;
		static constexpr auto Copernicium = Id::Cn;
		static constexpr auto Nihonium = Id::Nh;
		static constexpr auto Flerovium = Id::Fl;
		static constexpr auto Moscovium = Id::Mc;
		static constexpr auto Livermorium = Id::Lv;
		static constexpr auto Tennessine = Id::Ts;
		static constexpr auto Oganesson = Id::Og;

		static constexpr auto H = Id::H;
		static constexpr auto He = Id::He;
		static constexpr auto Li = Id::Li;
		static constexpr auto Be = Id::Be;
		static constexpr auto B = Id::B;
		static constexpr auto C = Id::C;
		static constexpr auto N = Id::N;
		static constexpr auto O = Id::O;
		static constexpr auto F = Id::F;
		static constexpr auto Ne = Id::Ne;
		static constexpr auto Na = Id::Na;
		static constexpr auto Mg = Id::Mg;
		static constexpr auto Al = Id::Al;
		static constexpr auto Si = Id::Si;
		static constexpr auto P = Id::P;
		static constexpr auto S = Id::S;
		static constexpr auto Cl = Id::Cl;
		static constexpr auto Ar = Id::Ar;
		static constexpr auto K = Id::K;
		static constexpr auto Ca = Id::Ca;
		static constexpr auto Sc = Id::Sc;
		static constexpr auto Ti = Id::Ti;
		static constexpr auto V = Id::V;
		static constexpr auto Cr = Id::Cr;
		static constexpr auto Mn = Id::Mn;
		static constexpr auto Fe = Id::Fe;
		static constexpr auto Co = Id::Co;
		static constexpr auto Ni = Id::Ni;
		static constexpr auto Cu = Id::Cu;
		static constexpr auto Zn = Id::Zn;
		static constexpr auto Ga = Id::Ga;
		static constexpr auto Ge = Id::Ge;
		static constexpr auto As = Id::As;
		static constexpr auto Se = Id::Se;
		static constexpr auto Br = Id::Br;
		static constexpr auto Kr = Id::Kr;
		static constexpr auto Rb = Id::Rb;
		static constexpr auto Sr = Id::Sr;
		static constexpr auto Y = Id::Y;
		static constexpr auto Zr = Id::Zr;
		static constexpr auto Nb = Id::Nb;
		static constexpr auto Mo = Id::Mo;
		static constexpr auto Tc = Id::Tc;
		static constexpr auto Ru = Id::Ru;
		static constexpr auto Rh = Id::Rh;
		static constexpr auto Pd = Id::Pd;
		static constexpr auto Ag = Id::Ag;
		static constexpr auto Cd = Id::Cd;
		static constexpr auto In = Id::In;
		static constexpr auto Sn = Id::Sn;
		static constexpr auto Sb = Id::Sb;
		static constexpr auto Te = Id::Te;
		static constexpr auto I = Id::I;
		static constexpr auto Xe = Id::Xe;
		static constexpr auto Cs = Id::Cs;
		static constexpr auto Ba = Id::Ba;
		static constexpr auto La = Id::La;
		static constexpr auto Ce = Id::Ce;
		static constexpr auto Pr = Id::Pr;
		static constexpr auto Nd = Id::Nd;
		static constexpr auto Pm = Id::Pm;
		static constexpr auto Sm = Id::Sm;
		static constexpr auto Eu = Id::Eu;
		static constexpr auto Gd = Id::Gd;
		static constexpr auto Tb = Id::Tb;
		static constexpr auto Dy = Id::Dy;
		static constexpr auto Ho = Id::Ho;
		static constexpr auto Er = Id::Er;
		static constexpr auto Tm = Id::Tm;
		static constexpr auto Yb = Id::Yb;
		static constexpr auto Lu = Id::Lu;
		static constexpr auto Hf = Id::Hf;
		static constexpr auto Ta = Id::Ta;
		static constexpr auto W = Id::W;
		static constexpr auto Re = Id::Re;
		static constexpr auto Os = Id::Os;
		static constexpr auto Ir = Id::Ir;
		static constexpr auto Pt = Id::Pt;
		static constexpr auto Au = Id::Au;
		static constexpr auto Hg = Id::Hg;
		static constexpr auto Tl = Id::Tl;
		static constexpr auto Pb = Id::Pb;
		static constexpr auto Bi = Id::Bi;
		static constexpr auto Po = Id::Po;
		static constexpr auto At = Id::At;
		static constexpr auto Rn = Id::Rn;
		static constexpr auto Fr = Id::Fr;
		static constexpr auto Ra = Id::Ra;
		static constexpr auto Ac = Id::Ac;
		static constexpr auto Th = Id::Th;
		static constexpr auto Pa = Id::Pa;
		static constexpr auto U = Id::U;
		static constexpr auto Np = Id::Np;
		static constexpr auto Pu = Id::Pu;
		static constexpr auto Am = Id::Am;
		static constexpr auto Cm = Id::Cm;
		static constexpr auto Bk = Id::Bk;
		static constexpr auto Cf = Id::Cf;
		static constexpr auto Es = Id::Es;
		static constexpr auto Fm = Id::Fm;
		static constexpr auto Md = Id::Md;
		static constexpr auto No = Id::No;
		static constexpr auto Lr = Id::Lr;
		static constexpr auto Rf = Id::Rf;
		static constexpr auto Db = Id::Db;
		static constexpr auto Sg = Id::Sg;
		static constexpr auto Bh = Id::Bh;
		static constexpr auto Hs = Id::Hs;
		static constexpr auto Mt = Id::Mt;
		static constexpr auto Ds = Id::Ds;
		static constexpr auto Rg = Id::Rg;
		static constexpr auto Cn = Id::Cn;
		static constexpr auto Nh = Id::Nh;
		static constexpr auto Fl = Id::Fl;
		static constexpr auto Mc = Id::Mc;
		static constexpr auto Lv = Id::Lv;
		static constexpr auto Ts = Id::Ts;
		static constexpr auto Og = Id::Og;
#pragma endregion

		Element() = default;


		// ReSharper disable once CppNonExplicitConvertingConstructor
		constexpr Element(const Id id)  // NOLINT(*-explicit-constructor)
		    : m_Id(id)
		{
			/* NO CODE */
		}


		explicit constexpr Element(const int atomicNumber) : m_Id(static_cast<Id>(atomicNumber))
		{
			assert(atomicNumber >= MinAtomicNumber && atomicNumber <= MaxAtomicNumber && "Atomic number out of bound");
		}


		constexpr int AtomicNumber() const
		{
			return static_cast<int>(m_Id);
		}


		constexpr const char* CName() const
		{
			return NameOfElement(AtomicNumber());
		}


		constexpr const char* CSymbol() const
		{
			return SymbolOfElement(AtomicNumber());
		}


		std::string Name() const
		{
			return CName();
		}


		std::string Symbol() const
		{
			return CSymbol();
		}


		constexpr const char* ToCString() const
		{
			return CSymbol();
		}


		std::string ToString() const
		{
			return Symbol();
		}


		friend std::ostream& operator<<(std::ostream& lhs, const Element element)
		{
			return lhs << element.ToString();
		}


		static Element Name2Element(const std::string& identifier)
		{
			for (int i = 0; i < PeriodicTableElementCount; i++)
			{
				if (SecUtility::IsCaseInsensitiveStringEqual(identifier, NameOfElement(i)))
				{
					return Element(i);
				}
			}

			throw std::runtime_error("Could not find symbol \"" + identifier + "\"");
		}


		static Element Symbol2Element(const std::string& identifier)
		{
			for (int i = 0; i < PeriodicTableElementCount; i++)
			{
				if (strcmp(identifier.c_str(), SymbolOfElement(i)) == 0)
				{
					return Element(i);
				}
			}

			throw std::runtime_error("Could not find symbol \"" + identifier + "\"");
		}


		static Element FuzzySymbol2Element(const std::string& identifier)
		{
			for (int i = 1; i < PeriodicTableElementCount; i++)
			{
				if (SecUtility::IsCaseInsensitiveStringEqual(identifier, SymbolOfElement(i)))
				{
					return Element(i);
				}
			}

			throw std::runtime_error("Could not find case insensitive symbol \"" + identifier + "\"");
		}


		static Element ToElement(const std::string& identifier) noexcept
		{
			if (identifier.size() <= 2)
			{
				return Symbol2Element(identifier);
			}

			return Name2Element(identifier);
		}


		CONSTEXPR23 const AtomicElectronConfiguration& MadelungElectronConfiguration() const
		{
			return MadelungElectronConfigurations()[static_cast<int>(m_Id)];
		}


		CONSTEXPR23 const AtomicElectronConfiguration& ElectronConfiguration() const
		{
			return ElectronConfigurations()[static_cast<int>(m_Id)];
		}


		friend std::istream& operator>>(std::istream& is, Element& out_element) noexcept
		{
			std::string symbol;
			is >> symbol;
			out_element = Symbol2Element(symbol);
			return is;
		}


		constexpr int Period() const
		{
			return Period(AtomicNumber());
		}

		constexpr int Group() const
		{
			return PeriodicTableGroupNumber(AtomicNumber());
		}

		constexpr bool IsNeutron() const
		{
			return m_Id == Neutron;
		}

		constexpr bool IsHydrogen() const
		{
			return m_Id == H;
		}

		constexpr bool IsCarbon() const
		{
			return m_Id == C;
		}

		constexpr bool IsNitrogen() const
		{
			return m_Id == N;
		}

		constexpr bool IsOxygen() const
		{
			return m_Id == O;
		}

		constexpr bool IsPhosphorus() const
		{
			return m_Id == P;
		}

		constexpr bool IsSulfur() const
		{
			return m_Id == S;
		}

		constexpr bool IsAlkaliMetal() const
		{
			return Group() == 1 && !IsHydrogen();
		}

		constexpr bool IsAlkaliEarth() const
		{
			return Group() == 2;
		}

		constexpr bool IsHalogen() const
		{
			return Group() == 17;
		}

		constexpr bool IsNobleGas() const
		{
			return Group() == 18;
		}

		constexpr bool IsFromViiiGroup() const
		{
			switch (AtomicNumber())
			{
				case 26:
				case 27:
				case 28:
				case 44:
				case 45:
				case 46:
				case 76:
				case 77:
				case 78:
				case 108:
				case 109:
				case 110:
					return true;
				default:
					return false;
			}
		}

		constexpr bool IsLanthanide() const
		{
			return AtomicNumber() >= 57 && AtomicNumber() <= 71;
		}

		constexpr bool IsActinide() const
		{
			return AtomicNumber() >= 89 && AtomicNumber() <= 103;
		}

		constexpr bool IsTransitionMetal() const
		{
			return IsFromDBlock() && AtomicNumber() != 71 && AtomicNumber() != 103;
		}

		constexpr bool IsBasicMetal() const
		{
			switch (AtomicNumber())
			{
				case 13:
				case 31:
				case 49:
				case 50:
				case 81:
				case 82:
				case 83:
					return true;
				default:
					return AtomicNumber() >= 113 && AtomicNumber() <= 116;
			}
		}

		constexpr bool IsMetalloid() const
		{
			switch (AtomicNumber())
			{
				case 5:
				case 14:
				case 32:
				case 33:
				case 51:
				case 52:
				case 84:
					return true;
				default:
					return false;
			}
		}

		constexpr bool IsNonmetal() const
		{
			switch (AtomicNumber())
			{
				case 6:
				case 7:
				case 8:
				case 15:
				case 16:
				case 34:
					return true;
				default:
					return false;
			}
		}

		// with 103Lr exception
		constexpr int CharacteristicOrbitalAzimuthalQuantumNumber() const
		{
			return m_Id == Id::Lr ? 1 : CharacteristicOrbitalAzimuthalQuantumNumber(AtomicNumber());
		}

		constexpr ElectronicSubShell CharacteristicShell() const
		{
			if (const int l = CharacteristicOrbitalAzimuthalQuantumNumber(); l == 0 || l == 1)
			{
				return {Period(), l};
			}
			else
			{
				return {Period() - (l - 1), l};
			}
		}

		constexpr bool IsFromSBlock() const
		{
			return IsFromSBlock(AtomicNumber());
		}

		constexpr bool IsFromPBlock() const
		{
			return IsFromPBlock(AtomicNumber());
		}

		constexpr bool IsFromDBlock() const
		{
			return IsFromDBlock(AtomicNumber());
		}

		constexpr bool IsFromFBlock() const
		{
			return IsFromFBlock(AtomicNumber());
		}

		constexpr Element NobleGasFromPeriod() const
		{
			return NobleGasFromPeriod(Period());
		}

		static constexpr Element NobleGasFromPeriod(const int period)
		{
			switch (period)
			{
				case 1:
					return Element(2);
				case 2:
					return Element(10);
				case 3:
					return Element(18);
				case 4:
					return Element(36);
				case 5:
					return Element(54);
				case 6:
					return Element(86);
				case 7:
					return Element(118);
				default:
					return Element(0);
			}
		}

		// 103Lr was located at f-block but have outermost shell of 7p.
		constexpr ElectronicSubShell OuterMostShell() const
		{
			return OuterMostShell(AtomicNumber());
		}

		CONSTEXPR23 ElectronicSubShell MadelungOuterMostShell() const
		{
			return MadelungOuterMostShell(AtomicNumber());
		}

		CONSTEXPR23 double MadelungSlaterEffectiveCoreCharge() const
		{
			return MadelungSlaterEffectiveCoreCharges()[AtomicNumber()];
		}

		CONSTEXPR23 double SlaterEffectiveCoreCharge() const
		{
			return SlaterEffectiveCoreCharges()[AtomicNumber()];
		}

		static constexpr double SlaterEffectiveCoreCharge(
		        const Element element,
		        const AtomicElectronConfiguration& electronConfig,
		        const ElectronicSubShell outerMostShell = AtomicElectronConfiguration::MaxElectronShell)
		{
			return SlaterEffectiveCoreCharge(element.AtomicNumber(), electronConfig, outerMostShell);
		}

		// should we include this?
		CONSTEXPR23 double Ghosh2008AtomicRadius() const
		{
			// Equation taken from 10.1016/j.theochem.2008.06.020
			const auto z = AtomicNumber() >= 55 && AtomicNumber() < 103 ? SlaterEffectiveCoreCharge()
			                                                            : MadelungSlaterEffectiveCoreCharge();
			const auto n = MadelungOuterMostShell().SlaterEffectivePrincipalQuantumNumber();
			return n * n / z;
		}

		// should we include this?
		CONSTEXPR23 double Ghosh2002AtomicRadius() const
		{
			// Equation taken from 10.3390/i3020087
			const auto z = AtomicNumber() >= 55 && AtomicNumber() < 103 ? SlaterEffectiveCoreCharge()
			                                                            : MadelungSlaterEffectiveCoreCharge();
			const auto n = MadelungOuterMostShell().SlaterEffectivePrincipalQuantumNumber();
			return n * Period() / z;
		}

		// should we include this?
		CONSTEXPR23 double Ghosh2009AtomicChemicalHardness() const
		{
			// Parameters taken from  10.1002/qua.22202
			// ReSharper disable once CppDFAUnusedValue
			double a = 0;  // so, constexpr requires that?
			// ReSharper disable once CppDFAUnusedValue
			double b = 0;  // so, constexpr requires that?

			switch (Period())
			{
				case 1:
					a = 0.6421;
					b = -2.3061;
					break;
				case 2:
					a = 0.50752;
					b = 0.13044;
					break;
				case 3:
					a = 0.580433;
					b = 0.513833;
					break;
				case 4:
					a = 0.667738;
					b = 0.867338;
					break;
				case 5:
					a = 0.754255;
					b = 0.709427;
					break;
				case 6:
					a = 0.47296;
					b = 0.1196;
					break;
				case 7:
					a = 0.61388;
					b = -0.00546;
					break;
				default:
					throw std::runtime_error("Ghosh2009 parameters is only available for first 7 periods");
			}

			return a * 7.2 / SecUtility::UnitOfMeasurement::BohrRadius2Angstrom(Ghosh2008AtomicRadius()) + b;
		}


	public:
		friend constexpr bool operator<(const Element lhs, const Element rhs)
		{
			return lhs.AtomicNumber() < rhs.AtomicNumber();
		}

		friend constexpr bool operator>(const Element lhs, const Element rhs)
		{
			return lhs.AtomicNumber() > rhs.AtomicNumber();
		}

		friend constexpr bool operator==(const Element lhs, const Element rhs)
		{
			return lhs.AtomicNumber() == rhs.AtomicNumber();
		}

		friend constexpr bool operator!=(const Element lhs, const Element rhs)
		{
			return lhs.AtomicNumber() != rhs.AtomicNumber();
		}


	private:
		Id m_Id = Id::Neutron;


	private:
		static constexpr int Period(const int atomicNumber)
		{
			if (atomicNumber <= 0)
			{
				return 0;
			}
			if (atomicNumber <= 2)
			{
				return 1;
			}
			if (atomicNumber <= 10)
			{
				return 2;
			}
			if (atomicNumber <= 18)
			{
				return 3;
			}
			if (atomicNumber <= 36)
			{
				return 4;
			}
			if (atomicNumber <= 54)
			{
				return 5;
			}
			if (atomicNumber <= 86)
			{
				return 6;
			}
			if (atomicNumber <= 118)
			{
				return 7;
			}

			// For the love of RedHat distribution of GNU 7.3.1 compiler...
#if defined(_MSC_VER)
			__assume(false);
#else
			__builtin_unreachable();
#endif
		}

		static constexpr bool IsFromSBlock(const int atomicNumber)
		{
			return CharacteristicOrbitalAzimuthalQuantumNumber(atomicNumber) == 0;
		}

		// without 103Lr exception
		static constexpr bool IsFromPBlock(const int atomicNumber)
		{
			return CharacteristicOrbitalAzimuthalQuantumNumber(atomicNumber) == 1;
		}

		// without 103Lr exception
		static constexpr bool IsFromDBlock(const int atomicNumber)
		{
			return CharacteristicOrbitalAzimuthalQuantumNumber(atomicNumber) == 2;
		}

		static constexpr bool IsFromFBlock(const int atomicNumber)
		{
			return CharacteristicOrbitalAzimuthalQuantumNumber(atomicNumber) == 3;
		}

		// with 103Lr exception
		static constexpr ElectronicSubShell OuterMostShell(const int atomicNumber)
		{
			return {Period(atomicNumber), IsFromPBlock(atomicNumber) || atomicNumber == 103 ? 1 : 0};
		}

		static constexpr ElectronicSubShell MadelungOuterMostShell(const int atomicNumber)
		{
			return {Period(atomicNumber), IsFromPBlock(atomicNumber) ? 1 : 0};
		}

		static constexpr double SlaterEffectiveCoreCharge(
		        const int atomicNumber,
		        const AtomicElectronConfiguration& electronConfig,
		        ElectronicSubShell outerMostShell = AtomicElectronConfiguration::MaxElectronShell)
		{
			double screeningConst = 0;

			while (electronConfig[outerMostShell] == 0 && outerMostShell.IsValid())
			{
				outerMostShell = ElectronicSubShell(outerMostShell.UnderlyingId() - 1);
			}

			auto currentShell = outerMostShell;

			const auto AdvanceToNextShell = [&currentShell, &electronConfig]
			{
				do
				{
					--currentShell;
					if (!currentShell.IsValid())
					{
						break;
					}
				} while (electronConfig[currentShell] == 0);
			};

			if (outerMostShell.IsPrincipal())
			{
				screeningConst += (electronConfig[outerMostShell] - 1) * 0.35;
				AdvanceToNextShell();
				if (currentShell.PrincipalQuantumNumber() == outerMostShell.PrincipalQuantumNumber())
				{
					screeningConst += electronConfig[currentShell] * 0.35;
					AdvanceToNextShell();
				}
			}
			else if (outerMostShell.IsSharp())
			{
				screeningConst += (electronConfig[outerMostShell] - 1)
				                  * (outerMostShell.PrincipalQuantumNumber() == 1 ? 0.30 : 0.35);
				AdvanceToNextShell();
			}
			else if (outerMostShell.IsDiffuse() || outerMostShell.IsFundamental())
			{
				screeningConst += (electronConfig[outerMostShell] - 1) * 1.00;
				AdvanceToNextShell();
			}

			if (outerMostShell.AzimuthalQuantumNumber().Value() >= 4)
			{
				std::cerr << "Unsupported element for evaluate Slater effective core charge" << std::endl;
				std::terminate();
			}


			while (currentShell.IsValid())
			{
				double screeningContribution = electronConfig[currentShell];

				if (currentShell.PrincipalQuantumNumber() == outerMostShell.PrincipalQuantumNumber() - 1)
				{
					screeningContribution *= currentShell == 5_Diffuse && outerMostShell.IsSharp() ? 0.35 : 0.85;
				}
				else if (currentShell.IsFundamental()
				         && currentShell.PrincipalQuantumNumber() == outerMostShell.PrincipalQuantumNumber() - 2
				         && !outerMostShell.IsPrincipal())
				{
					screeningContribution *= 0.35;
				}

				screeningConst += screeningContribution;

				AdvanceToNextShell();
			}

			return atomicNumber - screeningConst;
		}


		static constexpr const char* NameOfElement(const int atomicNumber)
		{
			return PeriodicTableNames[atomicNumber];
		}
		static constexpr std::array<const char*, PeriodicTableElementCount> PeriodicTableNames = {
		        "Neutron",    "Hydrogen",     "Helium",        "Lithium",     "Beryllium",   "Boron",
		        "Carbon",     "Nitrogen",     "Oxygen",        "Fluorine",    "Neon",        "Sodium",
		        "Magnesium",  "Aluminum",     "Silicon",       "Phosphorus",  "Sulfur",      "Chlorine",
		        "Argon",      "Potassium",    "Calcium",       "Scandium",    "Titanium",    "Vanadium",
		        "Chromium",   "Manganese",    "Iron",          "Cobalt",      "Nickel",      "Copper",
		        "Zinc",       "Gallium",      "Germanium",     "Arsenic",     "Selenium",    "Bromine",
		        "Krypton",    "Rubidium",     "Strontium",     "Yttrium",     "Zirconium",   "Niobium",
		        "Molybdenum", "Technetium",   "Ruthenium",     "Rhodium",     "Palladium",   "Silver",
		        "Cadmium",    "Indium",       "Tin",           "Antimony",    "Tellurium",   "Iodine",
		        "Xenon",      "Cesium",       "Barium",        "Lanthanum",   "Cerium",      "Praseodymium",
		        "Neodymium",  "Promethium",   "Samarium",      "Europium",    "Gadolinium",  "Terbium",
		        "Dysprosium", "Holmium",      "Erbium",        "Thulium",     "Ytterbium",   "Lutetium",
		        "Hafnium",    "Tantalum",     "Tungsten",      "Rhenium",     "Osmium",      "Iridium",
		        "Platinum",   "Gold",         "Mercury",       "Thallium",    "Lead",        "Bismuth",
		        "Polonium",   "Astatine",     "Radon",         "Francium",    "Radium",      "Actinium",
		        "Thorium",    "Protactinium", "Uranium",       "Neptunium",   "Plutonium",   "Americium",
		        "Curium",     "Berkelium",    "Californium",   "Einsteinium", "Fermium",     "Mendelevium",
		        "Nobelium",   "Lawrencium",   "Rutherfordium", "Dubnium",     "Seaborgium",  "Bohrium",
		        "Hassium",    "Meitnerium",   "Darmstadtium",  "Roentgenium", "Copernicium", "Nihonium",
		        "Flerovium",  "Moscovium",    "Livermorium",   "Tennessine",  "Oganesson"};


		static constexpr const char* SymbolOfElement(const int atomicNumber)
		{
			return PeriodicTableSymbols[atomicNumber];
		}
		static constexpr std::array<const char*, PeriodicTableElementCount> PeriodicTableSymbols = {
		        "n",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",
		        "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
		        "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
		        "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
		        "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
		        "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
		        "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};


		static constexpr int PeriodicTableGroupNumber(const int atomicNumber)
		{
			return PeriodicTableGroupNumbers[atomicNumber];
		}
		static constexpr std::array<int, PeriodicTableElementCount> PeriodicTableGroupNumbers = {
		        0, /* Neutron */
		        1,  18, 1,  2,  13, 14, 15, 16, 17, 18, 1,  2,  13, 14, 15, 16, 17, 18, 1,  2,  3,  4,  5,  6,
		        7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
		        13, 14, 15, 16, 17, 18, 1,  2,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 4,
		        5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 1,  2,  -2, -2, -2, -2, -2, -2, -2, -2,
		        -2, -2, -2, -2, -2, -2, -2, 4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18};

		// without 103Lr exception
		static constexpr int CharacteristicOrbitalAzimuthalQuantumNumber(const int atomicNumber)
		{
			return CharacteristicOrbitalAzimuthalQuantumNumbers[atomicNumber];
		}
		static constexpr std::array<int, PeriodicTableElementCount> CharacteristicOrbitalAzimuthalQuantumNumbers = {
		        0, /* Neutron */
		        0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
		        1, 1, 1, 1, 1, 1, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 3, 3, 3, 3,
		        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 3, 3,
		        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1};


		static constexpr std::array<int, 7> NobleGasAtomicNumbers = {2, 10, 18, 36, 54, 86, 118};

		static CONSTEXPR23 const std::array<AtomicElectronConfiguration, PeriodicTableElementCount>&
		MadelungElectronConfigurations()
		{
			static const std::array<AtomicElectronConfiguration, PeriodicTableElementCount> Data = []
			{
				std::array<AtomicElectronConfiguration, PeriodicTableElementCount> arr{};

				for (int z = 0; z < PeriodicTableElementCount; z++)
				{
					arr[z] = AtomicElectronConfiguration(z);
				}

				return arr;
			}();

			return Data;
		}

		static CONSTEXPR23 const std::array<AtomicElectronConfiguration, PeriodicTableElementCount>&
		ElectronConfigurations()
		{
			static const std::array<AtomicElectronConfiguration, PeriodicTableElementCount> Data = []
			{
				auto arr = MadelungElectronConfigurations();
				// 3d transition metal exceptions
				arr[24].MoveElectron(4_Sharp, 3_Diffuse);
				arr[29].MoveElectron(4_Sharp, 3_Diffuse);
				// 4d transition metal exceptions
				arr[41].MoveElectron(5_Sharp, 4_Diffuse);
				arr[42].MoveElectron(5_Sharp, 4_Diffuse);
				arr[44].MoveElectron(5_Sharp, 4_Diffuse);
				arr[45].MoveElectron(5_Sharp, 4_Diffuse);
				arr[46].MoveElectrons(5_Sharp, 4_Diffuse, 2);
				arr[47].MoveElectron(5_Sharp, 4_Diffuse);
				// 5d transition metal exceptions
				arr[78].MoveElectron(6_Sharp, 5_Diffuse);
				arr[79].MoveElectron(6_Sharp, 5_Diffuse);
				// Lanthanide exceptions
				arr[57].MoveElectron(4_Fundamental, 5_Diffuse);
				arr[58].MoveElectron(4_Fundamental, 5_Diffuse);
				arr[64].MoveElectron(4_Fundamental, 5_Diffuse);
				// Actinide exceptions
				arr[89].MoveElectron(5_Fundamental, 6_Diffuse);
				arr[90].MoveElectrons(5_Fundamental, 6_Diffuse, 2);
				arr[91].MoveElectron(5_Fundamental, 6_Diffuse);
				arr[92].MoveElectron(5_Fundamental, 6_Diffuse);
				arr[93].MoveElectron(5_Fundamental, 6_Diffuse);
				arr[103].MoveElectron(6_Diffuse, 7_Principal);
				return arr;
			}();

			return Data;
		}

		static CONSTEXPR23 const std::array<double, PeriodicTableElementCount>& MadelungSlaterEffectiveCoreCharges()
		{
			static const std::array<double, PeriodicTableElementCount> Data = []
			{
				std::array<double, PeriodicTableElementCount> charges{};
				for (int z = 1; z < PeriodicTableElementCount; z++)  // Yes, index starts from 1, not 0
				{
					charges[z] = SlaterEffectiveCoreCharge(
					        z, MadelungElectronConfigurations()[z], MadelungOuterMostShell(z));
				}
				return charges;
			}();

			return Data;
		}

		static CONSTEXPR23 const std::array<double, PeriodicTableElementCount>& SlaterEffectiveCoreCharges()
		{
			static const std::array<double, PeriodicTableElementCount> Data = []
			{
				std::array<double, PeriodicTableElementCount> charges{};

				for (int z = 1; z < PeriodicTableElementCount; z++)  // Yes, index starts from 1, not 0
				{
					charges[z] = SlaterEffectiveCoreCharge(z, ElectronConfigurations()[z], OuterMostShell(z));
				}

				return charges;
			}();

			return Data;
		}
	};
}  // namespace SecChem

template <>
struct std::hash<SecChem::Element>
{
	std::size_t operator()(const SecChem::Element& element) const noexcept
	{
		return std::hash<int>{}(element.AtomicNumber());
	}
};
