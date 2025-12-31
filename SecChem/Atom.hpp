//
// Created by Andy on 11/22/2025.
//

#pragma once

#include <Eigen/Dense>
#include <cinttypes>

#include "Utility/IEquatableWithTolerance.hpp"

#include "Element.hpp"

namespace SecChem
{
	enum class AtomTag : std::uint16_t
	{
		None = 0,
		GaussianFiniteNuclear = 1 << 0,
		ThomasFermiFiniteNuclear = 1 << 1,
		Frozen = 1 << 2,
		Link = 1 << 3,
		Buffer = 1 << 4,
	};

	constexpr AtomTag operator~(const AtomTag tag) noexcept
	{
		return static_cast<AtomTag>(~static_cast<std::uint16_t>(tag));
	}

	constexpr AtomTag operator|(const AtomTag lhs, const AtomTag rhs) noexcept
	{
		return static_cast<AtomTag>(static_cast<std::uint16_t>(lhs) | static_cast<std::uint16_t>(rhs));
	}

	constexpr AtomTag operator&(const AtomTag lhs, const AtomTag rhs) noexcept
	{
		return static_cast<AtomTag>(static_cast<std::uint16_t>(lhs) & static_cast<std::uint16_t>(rhs));
	}

	class Atom;
}  // namespace SecChem

template <>
struct SecUtility::Traits<SecChem::Atom>
{
	static constexpr auto DefaultEqualityComparisonTolerance = 1e-12;
};

namespace SecChem
{
	class Atom : public SecUtility::IEquatableWithTolerance<Atom>
	{
		friend IEquatableWithTolerance;

	public:
		Atom(const Element element, const Eigen::Vector3d& position, const AtomTag tag = {}) noexcept
		    : m_Element(element), m_Position(position), m_Tags(tag), m_Mass(element.Mass()),
		      m_NuclearRadius(element.NuclearRadius())
		{
			/* NO CODE */
		}

		static Atom AtomWithMass(const Element element,
		                         const Eigen::Vector3d& position,
		                         const double mass,
		                         const AtomTag tag = {}) noexcept
		{
			return Atom{element, position, tag, mass, element.NuclearRadius()};
		}

		static Atom AtomWithNuclearRadius(const Element element,
		                                  const Eigen::Vector3d& position,
		                                  const double nuclearRadius,
		                                  const AtomTag tag = {}) noexcept
		{
			return Atom{element, position, tag, element.Mass(), nuclearRadius};
		}

		static Atom AtomWithMassAndNuclearRadius(const Element element,
		                                         const Eigen::Vector3d& position,
		                                         const double mass,
		                                         const double nuclearRadius,
		                                         const AtomTag tag = {}) noexcept
		{
			return Atom{element, position, tag, mass, nuclearRadius};
		}

		constexpr Element Element() const noexcept
		{
			return m_Element;
		}

		constexpr auto AtomicNumber() const noexcept
		{
			return m_Element.AtomicNumber();
		}

		constexpr const Eigen::Vector3d& Position() const noexcept
		{
			return m_Position;
		}

		constexpr double Mass() const noexcept
		{
			return m_Mass;
		}

		constexpr double NuclearRadius() const noexcept
		{
			return m_NuclearRadius;
		}

		constexpr bool IsWithGaussianFiniteNuclear() const noexcept
		{
			return static_cast<bool>(m_Tags & AtomTag::GaussianFiniteNuclear);
		}

		constexpr bool IsWithThomasFermiFiniteNuclear() const noexcept
		{
			return static_cast<bool>(m_Tags & AtomTag::ThomasFermiFiniteNuclear);
		}

		constexpr bool IsWithFiniteNuclear() const noexcept
		{
			return static_cast<bool>(m_Tags & (AtomTag::GaussianFiniteNuclear | AtomTag::ThomasFermiFiniteNuclear));
		}

		constexpr bool IsFrozen() const noexcept
		{
			return static_cast<int>(m_Tags) & static_cast<int>(AtomTag::Frozen);
		}

		bool DistanceTo(const Atom other) const noexcept
		{
			return (m_Position - other.m_Position).norm();
		}


	private:
		Atom(const class Element element,
		     const Eigen::Vector3d& position,
		     const AtomTag tag,
		     const double mass,
		     const double nuclearRadius) noexcept
		    : m_Element(element), m_Position(position), m_Tags(tag), m_Mass(mass), m_NuclearRadius(nuclearRadius)
		{
			assert(!(static_cast<bool>(tag & AtomTag::GaussianFiniteNuclear)
			         && static_cast<bool>(tag & AtomTag::ThomasFermiFiniteNuclear)));
		}

		bool EqualsTo_Impl(const Atom other, const double tolerance) const noexcept
		{
			return m_Element == other.m_Element && m_Tags == other.m_Tags
			       && std::abs(m_Mass - other.Mass()) <= tolerance
			       && std::abs(m_NuclearRadius - other.NuclearRadius()) <= tolerance
			       && (m_Position - other.m_Position).norm() <= tolerance;
		}

		class Element m_Element = Element::Neutron;
		AtomTag m_Tags = AtomTag::None;
		Eigen::Vector3d m_Position = {};
		double m_Mass = 0;
		double m_NuclearRadius = 0;
	};
}  // namespace SecChem
