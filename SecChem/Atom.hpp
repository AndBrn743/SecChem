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

	constexpr AtomTag operator|(const AtomTag lhs, const AtomTag rhs) noexcept
	{
		return static_cast<AtomTag>(static_cast<std::uint16_t>(lhs) | static_cast<std::uint16_t>(rhs));
	}

	constexpr void operator|=(AtomTag& lhs, const AtomTag rhs) noexcept
	{
		lhs = static_cast<AtomTag>(static_cast<std::uint16_t>(lhs) | static_cast<std::uint16_t>(rhs));
	}

	constexpr AtomTag operator&(const AtomTag lhs, const AtomTag rhs) noexcept
	{
		return static_cast<AtomTag>(static_cast<std::uint16_t>(lhs) & static_cast<std::uint16_t>(rhs));
	}

	constexpr AtomTag AllAtomTags = AtomTag::GaussianFiniteNuclear | AtomTag::ThomasFermiFiniteNuclear | AtomTag::Frozen
	                                | AtomTag::Link | AtomTag::Buffer;

	constexpr AtomTag operator~(AtomTag tag) noexcept
	{
		return static_cast<AtomTag>(static_cast<std::uint16_t>(AllAtomTags) ^ static_cast<std::uint16_t>(tag));
	}

	constexpr bool IsValid(const AtomTag tag) noexcept
	{
		if ((tag & ~AllAtomTags) != AtomTag::None)
		{
			return false;
		}

		if (static_cast<bool>(tag & AtomTag::GaussianFiniteNuclear)
		    && static_cast<bool>(tag & AtomTag::ThomasFermiFiniteNuclear))
		{
			return false;
		}

		return true;
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
		Atom(const Element element, const Eigen::Vector3d& position, const AtomTag tags = {})
		    : Atom(element, position, tags, element.Mass(), element.NuclearRadius())
		{
			/* NO CODE */
		}

		static Atom AtomWithMass(const Element element,
		                         const Eigen::Vector3d& position,
		                         const double mass,
		                         const AtomTag tags = {})
		{
			return Atom{element, position, tags, mass, element.NuclearRadius()};
		}

		static Atom AtomWithNuclearRadius(const Element element,
		                                  const Eigen::Vector3d& position,
		                                  const double nuclearRadius,
		                                  const AtomTag tag = {})
		{
			return Atom{element, position, tag, element.Mass(), nuclearRadius};
		}

		static Atom AtomWithMassAndNuclearRadius(const Element element,
		                                         const Eigen::Vector3d& position,
		                                         const double mass,
		                                         const double nuclearRadius,
		                                         const AtomTag tag = {})
		{
			return Atom{element, position, tag, mass, nuclearRadius};
		}

		// GCC require it
		// ReSharper disable once CppRedundantElaboratedTypeSpecifier
		constexpr class Element Element() const noexcept
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

		void SetPosition(const Eigen::Vector3d& newPosition) noexcept
		{
			m_Position = newPosition;
		}

		void Translate(const Eigen::Vector3d& delta) noexcept
		{
			m_Position += delta;
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
			return static_cast<bool>(m_Tags & AtomTag::Frozen);
		}

		double DistanceTo(const Atom other) const noexcept
		{
			return (m_Position - other.m_Position).norm();
		}

		void RemoveTags(const AtomTag tags) noexcept
		{
			m_Tags = m_Tags & ~tags;
		}

		void AddTags(const AtomTag tags)
		{
			if (!IsValid(tags))
			{
				throw std::invalid_argument("Atom tag " + std::to_string(static_cast<int>(m_Tags)) + " is invalid");
			}

			if (IsWithThomasFermiFiniteNuclear() && static_cast<bool>(tags & AtomTag::GaussianFiniteNuclear))
			{
				RemoveTags(AtomTag::ThomasFermiFiniteNuclear);
			}
			else if (IsWithGaussianFiniteNuclear() && static_cast<bool>(tags & AtomTag::ThomasFermiFiniteNuclear))
			{
				RemoveTags(AtomTag::GaussianFiniteNuclear);
			}

			m_Tags = m_Tags | tags;
		}


	private:
		Atom(const class Element element,
		     const Eigen::Vector3d& position,
		     const AtomTag tags,
		     const double mass,
		     const double nuclearRadius)
		    : m_Element(element), m_Tags(tags), m_Position(position), m_Mass(mass), m_NuclearRadius(nuclearRadius)
		{
			if (!IsValid(tags))
			{
				throw std::invalid_argument("Atom tag " + std::to_string(static_cast<int>(tags)) + " is invalid");
			}
		}

		bool EqualsTo_Impl(const Atom& other, const double tolerance) const noexcept
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
