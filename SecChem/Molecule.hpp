// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Andy Brown

#pragma once

#include <cassert>
#include <numeric>
#include <utility>
#include <vector>

#include "Utility/IEquatableWithTolerance.hpp"

#include "Atom.hpp"

namespace SecChem
{
	template <typename Derived>
	class IMolecule;

	class Molecule;

	class SharedMolecule;

	template <typename>
	class Builder;
}  // namespace SecChem

template <>
struct SecUtility::Traits<SecChem::Molecule>
{
	static constexpr auto DefaultEqualityComparisonTolerance = 1e-12;
};

namespace SecChem
{
	template <typename Derived>
	class IMolecule
	{
	public:
		std::size_t AtomCount() const noexcept
		{
			return std::distance(cbegin(), cend());
		}

		std::size_t TotalNuclearCharge() const noexcept
		{
			return std::accumulate(cbegin(),
			                       cend(),
			                       std::size_t{0},
			                       [](const std::size_t acc, const Atom& atom)
			                       { return acc + static_cast<std::size_t>(atom.AtomicNumber()); });
		}

		double Mass() const noexcept
		{
			return std::accumulate(
			        cbegin(), cend(), double{0}, [](const double acc, const Atom& atom) { return acc + atom.Mass(); });
		}

		double NuclearRepulsionEnergy() const noexcept
		{
			double e = 0;

			for (auto it1 = cbegin(); it1 != cend(); ++it1)
			{
				const auto& atom1 = *it1;

				for (auto it2 = cbegin(); it2 != it1; ++it2)
				{
					const auto& atom2 = *it2;

					e += atom1.Element.AtomicNumber() * atom2.Element.AtomicNumber() / atom1.DistanceTo(atom2);
				}
			}

			return e;
		}

		const Atom& operator[](const std::size_t i) const noexcept
		{
			assert(i < AtomCount());
			return begin()[i];
		}

		Atom& operator[](const std::size_t i) noexcept
		{
			assert(i < AtomCount());
			return begin()[i];
		}

		std::size_t IndexOf(const Atom& atom) const
		{
			if (&atom >= cbegin() && &atom < cend())
			{
				return std::distance(static_cast<const Atom*>(cbegin()), &atom);
			}

			if (const auto it = std::find(cbegin(), cend(), atom); it != cend())
			{
				return std::distance(cbegin(), it);
			}

			throw std::out_of_range("Atom not found");
		}

		bool Contains(const Atom& atom) const noexcept
		{
			if (&atom >= cbegin() && &atom < cend())
			{
				return true;
			}

			return std::find(cbegin(), cend(), atom) != cend();
		}

		auto begin() noexcept
		{
			return AsDerived().begin_Impl();
		}

		auto end() noexcept
		{
			return AsDerived().end_Impl();
		}

		auto begin() const noexcept
		{
			return AsDerived().cbegin_Impl();
		}

		auto end() const noexcept
		{
			return AsDerived().cend_Impl();
		}

		auto cbegin() const noexcept
		{
			return AsDerived().cbegin_Impl();
		}

		auto cend() const noexcept
		{
			return AsDerived().cend_Impl();
		}

	private:
		constexpr const Derived& AsDerived() const noexcept
		{
			return *static_cast<const Derived*>(this);
		}

		constexpr Derived& AsDerived() noexcept
		{
			return *static_cast<Derived*>(this);
		}

		/* CRTP PURE VIRTUAL */ auto begin_Impl() noexcept = delete;
		/* CRTP PURE VIRTUAL */ auto end_Impl() noexcept = delete;
		/* CRTP PURE VIRTUAL */ auto cbegin_Impl() const noexcept = delete;
		/* CRTP PURE VIRTUAL */ auto cend_Impl() const noexcept = delete;
	};

	class Molecule : public IMolecule<Molecule>
	{
		friend IMolecule;
		friend SharedMolecule;

	public:
		Molecule() = default;

		template <typename InputIterator>
		Molecule(InputIterator atomsBegin, const InputIterator atomsEnd) : m_Atoms(atomsBegin, atomsEnd)
		{
			/* NO CODE */
		}

		Molecule(const std::initializer_list<Atom> atoms) : m_Atoms(atoms)
		{
			/* NO CODE */
		}

		explicit Molecule(std::vector<Atom> atoms) : m_Atoms(std::move(atoms))
		{
			/* NO CODE */
		}

		explicit Molecule(const SharedMolecule& sharedMolecule);


	private:
		/* CRTP OVERRIDE */ auto begin_Impl() noexcept
		{
			return m_Atoms.data();
		}

		/* CRTP OVERRIDE */ auto end_Impl() noexcept
		{
			return m_Atoms.data() + m_Atoms.size();
		}

		/* CRTP OVERRIDE */ auto cbegin_Impl() const noexcept
		{
			return m_Atoms.data();
		}

		/* CRTP OVERRIDE */ auto cend_Impl() const noexcept
		{
			return m_Atoms.data() + m_Atoms.size();
		}

		std::vector<Atom> m_Atoms;
	};

	class SharedMolecule : public IMolecule<SharedMolecule>
	{
		friend IMolecule;
		friend Molecule;

	public:
		template <typename InputIterator>
		SharedMolecule(InputIterator atomsBegin, const InputIterator atomsEnd)
		    : m_AtomsPtr(std::make_shared<std::vector<Atom>>(atomsBegin, atomsEnd))
		{
			/* NO CODE */
		}

		SharedMolecule(const std::initializer_list<Atom> atoms) : m_AtomsPtr(std::make_shared<std::vector<Atom>>(atoms))
		{
			/* NO CODE */
		}

		explicit SharedMolecule(std::vector<Atom> atoms)
		    : m_AtomsPtr(std::make_shared<std::vector<Atom>>(std::move(atoms)))
		{
			/* NO CODE */
		}

		explicit SharedMolecule(const Molecule& molecule)
		    : m_AtomsPtr(std::make_shared<std::vector<Atom>>(molecule.m_Atoms))
		{
			/* NO CODE */
		}

		explicit SharedMolecule(Molecule&& molecule)
		    : m_AtomsPtr(std::make_shared<std::vector<Atom>>(std::move(molecule).m_Atoms))
		{
			/* NO CODE */
		}

	private:
		// ReSharper disable once CppMemberFunctionMayBeConst
		/* CRTP OVERRIDE */ auto begin_Impl() noexcept
		{
			return m_AtomsPtr->data();
		}

		// ReSharper disable once CppMemberFunctionMayBeConst
		/* CRTP OVERRIDE */ auto end_Impl() noexcept
		{
			return m_AtomsPtr->data() + m_AtomsPtr->size();
		}

		/* CRTP OVERRIDE */ auto cbegin_Impl() const noexcept
		{
			return static_cast<const Atom*>(m_AtomsPtr->data());
		}

		/* CRTP OVERRIDE */ auto cend_Impl() const noexcept
		{
			return static_cast<const Atom*>(m_AtomsPtr->data() + m_AtomsPtr->size());
		}

		std::shared_ptr<std::vector<Atom>> m_AtomsPtr;
	};

	inline Molecule::Molecule(const SharedMolecule& sharedMolecule)
	    : m_Atoms(sharedMolecule.cbegin(), sharedMolecule.cend())
	{
		/* NO CODE */
	}
}  // namespace SecChem
