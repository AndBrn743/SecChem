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
	class AbstractMolecule;

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
	/// AbstractMolecule must use contiguous Atom storage
	template <typename Derived>
	class AbstractMolecule
	{
	public:
		std::size_t AtomCount() const noexcept
		{
			return static_cast<const Derived*>(this)->AtomCount_Impl();
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
#if __cpp_lib_to_address >= 201711L
			const auto first = std::to_address(cbegin());
#else
			const auto first = static_cast<const Atom*>(cbegin().operator->());
#endif
			static_assert(std::is_same_v<decltype(first), const Atom* const>);
			const auto last = first + AtomCount();
			const auto ptr = std::addressof(atom);

			if (ptr >= first && ptr < last)
			{
				return static_cast<std::size_t>(ptr - first);
			}

			if (const auto it = std::find(cbegin(), cend(), atom); it != cend())
			{
				return std::distance(cbegin(), it);
			}

			throw std::out_of_range("Atom not found");
		}

		template <typename CountingFunction>
		double CoordinationNumberOfAtomAtIndex(const std::size_t atomIndex, CountingFunction countingFunction) const
		{
			double cn = 0;

			for (std::size_t i = 0; i < AtomCount(); i++)
			{
				if (i != atomIndex)
				{
					cn += countingFunction((*this)[i], (*this)[atomIndex]);
				}
			}

			return cn;
		}

		template <typename CountingFunction>
		double CoordinationNumberOf(const Atom& atom, CountingFunction countingFunction) const
		{
			return CoordinationNumberOfAtomAtIndex(IndexOf(atom), countingFunction);
		}

		template <typename CountingFunction>
		std::vector<double> GenerateCoordinationNumbers(CountingFunction countingFunction) const
		{
			std::vector<double> coordinationNumbers(AtomCount());

			for (std::size_t i = 0; i < AtomCount(); i++)
			{
				for (std::size_t j = 0; j < i; j++)
				{
					const auto cn = countingFunction((*this)[i], (*this)[j]);
					coordinationNumbers[i] += cn;
					coordinationNumbers[j] += cn;
				}
			}

			return coordinationNumbers;
		}

		bool Contains(const Atom& atom) const noexcept
		{
#if __cpp_lib_to_address >= 201711L
			const auto first = std::to_address(cbegin());
#else
			const auto first = static_cast<const Atom*>(cbegin().operator->());
#endif
			static_assert(std::is_same_v<decltype(first), const Atom* const>);
			const auto last = first + AtomCount();
			const auto ptr = &atom;

			if (ptr >= first && ptr < last)
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

		/* CRTP VIRTUAL */ std::size_t AtomCount_Impl() const noexcept
		{
			return std::distance(cbegin(), cend());
		}

		/* CRTP PURE VIRTUAL */ auto begin_Impl() noexcept = delete;
		/* CRTP PURE VIRTUAL */ auto end_Impl() noexcept = delete;
		/* CRTP PURE VIRTUAL */ auto cbegin_Impl() const noexcept = delete;
		/* CRTP PURE VIRTUAL */ auto cend_Impl() const noexcept = delete;
	};

	class Molecule : public AbstractMolecule<Molecule>
	{
		friend AbstractMolecule;
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
			return m_Atoms.begin();
		}

		/* CRTP OVERRIDE */ auto end_Impl() noexcept
		{
			return m_Atoms.end();
		}

		/* CRTP OVERRIDE */ auto cbegin_Impl() const noexcept
		{
			return m_Atoms.cbegin();
		}

		/* CRTP OVERRIDE */ auto cend_Impl() const noexcept
		{
			return m_Atoms.cend();
		}

		std::vector<Atom> m_Atoms;
	};

	class SharedMolecule : public AbstractMolecule<SharedMolecule>
	{
		friend AbstractMolecule;
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
			return m_AtomsPtr->begin();
		}

		// ReSharper disable once CppMemberFunctionMayBeConst
		/* CRTP OVERRIDE */ auto end_Impl() noexcept
		{
			return m_AtomsPtr->end();
		}

		/* CRTP OVERRIDE */ auto cbegin_Impl() const noexcept
		{
			return m_AtomsPtr->cbegin();
		}

		/* CRTP OVERRIDE */ auto cend_Impl() const noexcept
		{
			return m_AtomsPtr->cend();
		}

		std::shared_ptr<std::vector<Atom>> m_AtomsPtr;
	};

	inline Molecule::Molecule(const SharedMolecule& sharedMolecule)
	    : m_Atoms(sharedMolecule.cbegin(), sharedMolecule.cend())
	{
		/* NO CODE */
	}
}  // namespace SecChem
