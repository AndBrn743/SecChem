// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Andy Brown

#pragma once
#if !defined(SECCHEM_GAUSSIAN_BASIS_SET_INTERNAL)
#error Do not include internal header files directly
#endif

#include <sstream>
#include <string>

#include "../../Molecule.hpp"


namespace SecChem::Gaussian
{
	class MolecularBasisSet
	{
		friend Builder<MolecularBasisSet>;

		using ElementaryBasis = const std::vector<BasisSet::Gaussian::AngularMomentumBlock>;

		// The raw pointer must points to the data member of SharedBasisSetLibrary. SharedBasisSetLibrary will
		// provide stable iterator, pointer, and reference and prevent dangling
		using ElementaryBasisPtr = const ElementaryBasis*;

	public:
		const BasisSet::Gaussian::SharedBasisSetLibrary& SharedBasisSetLibrary() const noexcept
		{
			return m_Library;
		}

		const SharedMolecule& Molecule() const noexcept
		{
			return m_Molecule;
		}

		const ElementaryBasis& ElementaryBasisAt(const std::size_t index) const
		{
			if (index >= m_BasisAssignments.size())
			{
				throw std::out_of_range("Basis assignment index out of range");
			}

			return *m_BasisAssignments[index];
		}

		const ElementaryBasis& ElementaryBasisOf(const Atom& atom) const
		{
			return *m_BasisAssignments[m_Molecule.IndexOf(atom)];
		}

		std::size_t UniqueElementaryBasisCount() const noexcept
		{
			return m_SubShellSegmentationTableOfEachElementaryBasis.size();
		}

		// this method DOES cache
		Eigen::Index PrimitiveSphericalOrbitalCountOf(const Atom& atom) const
		{
			const auto index = Molecule().IndexOf(atom);
			return m_SegmentationTableOfPrimitiveSphericalOrbitals[index + 1]
			       - m_SegmentationTableOfPrimitiveSphericalOrbitals[index];
		}

		// this method does NOT cache
		Eigen::Index PrimitiveCartesianOrbitalCountOf(const Atom& atom) const
		{
			return CountOfSomeKindOfOrbitalOf(atom,
			                                  [](const BasisSet::Gaussian::AngularMomentumBlock& amb)
			                                  { return amb.PrimitiveCartesianOrbitalCount(); });
		}

		// this method DOES cache
		Eigen::Index ContractedSphericalOrbitalCountOf(const Atom& atom) const
		{
			const auto index = Molecule().IndexOf(atom);
			return m_SegmentationTableOfContractedSphericalOrbitals[index + 1]
			       - m_SegmentationTableOfContractedSphericalOrbitals[index];
		}

		// this method does NOT cache
		Eigen::Index ContractedCartesianOrbitalCountOf(const Atom& atom) const
		{
			return CountOfSomeKindOfOrbitalOf(atom,
			                                  [](const BasisSet::Gaussian::AngularMomentumBlock& amb)
			                                  { return amb.ContractedCartesianOrbitalCount(); });
		}

		// this method DOES cache
		Eigen::Index PrimitiveSphericalOrbitalCount() const noexcept
		{
			return m_SegmentationTableOfPrimitiveSphericalOrbitals.back();
		}

		// this method does NOT cache
		Eigen::Index PrimitiveCartesianOrbitalCountOf() const noexcept
		{
			return CountOfSomeKindOfOrbital([](const BasisSet::Gaussian::AngularMomentumBlock& amb)
			                                { return amb.PrimitiveCartesianOrbitalCount(); });
		}

		// this method DOES cache
		Eigen::Index ContractedSphericalOrbitalCount() const noexcept
		{
			return m_SegmentationTableOfContractedSphericalOrbitals.back();
		}

		// this method does NOT cache
		Eigen::Index ContractedCartesianOrbitalCountOf() const noexcept
		{
			return CountOfSomeKindOfOrbital([](const BasisSet::Gaussian::AngularMomentumBlock& amb)
			                                { return amb.ContractedCartesianOrbitalCount(); });
		}

		// this method does NOT cache
		Eigen::Index PrimitiveSubShellCountOf(const Atom& atom) const
		{
			const auto& basis = ElementaryBasisOf(atom);
			return std::accumulate(basis.cbegin(),
			                       basis.cend(),
			                       Eigen::Index{0},
			                       [](const Eigen::Index acc, const BasisSet::Gaussian::AngularMomentumBlock& amb)
			                       { return acc + amb.PrimitiveShellCount(); });
		}

		// this method does NOT cache
		Eigen::Index ContractedSubShellCountOf(const Atom& atom) const
		{
			const auto& basis = ElementaryBasisOf(atom);
			return std::accumulate(basis.cbegin(),
			                       basis.cend(),
			                       Eigen::Index{0},
			                       [](const Eigen::Index acc, const BasisSet::Gaussian::AngularMomentumBlock& amb)
			                       { return acc + amb.ContractedShellCount(); });
		}

		Eigen::Index AtomIndexFromContractedSphericalOrbital(const Eigen::Index orbitalIndex) const noexcept
		{
			assert(orbitalIndex >= 0 && orbitalIndex <= m_SegmentationTableOfContractedSphericalOrbitals.back());

			return std::distance(m_SegmentationTableOfContractedSphericalOrbitals.cbegin(),
			                     std::upper_bound(m_SegmentationTableOfContractedSphericalOrbitals.cbegin(),
			                                      m_SegmentationTableOfContractedSphericalOrbitals.cend(),
			                                      orbitalIndex))
			       - 1;
		}

		const Atom& AtomOfContractedSphericalOrbital(const Eigen::Index orbitalIndex) const noexcept
		{
			return Molecule()[AtomIndexFromContractedSphericalOrbital(orbitalIndex)];
		}

		Eigen::Index IndexOfAtomOfPrimitiveSphericalOrbital(const Eigen::Index orbitalIndex) const noexcept
		{
			assert(orbitalIndex >= 0 && orbitalIndex <= m_SegmentationTableOfPrimitiveSphericalOrbitals.back());

			return std::distance(m_SegmentationTableOfPrimitiveSphericalOrbitals.cbegin(),
			                     std::upper_bound(m_SegmentationTableOfPrimitiveSphericalOrbitals.cbegin(),
			                                      m_SegmentationTableOfPrimitiveSphericalOrbitals.cend(),
			                                      orbitalIndex));
		}

		const Atom& AtomOfPrimitiveSphericalOrbital(const Eigen::Index orbitalIndex) const noexcept
		{
			return Molecule()[IndexOfAtomOfPrimitiveSphericalOrbital(orbitalIndex)];
		}

		ElectronicSubShell AtomicSubShellFromContractedSphericalOrbital(const Eigen::Index orbitalIndex) const noexcept
		{
			const auto atomIndex = AtomIndexFromContractedSphericalOrbital(orbitalIndex);
			return AtomicSubShellFromContractedSphericalOrbital(orbitalIndex, atomIndex);
		}

		std::pair<const Atom&, ElectronicSubShell> AtomAndSubShellFromContractedSphericalOrbital(
		        const Eigen::Index orbitalIndex) const noexcept
		{
			const auto atomIndex = AtomIndexFromContractedSphericalOrbital(orbitalIndex);
			return std::pair<const Atom&, ElectronicSubShell>{
			        Molecule()[atomIndex], AtomicSubShellFromContractedSphericalOrbital(orbitalIndex, atomIndex)};
		}

		std::pair<Eigen::Index, ElectronicSubShell> AtomIndexAndSubShellFromContractedSphericalOrbital(
		        const Eigen::Index orbitalIndex) const noexcept
		{
			const auto atomIndex = AtomIndexFromContractedSphericalOrbital(orbitalIndex);
			return std::pair{atomIndex, AtomicSubShellFromContractedSphericalOrbital(orbitalIndex, atomIndex)};
		}


	private:
		MolecularBasisSet(const BasisSet::Gaussian::SharedBasisSetLibrary& library,
		                  const SharedMolecule& molecule,
		                  std::vector<ElementaryBasisPtr> basisAssignments)
		    : m_Library(library), m_Molecule(molecule), m_BasisAssignments(std::move(basisAssignments)),
		      m_SegmentationTableOfPrimitiveSphericalOrbitals(CreateSegmentationTableOfSomeKindOfOrbitals(
		              [](const BasisSet::Gaussian::AngularMomentumBlock& amb)
		              { return amb.PrimitiveSphericalOrbitalCount(); })),
		      m_SegmentationTableOfContractedSphericalOrbitals(CreateSegmentationTableOfSomeKindOfOrbitals(
		              [](const BasisSet::Gaussian::AngularMomentumBlock& amb)
		              { return amb.ContractedSphericalOrbitalCount(); })),
		      m_SubShellSegmentationTableOfEachElementaryBasis(CreateSubShellSegmentationTableForEachElementaryBasis())
		{
			assert(molecule.AtomCount() == m_BasisAssignments.size());
		}

		template <typename OrbitalKindSelector>
		std::vector<Eigen::Index> CreateSegmentationTableOfSomeKindOfOrbitals(OrbitalKindSelector select) const noexcept
		{
			std::vector<Eigen::Index> segTable(m_BasisAssignments.size() + 1, 0);

			Eigen::Index offset = 0;
			std::transform(m_BasisAssignments.cbegin(),
			               m_BasisAssignments.cend(),
			               std::next(segTable.begin()),
			               [&offset, select](const ElementaryBasisPtr basisPtr)
			               {
				               const auto& basis = *basisPtr;
				               offset += std::accumulate(basis.cbegin(),
				                                         basis.cend(),
				                                         Eigen::Index{0},
				                                         [select](const Eigen::Index acc,
				                                                  const BasisSet::Gaussian::AngularMomentumBlock& amb)
				                                         { return acc + select(amb); });
				               return offset;
			               });

			return segTable;
		}

		using SubShellSegmentationTableOfEachElementaryBasis =
		        std::vector<std::pair<ElementaryBasisPtr, std::vector<Eigen::Index>>>;

		SubShellSegmentationTableOfEachElementaryBasis CreateSubShellSegmentationTableForEachElementaryBasis() const
		{
			std::vector<ElementaryBasisPtr> assignments = m_BasisAssignments;
			std::sort(assignments.begin(), assignments.end());
			const auto uniqueElementaryBasisCount =
			        std::distance(assignments.begin(), std::unique(assignments.begin(), assignments.end()));

			SubShellSegmentationTableOfEachElementaryBasis segmentationTable(uniqueElementaryBasisCount);
			std::transform(assignments.cbegin(),
			               assignments.cbegin() + uniqueElementaryBasisCount,
			               segmentationTable.begin(),
			               [](const ElementaryBasisPtr basisPtr)
			               {
				               const auto& basis = *basisPtr;
				               const auto segCount = basis.back().AngularMomentum().Value() + 1;

				               std::vector<Eigen::Index> segTable(segCount + 1, 0);
				               Eigen::Index offset = 0;
				               for (const auto& amb : basis)
				               {
					               // offset += amb.HasOrbital() ? amb.ContractedShellCount() : 0;
					               offset += amb.HasOrbital() ? amb.ContractedSphericalOrbitalCount() : 0;
					               segTable[amb.AngularMomentum().Value() + 1] = offset;
				               }
				               for (std::size_t i = 2; i < segTable.size(); ++i)  // `i` starts from 2, not 1, not 0
				               {
					               segTable[i] = std::max(segTable[i - 1], segTable[i]);
				               }

				               return std::pair{basisPtr, segTable};
			               });
			return segmentationTable;
		}

		template <typename OrbitalKindSelector>
		Eigen::Index CountOfSomeKindOfOrbitalOf(const Atom& atom, OrbitalKindSelector select) const
		{
			const auto& basis = ElementaryBasisOf(atom);
			return std::accumulate(basis.cbegin(),
			                       basis.cend(),
			                       Eigen::Index{0},
			                       [select](const Eigen::Index acc, const BasisSet::Gaussian::AngularMomentumBlock& amb)
			                       { return acc + select(amb); });
		}

		template <typename OrbitalKindSelector>
		Eigen::Index CountOfSomeKindOfOrbital(OrbitalKindSelector select) const noexcept
		{
			return std::accumulate(m_BasisAssignments.cbegin(),
			                       m_BasisAssignments.cend(),
			                       Eigen::Index{0},
			                       [select](const Eigen::Index acc, const ElementaryBasisPtr basisPtr)
			                       {
				                       const auto& basis = *basisPtr;
				                       return acc
				                              + std::accumulate(
				                                      basis.cbegin(),
				                                      basis.cend(),
				                                      Eigen::Index{0},
				                                      [select](const Eigen::Index innerAcc,
				                                               const BasisSet::Gaussian::AngularMomentumBlock& amb)
				                                      { return innerAcc + select(amb); });
			                       });
		}

		/// <c>atomIndex</c> must be resolved by IndexOfAtomOfContractedSphericalOrbital(orbitalIndex).
		/// the behavior is undefined otherwise.
		/// do not make this method public
		ElectronicSubShell AtomicSubShellFromContractedSphericalOrbital(const Eigen::Index orbitalIndex,
		                                                                const Eigen::Index atomIndex) const noexcept
		{
			const auto atomicOffset = m_SegmentationTableOfContractedSphericalOrbitals[atomIndex];
			const auto* basisPtr = m_BasisAssignments[atomIndex];

			const auto subshellSegTable =
			        std::lower_bound(m_SubShellSegmentationTableOfEachElementaryBasis.cbegin(),
			                         m_SubShellSegmentationTableOfEachElementaryBasis.cend(),
			                         basisPtr,
			                         [](const std::pair<ElementaryBasisPtr, std::vector<Eigen::Index>>& lhs,
			                            const ElementaryBasisPtr rhs) { return lhs.first < rhs; })
			                ->second;
			assert(subshellSegTable.size() >= 2);
			const auto segIterator = std::prev(std::upper_bound(
			        std::next(subshellSegTable.cbegin()), subshellSegTable.cend(), orbitalIndex - atomicOffset));
			const AzimuthalQuantumNumber angularMomentum{
			        static_cast<int>(std::distance(subshellSegTable.cbegin(), segIterator))};
			const auto principalQuantumNumber = static_cast<int>(orbitalIndex - *segIterator - atomicOffset)
			                                            / angularMomentum.MagneticQuantumNumberCount()
			                                    + 1 + angularMomentum.Value();

			return ElectronicSubShell{principalQuantumNumber, angularMomentum};
		}

		BasisSet::Gaussian::SharedBasisSetLibrary m_Library;
		SharedMolecule m_Molecule;
		std::vector<ElementaryBasisPtr> m_BasisAssignments;
		std::vector<Eigen::Index> m_SegmentationTableOfPrimitiveSphericalOrbitals;
		std::vector<Eigen::Index> m_SegmentationTableOfContractedSphericalOrbitals;
		SubShellSegmentationTableOfEachElementaryBasis m_SubShellSegmentationTableOfEachElementaryBasis;
	};
}  // namespace SecChem::Gaussian


template <>
class SecChem::Builder<SecChem::Gaussian::MolecularBasisSet>
{
public:
	explicit Builder(const BasisSet::Gaussian::SharedBasisSetLibrary& library) : m_Library(library)
	{
		/* NO CODE */
	}

	Builder& SetOrOverwriteGlobalDefaultBasisSetTo(const std::string& globalDefaultBasisSetName)
	{
		if (!m_Library.Has(globalDefaultBasisSetName))
		{
			throw std::runtime_error("Given basis set library does not contain basis set named \""
			                         + globalDefaultBasisSetName + '"');
		}

		m_WasGlobalDefaultBasisSetSet = true;

		for (const auto& [element, basis] : m_Library[globalDefaultBasisSetName])
		{
			m_DefaultBasisSet.emplace(element, &basis);
		}

		return *this;
	}

	Builder& SetOrOverwriteElementaryDefaultBasisSetTo(const std::string& defaultBasisSetName, const Element element)
	{
		if (!m_WasGlobalDefaultBasisSetSet)
		{
			throw std::runtime_error(
			        "The elementary default basis set can only be set after the global default was set");
		}

		if (!m_Library.Has(defaultBasisSetName))
		{
			throw std::runtime_error("Given basis set library does not contain basis set named \"" + defaultBasisSetName
			                         + '"');
		}

		if (!m_Library[defaultBasisSetName].Has(element))
		{
			throw std::runtime_error("Basis set named \"" + defaultBasisSetName + "\" does not contain "
			                         + element.ToString());
		}

		m_DefaultBasisSet[element] = &m_Library[defaultBasisSetName][element];

		return *this;
	}

	template <typename Parser>
	Gaussian::MolecularBasisSet BuildWith(Parser parser, std::istream& input)
	{
		std::vector<Atom> atoms{};
		std::vector<Gaussian::MolecularBasisSet::ElementaryBasisPtr> assignments{};

		for (std::string line; std::getline(input, line); /* NO CODE */)
		{
			const auto [atom, basisPtr] = parser.ParseLine(
			        line, std::as_const(atoms), std::as_const(m_Library), std::as_const(m_DefaultBasisSet));
			assert(basisPtr != nullptr);

			assignments.emplace_back(basisPtr);
			atoms.push_back(std::move(atom));
		}

		return {m_Library, SharedMolecule{std::move(atoms)}, std::move(assignments)};
	}

	template <typename Parser>
	Gaussian::MolecularBasisSet BuildWith(Parser parser, const std::string& input)
	{
		std::istringstream stream{input};
		return BuildWith(parser, stream);
	}

	template <typename Parser>
	Gaussian::MolecularBasisSet BuildWith(Parser parser, const char* input)
	{
		return BuildWith(parser, std::string{input});
	}


private:
	void SetDefaultBasisSet(const std::string& defaultBasisSetName)
	{
		for (const auto& [element, basis] : m_Library[defaultBasisSetName])
		{
			m_DefaultBasisSet.emplace(element, &basis);
		}
	}

	BasisSet::Gaussian::SharedBasisSetLibrary m_Library;
	std::unordered_map<Element, const std::vector<BasisSet::Gaussian::AngularMomentumBlock>*> m_DefaultBasisSet;
	bool m_WasGlobalDefaultBasisSetSet = false;
};
