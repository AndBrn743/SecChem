// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Andy Brown

#pragma once
#if !defined(SECCHEM_GAUSSIAN_BASIS_SET_INTERNAL)
#error Do not include internal header files directly
#endif

#include <range/v3/algorithm.hpp>
#include <range/v3/numeric.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/transform.hpp>
#include <sstream>
#include <string>
#include <utility>

#include <SecChem/Molecule.hpp>
#include <SecChem/Polyfill/flat_map.hpp>


namespace SecChem::BasisSet::Gaussian
{
	class MolecularBasisSet
	{
		friend Builder<MolecularBasisSet>;

	public:
		// GCC require it
		// ReSharper disable once CppRedundantQualifier
		const Gaussian::SharedBasisSetLibrary& SharedBasisSetLibrary() const noexcept
		{
			return m_Library;
		}

		const SharedMolecule& Molecule() const noexcept
		{
			return m_Molecule;
		}

		auto ElementaryBasisSets() const noexcept
		{
			return m_BasisAssignments
			       | ranges::views::transform([](const ElementaryBasisSet* ptr) -> const ElementaryBasisSet&
			                                  { return *ptr; });
		}

		const ElementaryBasisSet& ElementaryBasisAt(const std::size_t index) const
		{
			if (index >= m_BasisAssignments.size())
			{
				throw std::out_of_range("Basis assignment index out of range");
			}

			return *m_BasisAssignments[index];
		}

		const ElementaryBasisSet& ElementaryBasisOf(const Atom& atom) const
		{
			return *m_BasisAssignments[m_Molecule.IndexOf(atom)];
		}

		std::size_t UniqueElementaryBasisSetCount() const noexcept
		{
			return m_ComputedElementaryBasisInfoTable.size();
		}

		auto UniqueElementaryBasisSet() const noexcept
		{
			// due to the limitation of range-v3 0.12.0 and std::range of C++20 (C++23's std::range is fine) we can only
			// use non-const flat_map. this should be fine, we just have to be *VERY* careful
			return const_cast<polyfill::flat_map<const ElementaryBasisSet*, ComputedElementaryBasisInfo>&>(
			               m_ComputedElementaryBasisInfoTable)
			       | ranges::views::transform([](const auto& kv) -> const ElementaryBasisSet&
			                                  { return *std::get<0>(kv); });
		}

		Eigen::Index ContractedSubShellCount() const noexcept
		{
			return std::accumulate(m_ComputedElementaryBasisInfoTable.cbegin(),
			                       m_ComputedElementaryBasisInfoTable.cend(),
			                       Eigen::Index{0},
			                       [](const Eigen::Index acc, const auto& kv)
			                       {
				                       const ComputedElementaryBasisInfo& info = std::get<1>(kv);
				                       return acc + info.ReferenceCount * info.ContractedSubShellCount;
			                       });
		}

		Eigen::Index PrimitiveSubShellCount() const noexcept
		{
			return std::accumulate(m_ComputedElementaryBasisInfoTable.cbegin(),
			                       m_ComputedElementaryBasisInfoTable.cend(),
			                       Eigen::Index{0},
			                       [](const Eigen::Index acc, const auto& kv)
			                       {
				                       const ComputedElementaryBasisInfo& info = std::get<1>(kv);
				                       return acc + info.ReferenceCount * info.PrimitiveSubShellCount;
			                       });
		}

		// this method DOES cache
		Eigen::Index PrimitiveSphericalOrbitalCountOf(const Atom& atom) const
		{
			const auto index = Molecule().IndexOf(atom);
			return m_PrimitiveSphericalOrbitalSegmentationTable[index + 1]
			       - m_PrimitiveSphericalOrbitalSegmentationTable[index];
		}

		// this method does NOT cache
		Eigen::Index PrimitiveCartesianOrbitalCountOf(const Atom& atom) const
		{
			return CountOfSomeKindOfOrbitalOf(atom, &AngularMomentumBlock::PrimitiveCartesianOrbitalCount);
		}

		// this method DOES cache
		Eigen::Index ContractedSphericalOrbitalCountOf(const Atom& atom) const
		{
			const auto index = Molecule().IndexOf(atom);
			return m_ContractedSphericalOrbitalSegmentationTable[index + 1]
			       - m_ContractedSphericalOrbitalSegmentationTable[index];
		}

		// this method does NOT cache
		Eigen::Index ContractedCartesianOrbitalCountOf(const Atom& atom) const
		{
			return CountOfSomeKindOfOrbitalOf(atom, &AngularMomentumBlock::ContractedCartesianOrbitalCount);
		}

		// this method DOES cache
		Eigen::Index PrimitiveSphericalOrbitalCount() const noexcept
		{
			return m_PrimitiveSphericalOrbitalSegmentationTable.back();
		}

		// this method does NOT cache
		Eigen::Index PrimitiveCartesianOrbitalCount() const noexcept
		{
			return CountOfSomeKindOfOrbital(&AngularMomentumBlock::PrimitiveCartesianOrbitalCount);
		}

		// this method DOES cache
		Eigen::Index ContractedSphericalOrbitalCount() const noexcept
		{
			return m_ContractedSphericalOrbitalSegmentationTable.back();
		}

		// this method does NOT cache
		Eigen::Index ContractedCartesianOrbitalCount() const noexcept
		{
			return CountOfSomeKindOfOrbital(&AngularMomentumBlock::ContractedCartesianOrbitalCount);
		}

		// this method does NOT cache
		Eigen::Index PrimitiveSubShellCountOf(const Atom& atom) const
		{
			return ranges::accumulate(ElementaryBasisOf(atom).AngularMomentumBlocks,
			                          Eigen::Index{0},
			                          std::plus<>{},
			                          &AngularMomentumBlock::PrimitiveShellCount);
		}

		// this method does NOT cache
		Eigen::Index ContractedSubShellCountOf(const Atom& atom) const
		{
			return ranges::accumulate(ElementaryBasisOf(atom).AngularMomentumBlocks,
			                          Eigen::Index{0},
			                          std::plus<>{},
			                          &AngularMomentumBlock::ContractedShellCount);
		}

		auto PrimitiveSubShellsOfAtomAt(const std::size_t index) const
		{
			return ElementaryBasisAt(index).AngularMomentumBlocks
			       | ranges::views::transform([](const AngularMomentumBlock& amb) { return amb.PrimitiveShells(); })
			       | ranges::views::join;
		}

		auto ContractedSubShellsOfAtomAt(const std::size_t index) const
		{
			return ElementaryBasisAt(index).AngularMomentumBlocks
			       | ranges::views::transform([](const AngularMomentumBlock& amb) { return amb.ContractedShells(); })
			       | ranges::views::join;
		}

		auto EcpOffsettedPrimitiveSubShellsOfAtomAt(const std::size_t index) const
		{
			const auto* basisPtr = m_BasisAssignments[index];
			const auto& ambs = basisPtr->AngularMomentumBlocks;

			const auto* ambsPtr = ambs.data();
			const auto* principalQuantumNumberOffsetsPtr =
			        m_ComputedElementaryBasisInfoTable.at(basisPtr).PrincipalQuantumNumberOffsetTable.data();

			return ambs
			       | ranges::views::transform(
			               [principalQuantumNumberOffsetsPtr, ambsPtr](const AngularMomentumBlock& amb)
			               { return amb.PrimitiveShells(principalQuantumNumberOffsetsPtr[&amb - ambsPtr]); })
			       | ranges::views::join;
		}

		auto EcpOffsettedContractedSubShellsOfAtomAt(const std::size_t index) const
		{
			const auto* basisPtr = m_BasisAssignments[index];
			const auto& ambs = basisPtr->AngularMomentumBlocks;

			const auto* ambsPtr = ambs.data();
			const auto* principalQuantumNumberOffsetsPtr =
			        m_ComputedElementaryBasisInfoTable.at(basisPtr).PrincipalQuantumNumberOffsetTable.data();

			return ambs
			       | ranges::views::transform(
			               [principalQuantumNumberOffsetsPtr, ambsPtr](const AngularMomentumBlock& amb)
			               { return amb.ContractedShells(principalQuantumNumberOffsetsPtr[&amb - ambsPtr]); })
			       | ranges::views::join;
		}

		auto PrimitiveSubShellsOf(const Atom& atom) const
		{
			return PrimitiveSubShellsOfAtomAt(m_Molecule.IndexOf(atom));
		}

		auto ContractedSubShellsOf(const Atom& atom) const
		{
			return ContractedSubShellsOfAtomAt(m_Molecule.IndexOf(atom));
		}

		auto EcpOffsettedPrimitiveSubShellsOf(const Atom& atom) const
		{
			return EcpOffsettedPrimitiveSubShellsOfAtomAt(m_Molecule.IndexOf(atom));
		}

		auto EcpOffsettedContractedSubShellsOf(const Atom& atom) const
		{
			return EcpOffsettedContractedSubShellsOfAtomAt(m_Molecule.IndexOf(atom));
		}

		auto ContractedSphericalOrbitalsFrom(const Atom& atom) const
		{
			const auto atomIndex = m_Molecule.IndexOf(atom);
			return ranges::views::iota(m_ContractedSphericalOrbitalSegmentationTable[atomIndex],
			                           m_ContractedSphericalOrbitalSegmentationTable[atomIndex + 1]);
		}

		auto ContractedSphericalOrbitalsFrom(const Atom& atom, const ElectronicSubShell shell) const
		{
			const auto atomIndex = m_Molecule.IndexOf(atom);
			const auto atomicOffset = m_ContractedSphericalOrbitalSegmentationTable[atomIndex];
			const ElementaryBasisSet* basisPtr = m_BasisAssignments[atomIndex];
			const auto azimuthalShellOffset =
			        m_ComputedElementaryBasisInfoTable.at(basisPtr)
			                .ContractedSphericalSubShellSegmentationTable[shell.AzimuthalQuantumNumber().Value()];
			const auto sphericalMagneticQuantumNumberCount = shell.MagneticQuantumNumberCount();
			const auto offset = atomicOffset + azimuthalShellOffset
			                    + sphericalMagneticQuantumNumberCount
			                              * (shell.PrincipalQuantumNumber()
			                                 - shell.AzimuthalQuantumNumber().MinPrincipalQuantumNumber());

			return ranges::views::iota(offset, offset + sphericalMagneticQuantumNumberCount);
		}

		template <typename VectorLike>
		auto ContractedSphericalOrbitalSegmentOf(VectorLike&& vector,
		                                         const Atom& atom,
		                                         const ElectronicSubShell shell) const
		{
			const auto atomIndex = m_Molecule.IndexOf(atom);
			const auto atomicOffset = m_ContractedSphericalOrbitalSegmentationTable[atomIndex];
			const ElementaryBasisSet* basisPtr = m_BasisAssignments[atomIndex];
			const auto azimuthalShellOffset =
			        m_ComputedElementaryBasisInfoTable.at(basisPtr)
			                .ContractedSphericalSubShellSegmentationTable[shell.AzimuthalQuantumNumber().Value()];
			const auto sphericalMagneticQuantumNumberCount = shell.MagneticQuantumNumberCount();
			const auto offset = atomicOffset + azimuthalShellOffset
			                    + sphericalMagneticQuantumNumberCount
			                              * (shell.PrincipalQuantumNumber()
			                                 - shell.AzimuthalQuantumNumber().MinPrincipalQuantumNumber());

			return vector.segment(offset, sphericalMagneticQuantumNumberCount);
		}

		template <typename VectorLike>
		auto ContractedSphericalOrbitalSegmentOf(VectorLike&& vector, const Atom& atom) const
		{
			const auto atomIndex = m_Molecule.IndexOf(atom);
			return vector.segment(m_ContractedSphericalOrbitalSegmentationTable[atomIndex],
			                      m_ContractedSphericalOrbitalSegmentationTable[atomIndex + 1]
			                              - m_ContractedSphericalOrbitalSegmentationTable[atomIndex]);
		}

		Eigen::Index AtomIndexFromContractedSphericalOrbital(const Eigen::Index orbitalIndex) const noexcept
		{
			assert(orbitalIndex >= 0 && orbitalIndex < m_ContractedSphericalOrbitalSegmentationTable.back());

			return std::distance(
			        m_ContractedSphericalOrbitalSegmentationTable.cbegin(),
			        std::prev(ranges::upper_bound(m_ContractedSphericalOrbitalSegmentationTable, orbitalIndex)));
		}

		const Atom& AtomFromContractedSphericalOrbital(const Eigen::Index orbitalIndex) const noexcept
		{
			return Molecule()[AtomIndexFromContractedSphericalOrbital(orbitalIndex)];
		}

		Eigen::Index AtomIndexFromPrimitiveSphericalOrbital(const Eigen::Index orbitalIndex) const noexcept
		{
			assert(orbitalIndex >= 0 && orbitalIndex < m_PrimitiveSphericalOrbitalSegmentationTable.back());

			return std::distance(
			        m_PrimitiveSphericalOrbitalSegmentationTable.cbegin(),
			        std::prev(ranges::upper_bound(m_PrimitiveSphericalOrbitalSegmentationTable, orbitalIndex)));
		}

		const Atom& AtomFromPrimitiveSphericalOrbital(const Eigen::Index orbitalIndex) const noexcept
		{
			return Molecule()[AtomIndexFromPrimitiveSphericalOrbital(orbitalIndex)];
		}

		ElectronicSubShell AtomicSubShellFromContractedSphericalOrbital(const Eigen::Index orbitalIndex) const noexcept
		{
			return AtomicSubShellFromContractedSphericalOrbital(orbitalIndex,
			                                                    AtomIndexFromContractedSphericalOrbital(orbitalIndex));
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
		MolecularBasisSet(Gaussian::SharedBasisSetLibrary library,
		                  const SharedMolecule& molecule,
		                  std::vector<const ElementaryBasisSet*> basisAssignments)
		    : m_Library(std::move(library)), m_Molecule(molecule), m_BasisAssignments(std::move(basisAssignments)),
		      m_PrimitiveSphericalOrbitalSegmentationTable(CreateSegmentationTableOfSomeKindOfOrbitals(
		              [](const AngularMomentumBlock& amb) { return amb.PrimitiveSphericalOrbitalCount(); })),
		      m_ContractedSphericalOrbitalSegmentationTable(CreateSegmentationTableOfSomeKindOfOrbitals(
		              [](const AngularMomentumBlock& amb) { return amb.ContractedSphericalOrbitalCount(); })),
		      m_ComputedElementaryBasisInfoTable(ComputedElementaryBasisInfo::CreateStatisticsFor(m_BasisAssignments))
		{
			assert(molecule.AtomCount() == m_BasisAssignments.size());
		}

		template <typename OrbitalCounter>
		std::vector<Eigen::Index> CreateSegmentationTableOfSomeKindOfOrbitals(
		        OrbitalCounter orbitalCountOf) const noexcept
		{
			std::vector<Eigen::Index> segTable(m_BasisAssignments.size() + 1, 0);

			Eigen::Index offset = 0;
			ranges::transform(m_BasisAssignments,
			                  std::next(segTable.begin()),
			                  [&offset, orbitalCountOf](const ElementaryBasisSet* basisPtr)
			                  {
				                  offset += ranges::accumulate(basisPtr->AngularMomentumBlocks,
				                                               Eigen::Index{0},
				                                               std::plus<>{},
				                                               orbitalCountOf);
				                  return offset;
			                  });

			return segTable;
		}

		struct ComputedElementaryBasisInfo
		{
			static polyfill::flat_map<const ElementaryBasisSet*, ComputedElementaryBasisInfo> CreateStatisticsFor(
			        const std::vector<const ElementaryBasisSet*>& basisAssignments)
			{
				polyfill::flat_map<const ElementaryBasisSet*, ComputedElementaryBasisInfo> statistics;

				for (const ElementaryBasisSet* basisPtr : basisAssignments)
				{
					auto& [referenceCount,
					       primitiveSphericalSubShellSegTable,
					       contractedSphericalSubShellSegTable,
					       primitiveSubShellCount,
					       contractedSubShellCount,
					       principalQuantumNumberOffsetTable] = statistics[basisPtr];

					referenceCount++;

					if (referenceCount == 1)
					{
						primitiveSphericalSubShellSegTable = CreateSubShellSegmentationTableFor(
						        *basisPtr,
						        [](const AngularMomentumBlock& amb) { return amb.PrimitiveSphericalOrbitalCount(); });

						contractedSphericalSubShellSegTable = CreateSubShellSegmentationTableFor(
						        *basisPtr,
						        [](const AngularMomentumBlock& amb) { return amb.ContractedSphericalOrbitalCount(); });

						primitiveSubShellCount = ranges::accumulate(basisPtr->AngularMomentumBlocks,
						                                            Eigen::Index{0},
						                                            std::plus<>{},
						                                            &AngularMomentumBlock::PrimitiveShellCount);

						contractedSubShellCount = ranges::accumulate(basisPtr->AngularMomentumBlocks,
						                                             Eigen::Index{0},
						                                             std::plus<>{},
						                                             &AngularMomentumBlock::ContractedShellCount);

						principalQuantumNumberOffsetTable = basisPtr->CreatePrincipalQuantumNumberOffsetTable();
					}
				}

				return statistics;
			}

			Eigen::Index ReferenceCount;
			std::vector<Eigen::Index> PrimitiveSphericalSubShellSegmentationTable;
			std::vector<Eigen::Index> ContractedSphericalSubShellSegmentationTable;
			Eigen::Index PrimitiveSubShellCount;
			Eigen::Index ContractedSubShellCount;
			std::vector<int> PrincipalQuantumNumberOffsetTable;
		};

		using SubShellSegmentationTableOfEachElementaryBasis =
		        std::vector<std::pair<const ElementaryBasisSet*, std::vector<Eigen::Index>>>;

		template <typename OrbitalCounter>
		static std::vector<Eigen::Index> CreateSubShellSegmentationTableFor(const ElementaryBasisSet& basis,
		                                                                    OrbitalCounter orbitalCountOf)
		{
			const auto segCount = basis.AngularMomentumBlocks.back().AngularMomentum().Value() + 1;

			std::vector<Eigen::Index> segTable(segCount + 1, 0);
			Eigen::Index offset = 0;

			for (const auto& amb : basis.AngularMomentumBlocks)
			{
				offset += orbitalCountOf(amb);
				segTable[amb.AngularMomentum().Value() + 1] = offset;
			}

			for (std::size_t i = 2; i < segTable.size(); ++i)  // `i` starts from 2, not 1, not 0
			{
				segTable[i] = std::max(segTable[i - 1], segTable[i]);
			}

			return segTable;
		}

		template <typename OrbitalCounter>
		Eigen::Index CountOfSomeKindOfOrbitalOf(const Atom& atom, OrbitalCounter orbitalCountOf) const
		{
			return ranges::accumulate(
			        ElementaryBasisOf(atom).AngularMomentumBlocks, Eigen::Index{0}, std::plus<>{}, orbitalCountOf);
		}

		template <typename OrbitalCounter>
		Eigen::Index CountOfSomeKindOfOrbital(OrbitalCounter orbitalCountOf) const noexcept
		{
			return ranges::accumulate(m_BasisAssignments,
			                          Eigen::Index{0},
			                          [orbitalCountOf](const Eigen::Index acc, const ElementaryBasisSet* basisPtr)
			                          {
				                          return acc
				                                 + ranges::accumulate(basisPtr->AngularMomentumBlocks,
				                                                      Eigen::Index{0},
				                                                      std::plus<>{},
				                                                      orbitalCountOf);
			                          });
		}

		/// <c>atomIndex</c> must be resolved by IndexOfAtomOfContractedSphericalOrbital(orbitalIndex).
		/// the behavior is undefined otherwise.
		/// do not make this method public
		ElectronicSubShell AtomicSubShellFromContractedSphericalOrbital(const Eigen::Index orbitalIndex,
		                                                                const Eigen::Index atomIndex) const noexcept
		{
			const auto atomicOffset = m_ContractedSphericalOrbitalSegmentationTable[atomIndex];
			const ElementaryBasisSet* basisPtr = m_BasisAssignments[atomIndex];

			const auto& subshellSegTable =
			        m_ComputedElementaryBasisInfoTable.at(basisPtr).ContractedSphericalSubShellSegmentationTable;
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

		Gaussian::SharedBasisSetLibrary m_Library;
		SharedMolecule m_Molecule;
		std::vector<const ElementaryBasisSet*> m_BasisAssignments{};
		std::vector<Eigen::Index> m_PrimitiveSphericalOrbitalSegmentationTable{};
		std::vector<Eigen::Index> m_ContractedSphericalOrbitalSegmentationTable{};
		polyfill::flat_map<const ElementaryBasisSet*, ComputedElementaryBasisInfo> m_ComputedElementaryBasisInfoTable{};
	};
}  // namespace SecChem::BasisSet::Gaussian


template <>
class SecChem::Builder<SecChem::BasisSet::Gaussian::MolecularBasisSet>
{
public:
	explicit Builder(BasisSet::Gaussian::SharedBasisSetLibrary library) : m_Library(std::move(library))
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
	BasisSet::Gaussian::MolecularBasisSet BuildWith(Parser parser, std::istream& input)
	{
		std::vector<Atom> atoms{};
		std::vector<const BasisSet::Gaussian::ElementaryBasisSet*> assignments{};

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
	BasisSet::Gaussian::MolecularBasisSet BuildWith(Parser parser, const std::string& input)
	{
		std::istringstream stream{input};
		return BuildWith(parser, stream);
	}

	template <typename Parser>
	BasisSet::Gaussian::MolecularBasisSet BuildWith(Parser parser, const char* input)
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
	std::unordered_map<Element, const BasisSet::Gaussian::ElementaryBasisSet*> m_DefaultBasisSet;
	bool m_WasGlobalDefaultBasisSetSet = false;
};
