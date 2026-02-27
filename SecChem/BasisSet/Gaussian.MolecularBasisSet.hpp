// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Andy Brown

#pragma once

#include <functional>
#include <range/v3/algorithm.hpp>
#include <range/v3/numeric.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/transform.hpp>
#include <sstream>
#include <string>
#include <utility>

#include <SecChem/Molecule.hpp>
#include <SecChem/Polyfill/flat_map.hpp>

#include "Gaussian.BasisSet.hpp"

namespace SecChem::BasisSet::Gaussian
{
	namespace Detail::MolecularBasisSet
	{
		enum class Contraction
		{
			Unspecified,
			Primitive,
			Contracted
		};

		enum class Representation
		{
			Unspecified,
			Cartesian,
			Spherical
		};

		template <Contraction C, Representation R>
		class BasisView;
	}  // namespace Detail::MolecularBasisSet

	class MolecularBasisSet
	{
		friend Builder<MolecularBasisSet>;
		friend Detail::MolecularBasisSet::BasisView<Detail::MolecularBasisSet::Contraction::Primitive,
		                                            Detail::MolecularBasisSet::Representation::Cartesian>;
		friend Detail::MolecularBasisSet::BasisView<Detail::MolecularBasisSet::Contraction::Primitive,
		                                            Detail::MolecularBasisSet::Representation::Spherical>;
		friend Detail::MolecularBasisSet::BasisView<Detail::MolecularBasisSet::Contraction::Contracted,
		                                            Detail::MolecularBasisSet::Representation::Cartesian>;
		friend Detail::MolecularBasisSet::BasisView<Detail::MolecularBasisSet::Contraction::Contracted,
		                                            Detail::MolecularBasisSet::Representation::Spherical>;

		friend Detail::MolecularBasisSet::BasisView<Detail::MolecularBasisSet::Contraction::Unspecified,
		                                            Detail::MolecularBasisSet::Representation::Cartesian>;
		friend Detail::MolecularBasisSet::BasisView<Detail::MolecularBasisSet::Contraction::Unspecified,
		                                            Detail::MolecularBasisSet::Representation::Spherical>;
		friend Detail::MolecularBasisSet::BasisView<Detail::MolecularBasisSet::Contraction::Primitive,
		                                            Detail::MolecularBasisSet::Representation::Unspecified>;
		friend Detail::MolecularBasisSet::BasisView<Detail::MolecularBasisSet::Contraction::Contracted,
		                                            Detail::MolecularBasisSet::Representation::Unspecified>;

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

		Detail::MolecularBasisSet::BasisView<Detail::MolecularBasisSet::Contraction::Primitive,
		                                     Detail::MolecularBasisSet::Representation::Unspecified>
		Primitive() const;
		Detail::MolecularBasisSet::BasisView<Detail::MolecularBasisSet::Contraction::Contracted,
		                                     Detail::MolecularBasisSet::Representation::Unspecified>
		Contracted() const;
		Detail::MolecularBasisSet::BasisView<Detail::MolecularBasisSet::Contraction::Unspecified,
		                                     Detail::MolecularBasisSet::Representation::Cartesian>
		Cartesian() const;
		Detail::MolecularBasisSet::BasisView<Detail::MolecularBasisSet::Contraction::Unspecified,
		                                     Detail::MolecularBasisSet::Representation::Spherical>
		Spherical() const;


	private:
		MolecularBasisSet(Gaussian::SharedBasisSetLibrary library,
		                  const SharedMolecule& molecule,
		                  std::vector<const ElementaryBasisSet*> basisAssignments)
		    : m_Library(std::move(library)), m_Molecule(molecule), m_BasisAssignments(std::move(basisAssignments)),
		      m_PrimitiveSphericalOrbitalSegmentationTable(
		              CreateSegmentationTableOfSomeKindOfOrbitals(&AzimuthalShell::PrimitiveSphericalOrbitalCount)),
		      m_ContractedSphericalOrbitalSegmentationTable(
		              CreateSegmentationTableOfSomeKindOfOrbitals(&AzimuthalShell::ContractedSphericalOrbitalCount)),
		      m_PrimitiveCartesianOrbitalSegmentationTable(
		              CreateSegmentationTableOfSomeKindOfOrbitals(&AzimuthalShell::PrimitiveCartesianOrbitalCount)),
		      m_ContractedCartesianOrbitalSegmentationTable(
		              CreateSegmentationTableOfSomeKindOfOrbitals(&AzimuthalShell::ContractedCartesianOrbitalCount)),
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
				                  offset += ranges::accumulate(
				                          basisPtr->AzimuthalShells, Eigen::Index{0}, std::plus<>{}, orbitalCountOf);
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
					       primitiveSphericalSubshellSegTable,
					       contractedSphericalSubshellSegTable,
					       primitiveCartesianSubshellSegTable,
					       contractedCartesianSubshellSegTable,
					       primitiveSubshellCount,
					       contractedSubshellCount,
					       principalQuantumNumberOffsetTable] = statistics[basisPtr];

					referenceCount++;

					if (referenceCount == 1)
					{
						primitiveSphericalSubshellSegTable = CreateSubshellSegmentationTableFor(
						        *basisPtr, &AzimuthalShell::PrimitiveSphericalOrbitalCount);

						contractedSphericalSubshellSegTable = CreateSubshellSegmentationTableFor(
						        *basisPtr, &AzimuthalShell::ContractedSphericalOrbitalCount);

						primitiveCartesianSubshellSegTable = CreateSubshellSegmentationTableFor(
						        *basisPtr, &AzimuthalShell::PrimitiveCartesianOrbitalCount);

						contractedCartesianSubshellSegTable = CreateSubshellSegmentationTableFor(
						        *basisPtr, &AzimuthalShell::ContractedCartesianOrbitalCount);

						primitiveSubshellCount = ranges::accumulate(basisPtr->AzimuthalShells,
						                                            Eigen::Index{0},
						                                            std::plus<>{},
						                                            &AzimuthalShell::PrimitiveShellCount);

						contractedSubshellCount = ranges::accumulate(basisPtr->AzimuthalShells,
						                                             Eigen::Index{0},
						                                             std::plus<>{},
						                                             &AzimuthalShell::ContractedShellCount);

						principalQuantumNumberOffsetTable = basisPtr->CreatePrincipalQuantumNumberOffsetTable();
					}
				}

				return statistics;
			}

			template <Detail::MolecularBasisSet::Contraction C, Detail::MolecularBasisSet::Representation R>
			const auto& SubshellSegmentationTable() const noexcept
			{
				static_assert((C == Detail::MolecularBasisSet::Contraction::Contracted
				               || C == Detail::MolecularBasisSet::Contraction::Primitive)
				              && (R == Detail::MolecularBasisSet::Representation::Cartesian
				                  || R == Detail::MolecularBasisSet::Representation::Spherical));

				if constexpr (C == Detail::MolecularBasisSet::Contraction::Contracted)
				{
					if constexpr (R == Detail::MolecularBasisSet::Representation::Cartesian)
					{
						return ContractedCartesianSubshellSegmentationTable;
					}
					else
					{
						return ContractedSphericalSubshellSegmentationTable;
					}
				}
				else
				{
					if constexpr (R == Detail::MolecularBasisSet::Representation::Cartesian)
					{
						return PrimitiveCartesianSubshellSegmentationTable;
					}
					else
					{
						return PrimitiveSphericalSubshellSegmentationTable;
					}
				}
			}


			template <Detail::MolecularBasisSet::Contraction C>
			const auto& SubshellCount() const noexcept
			{
				static_assert(C == Detail::MolecularBasisSet::Contraction::Contracted
				              || C == Detail::MolecularBasisSet::Contraction::Primitive);

				return C == Detail::MolecularBasisSet::Contraction::Contracted ? ContractedSubshellCount
				                                                               : PrimitiveSubshellCount;
			}

			Eigen::Index ReferenceCount;
			std::vector<Eigen::Index> PrimitiveSphericalSubshellSegmentationTable;
			std::vector<Eigen::Index> ContractedSphericalSubshellSegmentationTable;
			std::vector<Eigen::Index> PrimitiveCartesianSubshellSegmentationTable;
			std::vector<Eigen::Index> ContractedCartesianSubshellSegmentationTable;
			Eigen::Index PrimitiveSubshellCount;
			Eigen::Index ContractedSubshellCount;
			std::vector<int> PrincipalQuantumNumberOffsetTable;
		};

		using SubshellSegmentationTableOfEachElementaryBasis =
		        std::vector<std::pair<const ElementaryBasisSet*, std::vector<Eigen::Index>>>;

		template <typename OrbitalCounter>
		static std::vector<Eigen::Index> CreateSubshellSegmentationTableFor(const ElementaryBasisSet& basis,
		                                                                    OrbitalCounter orbitalCountOf)
		{
			const auto segCount = basis.AzimuthalShells.back().AngularMomentum().Value() + 1;

			std::vector<Eigen::Index> segTable(segCount + 1, 0);
			Eigen::Index offset = 0;

			for (const auto& amb : basis.AzimuthalShells)
			{
				offset += std::invoke(orbitalCountOf, amb);
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
			        ElementaryBasisOf(atom).AzimuthalShells, Eigen::Index{0}, std::plus<>{}, orbitalCountOf);
		}

		template <typename OrbitalCounter>
		Eigen::Index CountOfSomeKindOfOrbital(OrbitalCounter orbitalCountOf) const noexcept
		{
			return ranges::accumulate(m_BasisAssignments,
			                          Eigen::Index{0},
			                          [orbitalCountOf](const Eigen::Index acc, const ElementaryBasisSet* basisPtr)
			                          {
				                          return acc
				                                 + ranges::accumulate(basisPtr->AzimuthalShells,
				                                                      Eigen::Index{0},
				                                                      std::plus<>{},
				                                                      orbitalCountOf);
			                          });
		}

		template <Detail::MolecularBasisSet::Contraction C, Detail::MolecularBasisSet::Representation R>
		const std::vector<Eigen::Index>& OrbitalSegmentationTable() const noexcept
		{
			static_assert((C == Detail::MolecularBasisSet::Contraction::Contracted
			               || C == Detail::MolecularBasisSet::Contraction::Primitive)
			              && (R == Detail::MolecularBasisSet::Representation::Spherical
			                  || R == Detail::MolecularBasisSet::Representation::Cartesian));

			if constexpr (C == Detail::MolecularBasisSet::Contraction::Contracted)
			{
				if constexpr (R == Detail::MolecularBasisSet::Representation::Spherical)
				{
					return m_ContractedSphericalOrbitalSegmentationTable;
				}
				else
				{
					return m_ContractedCartesianOrbitalSegmentationTable;
				}
			}
			else
			{
				if constexpr (R == Detail::MolecularBasisSet::Representation::Spherical)
				{
					return m_PrimitiveSphericalOrbitalSegmentationTable;
				}
				else
				{
					return m_PrimitiveCartesianOrbitalSegmentationTable;
				}
			}
		}

		Gaussian::SharedBasisSetLibrary m_Library;
		SharedMolecule m_Molecule;
		std::vector<const ElementaryBasisSet*> m_BasisAssignments{};
		std::vector<Eigen::Index> m_PrimitiveSphericalOrbitalSegmentationTable{};
		std::vector<Eigen::Index> m_ContractedSphericalOrbitalSegmentationTable{};
		std::vector<Eigen::Index> m_PrimitiveCartesianOrbitalSegmentationTable{};
		std::vector<Eigen::Index> m_ContractedCartesianOrbitalSegmentationTable{};
		polyfill::flat_map<const ElementaryBasisSet*, ComputedElementaryBasisInfo> m_ComputedElementaryBasisInfoTable{};
	};


	namespace Detail::MolecularBasisSet
	{
		template <Contraction C, Representation R>
		class BasisView
		{
			static_assert((C == Contraction::Contracted || C == Contraction::Primitive)
			              && (R == Representation::Spherical || R == Representation::Cartesian));

			struct OrbitalSegmentDescription
			{
				Eigen::Index First;
				Eigen::Index Count;
			};

		public:
			constexpr explicit BasisView(const Gaussian::MolecularBasisSet& basis) noexcept : cr_BasisSet(basis)
			{
				/* NO CODE */
			}

			// ---- final operation ----

			std::size_t OrbitalCount() const
			{
				return cr_BasisSet.OrbitalSegmentationTable<C, R>().back();
			}

			std::size_t OrbitalCountOf(const Atom& atom) const
			{
				const auto index = cr_BasisSet.Molecule().IndexOf(atom);
				return cr_BasisSet.OrbitalSegmentationTable<C, R>()[index + 1]
				       - cr_BasisSet.OrbitalSegmentationTable<C, R>()[index];
			}

			Eigen::Index SubshellCount() const noexcept
			{
				return BasisView<C, Representation::Unspecified>{cr_BasisSet}.SubshellCount();
			}

			auto SubshellsOf(const Atom& atom) const
			{
				return BasisView<C, Representation::Unspecified>{cr_BasisSet}.SubshellsOf(atom);
			}

			auto EcpOffsettedSubshellsOf(const Atom& atom) const
			{
				return BasisView<C, Representation::Unspecified>{cr_BasisSet}.EcpOffsettedSubshellsOf(atom);
			}

			auto OrbitalsFrom(const Atom& atom) const
			{
				const auto atomIndex = cr_BasisSet.m_Molecule.IndexOf(atom);
				return ranges::views::iota(cr_BasisSet.OrbitalSegmentationTable<C, R>()[atomIndex],
				                           cr_BasisSet.OrbitalSegmentationTable<C, R>()[atomIndex + 1]);
			}

			auto OrbitalsFrom(const Atom& atom, const ElectronicSubshell shell) const
			{
				const auto atomIndex = cr_BasisSet.m_Molecule.IndexOf(atom);
				const auto atomicOffset = cr_BasisSet.OrbitalSegmentationTable<C, R>()[atomIndex];
				const ElementaryBasisSet* basisPtr = cr_BasisSet.m_BasisAssignments[atomIndex];
				const auto azimuthalShellOffset =
				        cr_BasisSet.m_ComputedElementaryBasisInfoTable.at(basisPtr)
				                .template SubshellSegmentationTable<C, R>()[shell.AzimuthalQuantumNumber().Value()];
				const auto magneticQuantumNumberCount = R == Representation::Spherical
				                                                ? shell.MagneticQuantumNumberCount()
				                                                : shell.CartesianMagneticQuantumNumberCount();
				const auto offset = atomicOffset + azimuthalShellOffset
				                    + magneticQuantumNumberCount
				                              * (shell.PrincipalQuantumNumber()
				                                 - shell.AzimuthalQuantumNumber().MinPrincipalQuantumNumber());

				return ranges::views::iota(offset, offset + magneticQuantumNumberCount);
			}


			Eigen::Index AtomIndexFromOrbital(const Eigen::Index orbitalIndex) const noexcept
			{
				assert(orbitalIndex >= 0 && orbitalIndex < (cr_BasisSet.OrbitalSegmentationTable<C, R>().back()));

				return std::distance(
				        cr_BasisSet.OrbitalSegmentationTable<C, R>().cbegin(),
				        std::prev(ranges::upper_bound(cr_BasisSet.OrbitalSegmentationTable<C, R>(), orbitalIndex)));
			}

			const Atom& AtomFromOrbital(const Eigen::Index orbitalIndex) const noexcept
			{
				return cr_BasisSet.m_Molecule[AtomIndexFromOrbital(orbitalIndex)];
			}

			ElectronicSubshell AtomicSubshellFromOrbital(const Eigen::Index orbitalIndex) const noexcept
			{
				return AtomicSubshellFromOrbital(orbitalIndex, AtomIndexFromOrbital(orbitalIndex));
			}

			std::pair<const Atom&, ElectronicSubshell> AtomAndSubshellFromOrbital(
			        const Eigen::Index orbitalIndex) const noexcept
			{
				const auto atomIndex = AtomIndexFromOrbital(orbitalIndex);
				return std::pair<const Atom&, ElectronicSubshell>{Molecule()[atomIndex],
				                                                  AtomicSubshellFromOrbital(orbitalIndex, atomIndex)};
			}

			std::pair<Eigen::Index, ElectronicSubshell> AtomIndexAndSubshellFromOrbital(
			        const Eigen::Index orbitalIndex) const noexcept
			{
				const auto atomIndex = AtomIndexFromOrbital(orbitalIndex);
				return std::pair{atomIndex, AtomicSubshellFromOrbital(orbitalIndex, atomIndex)};
			}

			// Returned Eigen views are valid only as long as the input vector lives.
			template <typename VectorLike>
			auto OrbitalSegmentOf(Eigen::MatrixBase<VectorLike>& vector, const Atom& atom) const
			{
				auto desc = OrbitalSegmentDescriptionOf(cr_BasisSet.m_Molecule.IndexOf(atom));
				return vector.segment(desc.First, desc.Count);
			}

			// Returned Eigen views are valid only as long as the input vector lives.
			template <typename VectorLike>
			auto OrbitalSegmentOf(const Eigen::MatrixBase<VectorLike>& vector, const Atom& atom) const
			{
				auto desc = OrbitalSegmentDescriptionOf(cr_BasisSet.m_Molecule.IndexOf(atom));
				return vector.segment(desc.First, desc.Count);
			}

			// Returned Eigen object will out-live the input vector
			template <typename VectorLike>
			auto OrbitalSegmentOf(Eigen::MatrixBase<VectorLike>&& vector, const Atom& atom) const
			{
				auto desc = OrbitalSegmentDescriptionOf(cr_BasisSet.m_Molecule.IndexOf(atom));
				return vector.segment(desc.First, desc.Count).eval();
			}

			// Returned Eigen views are valid only as long as the input vector lives.
			template <typename VectorLike>
			auto OrbitalSegmentOf(Eigen::MatrixBase<VectorLike>& vector,
			                      const Atom& atom,
			                      const ElectronicSubshell shell) const
			{
				auto desc = OrbitalSegmentDescriptionOf(cr_BasisSet.m_Molecule.IndexOf(atom), shell);
				return vector.segment(desc.First, desc.Count);
			}

			// Returned Eigen views are valid only as long as the input vector lives.
			template <typename VectorLike>
			auto OrbitalSegmentOf(const Eigen::MatrixBase<VectorLike>& vector,
			                      const Atom& atom,
			                      const ElectronicSubshell shell) const
			{
				auto desc = OrbitalSegmentDescriptionOf(cr_BasisSet.m_Molecule.IndexOf(atom), shell);
				return vector.segment(desc.First, desc.Count);
			}

			// Returned Eigen object will out-live the input vector
			template <typename VectorLike>
			auto OrbitalSegmentOf(Eigen::MatrixBase<VectorLike>&& vector,
			                      const Atom& atom,
			                      const ElectronicSubshell shell) const
			{
				auto desc = OrbitalSegmentDescriptionOf(cr_BasisSet.m_Molecule.IndexOf(atom), shell);
				return vector.segment(desc.First, desc.Count).eval();
			}


		private:
			ElectronicSubshell AtomicSubshellFromOrbital(const Eigen::Index orbitalIndex,
			                                             const Eigen::Index atomIndex) const noexcept
			{
				const auto atomicOffset = cr_BasisSet.OrbitalSegmentationTable<C, R>()[atomIndex];
				const ElementaryBasisSet* basisPtr = cr_BasisSet.m_BasisAssignments[atomIndex];

				const auto& subshellSegTable = cr_BasisSet.m_ComputedElementaryBasisInfoTable.at(basisPtr)
				                                       .template SubshellSegmentationTable<C, R>();
				assert(subshellSegTable.size() >= 2);
				const auto segIterator = std::prev(std::upper_bound(
				        std::next(subshellSegTable.cbegin()), subshellSegTable.cend(), orbitalIndex - atomicOffset));
				const AzimuthalQuantumNumber angularMomentum{
				        static_cast<int>(std::distance(subshellSegTable.cbegin(), segIterator))};
				const auto magneticQuantumNumberCount = R == Representation::Spherical
				                                                ? angularMomentum.MagneticQuantumNumberCount()
				                                                : angularMomentum.CartesianMagneticQuantumNumberCount();
				const auto principalQuantumNumber =
				        angularMomentum.MinPrincipalQuantumNumber()
				        + static_cast<int>(orbitalIndex - *segIterator - atomicOffset) / magneticQuantumNumberCount;

				return ElectronicSubshell{principalQuantumNumber, angularMomentum};
			}

			OrbitalSegmentDescription OrbitalSegmentDescriptionOf(const std::size_t atomIndex) const
			{
				return {cr_BasisSet.OrbitalSegmentationTable<C, R>()[atomIndex],
				        cr_BasisSet.OrbitalSegmentationTable<C, R>()[atomIndex + 1]
				                - cr_BasisSet.OrbitalSegmentationTable<C, R>()[atomIndex]};
			}

			OrbitalSegmentDescription OrbitalSegmentDescriptionOf(const std::size_t atomIndex,
			                                                      const ElectronicSubshell shell) const
			{
				const auto atomicOffset = cr_BasisSet.OrbitalSegmentationTable<C, R>()[atomIndex];
				const ElementaryBasisSet* basisPtr = cr_BasisSet.m_BasisAssignments[atomIndex];
				const auto azimuthalShellOffset =
				        cr_BasisSet.m_ComputedElementaryBasisInfoTable.at(basisPtr)
				                .template SubshellSegmentationTable<C, R>()[shell.AzimuthalQuantumNumber().Value()];
				const auto sphericalMagneticQuantumNumberCount = R == Representation::Spherical
				                                                         ? shell.MagneticQuantumNumberCount()
				                                                         : shell.CartesianMagneticQuantumNumberCount();
				const auto offset = atomicOffset + azimuthalShellOffset
				                    + sphericalMagneticQuantumNumberCount
				                              * (shell.PrincipalQuantumNumber()
				                                 - shell.AzimuthalQuantumNumber().MinPrincipalQuantumNumber());

				return {offset, sphericalMagneticQuantumNumberCount};
			}

		private:
			const Gaussian::MolecularBasisSet& cr_BasisSet;
		};

		template <>
		class BasisView<Contraction::Primitive, Representation::Unspecified>
		{
		public:
			constexpr explicit BasisView(const Gaussian::MolecularBasisSet& basis) noexcept : cr_BasisSet(basis)
			{
				/* NO CODE */
			}

			auto Cartesian() const
			{
				return BasisView<Contraction::Primitive, Representation::Cartesian>{cr_BasisSet};
			}

			auto Spherical() const
			{
				return BasisView<Contraction::Primitive, Representation::Spherical>{cr_BasisSet};
			}

			Eigen::Index SubshellCount() const noexcept
			{
				// we need to remove the const-ness for the limitation of range-v3 and C++20 range library.
				// we won't and shan't modify the map, but it won't compile unless const was removed
				return std::accumulate(
				        cr_BasisSet.m_ComputedElementaryBasisInfoTable.cbegin(),
				        cr_BasisSet.m_ComputedElementaryBasisInfoTable.cend(),
				        Eigen::Index{0},
				        [](const Eigen::Index acc, const auto& kv)
				        {
					        const auto& info = std::get<1>(kv);
					        return acc + info.ReferenceCount * info.template SubshellCount<Contraction::Primitive>();
				        });
			}

			auto SubshellsOf(const Atom& atom) const
			{
				return cr_BasisSet.ElementaryBasisAt(cr_BasisSet.m_Molecule.IndexOf(atom)).AzimuthalShells
				       | ranges::views::transform([](const AzimuthalShell& amb) { return amb.PrimitiveShells(); })
				       | ranges::views::join;
			}

			auto EcpOffsettedSubshellsOf(const Atom& atom) const
			{
				const auto* basisPtr = cr_BasisSet.m_BasisAssignments[cr_BasisSet.m_Molecule.IndexOf(atom)];
				const auto& ambs = basisPtr->AzimuthalShells;

				const auto* ambsPtr = ambs.data();
				const auto* principalQuantumNumberOffsetsPtr =
				        cr_BasisSet.m_ComputedElementaryBasisInfoTable.at(basisPtr)
				                .PrincipalQuantumNumberOffsetTable.data();

				return ambs
				       | ranges::views::transform(
				               [principalQuantumNumberOffsetsPtr, ambsPtr](const AzimuthalShell& amb)
				               { return amb.PrimitiveShells(principalQuantumNumberOffsetsPtr[&amb - ambsPtr]); })
				       | ranges::views::join;
			}


		private:
			const Gaussian::MolecularBasisSet& cr_BasisSet;
		};

		template <>
		class BasisView<Contraction::Contracted, Representation::Unspecified>
		{
		public:
			constexpr explicit BasisView(const Gaussian::MolecularBasisSet& basis) noexcept : cr_BasisSet(basis)
			{
				/* NO CODE */
			}

			auto Cartesian() const
			{
				return BasisView<Contraction::Contracted, Representation::Cartesian>{cr_BasisSet};
			}

			auto Spherical() const
			{
				return BasisView<Contraction::Contracted, Representation::Spherical>{cr_BasisSet};
			}

			Eigen::Index SubshellCount() const noexcept
			{
				// we need to remove the const-ness for the limitation of range-v3 and C++20 range library
				// we won't and shan't modify the map, but it won't compile unless const was removed
				return std::accumulate(
				        cr_BasisSet.m_ComputedElementaryBasisInfoTable.cbegin(),
				        cr_BasisSet.m_ComputedElementaryBasisInfoTable.cend(),
				        Eigen::Index{0},
				        [](const Eigen::Index acc, const auto& kv)
				        {
					        const auto& info = std::get<1>(kv);
					        return acc + info.ReferenceCount * info.template SubshellCount<Contraction::Contracted>();
				        });
			}

			auto SubshellsOf(const Atom& atom) const
			{
				return cr_BasisSet.ElementaryBasisAt(cr_BasisSet.m_Molecule.IndexOf(atom)).AzimuthalShells
				       | ranges::views::transform([](const AzimuthalShell& amb) { return amb.ContractedShells(); })
				       | ranges::views::join;
			}

			auto EcpOffsettedSubshellsOf(const Atom& atom) const
			{
				const auto* basisPtr = cr_BasisSet.m_BasisAssignments[cr_BasisSet.m_Molecule.IndexOf(atom)];
				const auto& ambs = basisPtr->AzimuthalShells;

				const auto* ambsPtr = ambs.data();
				const auto* principalQuantumNumberOffsetsPtr =
				        cr_BasisSet.m_ComputedElementaryBasisInfoTable.at(basisPtr)
				                .PrincipalQuantumNumberOffsetTable.data();

				return ambs
				       | ranges::views::transform(
				               [principalQuantumNumberOffsetsPtr, ambsPtr](const AzimuthalShell& amb)
				               { return amb.ContractedShells(principalQuantumNumberOffsetsPtr[&amb - ambsPtr]); })
				       | ranges::views::join;
			}


		private:
			const Gaussian::MolecularBasisSet& cr_BasisSet;
		};

		template <>
		class BasisView<Contraction::Unspecified, Representation::Cartesian>
		{
		public:
			constexpr explicit BasisView(const Gaussian::MolecularBasisSet& basis) noexcept : cr_BasisSet(basis)
			{
				/* NO CODE */
			}

			auto Primitive() const
			{
				return BasisView<Contraction::Primitive, Representation::Cartesian>{cr_BasisSet};
			}

			auto Contracted() const
			{
				return BasisView<Contraction::Contracted, Representation::Cartesian>{cr_BasisSet};
			}

		private:
			const Gaussian::MolecularBasisSet& cr_BasisSet;
		};

		template <>
		class BasisView<Contraction::Unspecified, Representation::Spherical>
		{
		public:
			constexpr explicit BasisView(const Gaussian::MolecularBasisSet& basis) noexcept : cr_BasisSet(basis)
			{
				/* NO CODE */
			}

			auto Primitive() const
			{
				return BasisView<Contraction::Primitive, Representation::Spherical>{cr_BasisSet};
			}

			auto Contracted() const
			{
				return BasisView<Contraction::Contracted, Representation::Spherical>{cr_BasisSet};
			}

		private:
			const Gaussian::MolecularBasisSet& cr_BasisSet;
		};

		template <>
		class BasisView<Contraction::Unspecified, Representation::Unspecified>
		{
		public:
			constexpr explicit BasisView(const Gaussian::MolecularBasisSet& basis) noexcept : cr_BasisSet(basis)
			{
				/* NO CODE */
			}

			auto Primitive() const
			{
				return BasisView<Contraction::Primitive, Representation::Unspecified>{cr_BasisSet};
			}

			auto Contracted() const
			{
				return BasisView<Contraction::Contracted, Representation::Unspecified>{cr_BasisSet};
			}

			auto Cartesian() const
			{
				return BasisView<Contraction::Unspecified, Representation::Cartesian>{cr_BasisSet};
			}

			auto Spherical() const
			{
				return BasisView<Contraction::Unspecified, Representation::Spherical>{cr_BasisSet};
			}

		private:
			const Gaussian::MolecularBasisSet& cr_BasisSet;
		};
	}  // namespace Detail::MolecularBasisSet


	inline Detail::MolecularBasisSet::BasisView<Detail::MolecularBasisSet::Contraction::Primitive,
	                                            Detail::MolecularBasisSet::Representation::Unspecified>
	MolecularBasisSet::Primitive() const
	{
		return Detail::MolecularBasisSet::BasisView<Detail::MolecularBasisSet::Contraction::Primitive,
		                                            Detail::MolecularBasisSet::Representation::Unspecified>{*this};
	}

	inline Detail::MolecularBasisSet::BasisView<Detail::MolecularBasisSet::Contraction::Contracted,
	                                            Detail::MolecularBasisSet::Representation::Unspecified>
	MolecularBasisSet::Contracted() const
	{
		return Detail::MolecularBasisSet::BasisView<Detail::MolecularBasisSet::Contraction::Contracted,
		                                            Detail::MolecularBasisSet::Representation::Unspecified>{*this};
	}

	inline Detail::MolecularBasisSet::BasisView<Detail::MolecularBasisSet::Contraction::Unspecified,
	                                            Detail::MolecularBasisSet::Representation::Cartesian>
	MolecularBasisSet::Cartesian() const
	{
		return Detail::MolecularBasisSet::BasisView<Detail::MolecularBasisSet::Contraction::Unspecified,
		                                            Detail::MolecularBasisSet::Representation::Cartesian>{*this};
	}

	inline Detail::MolecularBasisSet::BasisView<Detail::MolecularBasisSet::Contraction::Unspecified,
	                                            Detail::MolecularBasisSet::Representation::Spherical>
	MolecularBasisSet::Spherical() const
	{
		return Detail::MolecularBasisSet::BasisView<Detail::MolecularBasisSet::Contraction::Unspecified,
		                                            Detail::MolecularBasisSet::Representation::Spherical>{*this};
	}
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
