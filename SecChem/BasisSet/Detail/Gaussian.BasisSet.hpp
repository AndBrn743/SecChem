// SPDX-License-Identifier: MIT
// Copyright (c) 2025-2026 Andy Brown

#pragma once
#include "Gaussian.SemiLocalEcp.hpp"


#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/algorithm/mismatch.hpp>
#if !defined(SECCHEM_GAUSSIAN_BASIS_SET_INTERNAL)
#error Do not include internal header files directly
#endif

#include <range/v3/algorithm/all_of.hpp>
#include <range/v3/algorithm/find_if.hpp>
#include <range/v3/algorithm/find_if_not.hpp>
#include <range/v3/algorithm/is_sorted.hpp>
#include <range/v3/algorithm/none_of.hpp>
#include <range/v3/algorithm/remove_if.hpp>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/numeric/accumulate.hpp>

namespace SecChem::BasisSet::Gaussian
{
	class ElementaryBasisSet
	{
	private:
		template <typename T>
		static auto ConcatSetsOf(const std::vector<T>& sortedNonEmptyBlocks)
		{
			using Iterator = typename std::vector<T>::const_iterator;
			std::vector<std::pair<Iterator, Iterator>> concatSets;
			auto it0 = sortedNonEmptyBlocks.begin();
			while (it0 != sortedNonEmptyBlocks.end())
			{
				auto it1 = std::find_if_not(std::next(it0),
				                            sortedNonEmptyBlocks.end(),
				                            [l0 = it0->AngularMomentum()](const T& block)
				                            { return block.AngularMomentum() == l0; });

				if (it1 != std::next(it0))
				{
					concatSets.push_back({it0, it1});
					it0 = it1;
				}
				else
				{
					++it0;
				}
			}

			return concatSets;
		}

		static bool IsInStandardStorageOrder_OverloadSet(const AzimuthalShell& lhs,
		                                                 const AzimuthalShell& rhs)
		{
			if (lhs.AngularMomentum() != rhs.AngularMomentum())
			{
				return lhs.AngularMomentum() < rhs.AngularMomentum();
			}

			//----------------------

			if (lhs.IsNotEmpty() != rhs.IsNotEmpty())
			{
				return lhs.IsNotEmpty();
			}

			if (lhs.IsNotEmpty())
			{
				if (lhs.ExponentSet().size() != rhs.ExponentSet().size())
				{
					return lhs.ExponentSet().size() > rhs.ExponentSet().size();
				}

				return lhs.ExponentSet()[0] > rhs.ExponentSet()[0];
			}

			return false;
		}

		static bool IsInStandardStorageOrder_OverloadSet(const SemiLocalEcpChannel& lhs, const SemiLocalEcpChannel& rhs)
		{
			if (lhs.AngularMomentum() != rhs.AngularMomentum())
			{
				return lhs.AngularMomentum() < rhs.AngularMomentum();
			}

			if (lhs.GaussianTermCount() != rhs.GaussianTermCount())
			{
				return lhs.GaussianTermCount() > rhs.GaussianTermCount();
			}

			if (lhs.GaussianExponent(0) != rhs.GaussianExponent(0))
			{
				return lhs.GaussianExponent(0) > rhs.GaussianExponent(0);
			}

			if (lhs.Coefficient(0) != rhs.Coefficient(0))
			{
				return lhs.Coefficient(0) > rhs.Coefficient(0);
			}

			if (lhs.RExponent(0) != rhs.RExponent(0))
			{
				return lhs.RExponent(0) > rhs.RExponent(0);
			}

			assert(lhs == rhs);

			return true;
		}

	public:
		std::vector<AzimuthalShell> AzimuthalShells{};
		std::vector<SemiLocalEcpChannel> SemiLocalEcpChannels{};
		int EcpElectronCount{};

		static constexpr auto IsInStandardStorageOrder = [](const auto& lhs, const auto& rhs)
		{ return IsInStandardStorageOrder_OverloadSet(lhs, rhs); };

		bool IsInStandardRepresentation() const noexcept
		{
			if (!AzimuthalShells.empty()
			    && (ranges::any_of(AzimuthalShells, &AzimuthalShell::IsEmpty)
			        || !ranges::is_sorted(AzimuthalShells, IsInStandardStorageOrder)
			        || !ConcatSetsOf(AzimuthalShells).empty()))
			{
				return false;
			}

			if (!SemiLocalEcpChannels.empty()
			    && (ranges::any_of(SemiLocalEcpChannels, &SemiLocalEcpChannel::IsEmpty)
			        || !ranges::is_sorted(SemiLocalEcpChannels, IsInStandardStorageOrder)
			        || !ConcatSetsOf(SemiLocalEcpChannels).empty()))
			{
				return false;
			}

			return true;
		}

		void StandardizeRepresentation()
		{
			AzimuthalShells = StandardizedRepresentationOf(
			        AzimuthalShells, [](const AzimuthalShell& s) { return s.IsEmpty(); });

			SemiLocalEcpChannels = StandardizedRepresentationOf(
			        SemiLocalEcpChannels, [](const SemiLocalEcpChannel& p) { return p.IsEmpty(); });
		}

		template <typename T, typename IsEmpty>
		static std::vector<T> StandardizedRepresentationOf(std::vector<T> blocks, const IsEmpty isEmpty)
		{
			if (blocks.empty() || ranges::all_of(blocks, isEmpty))
			{
				return {};
			}

			if (blocks.size() == 1)
			{
				return isEmpty(blocks[0]) ? std::vector<T>{} : blocks;
			}

			blocks.erase(ranges::remove_if(blocks, isEmpty), blocks.end());
			ranges::sort(blocks, IsInStandardStorageOrder);

			const auto concatSets = ConcatSetsOf(blocks);
			if (concatSets.empty())
			{
				return blocks;
			}

			const auto concatenatedAmbCount =
			        blocks.size()
			        - ranges::accumulate(concatSets,
			                             Eigen::Index{0},
			                             [](const Eigen::Index acc, const auto& iteratorPair)
			                             { return acc + std::distance(iteratorPair.first, iteratorPair.second) - 1; });

			std::vector<T> concatenatedBlocks;
			concatenatedBlocks.reserve(concatenatedAmbCount);

			auto ambIterator = blocks.cbegin();
			auto concatSetIterator = concatSets.cbegin();
			while (ambIterator != blocks.cend())
			{
				if (concatSetIterator != concatSets.cend() && ambIterator == concatSetIterator->first)
				{
					concatenatedBlocks.emplace_back(T::Concat(concatSetIterator->first, concatSetIterator->second));
					ambIterator = concatSetIterator->second;
					++concatSetIterator;
				}
				else
				{
					concatenatedBlocks.emplace_back(*ambIterator);
					++ambIterator;
				}
			}

			return concatenatedBlocks;
		}

		bool EqualsTo(const ElementaryBasisSet& other, const Scalar tolerance = 1e-15) const noexcept
		{
			if (AzimuthalShells.size() != other.AzimuthalShells.size()
			    || SemiLocalEcpChannels.size() != other.SemiLocalEcpChannels.size())
			{
				return false;
			}

			return ranges::mismatch(AzimuthalShells,
			                        other.AzimuthalShells,
			                        [tolerance](const auto& lhs, const auto& rhs)
			                        { return lhs.EqualsTo(rhs, tolerance); })
			               .in1
			       == AzimuthalShells.cend();
		}

		bool NotEqualsTo(const ElementaryBasisSet& other, const Scalar tolerance = 1e-15) const noexcept
		{
			return !EqualsTo(other, tolerance);
		}

		bool operator==(const ElementaryBasisSet& other) const
		{
			return EqualsTo(other, 0);
		}

		bool operator!=(const ElementaryBasisSet& other) const
		{
			return !EqualsTo(other, 0);
		}

		/// <c>EcpElectronCount</c> must be 2, 4, 10, 12, 18, 28, 30, 36, 46, 60, 62, 68, 78, 92, or 110.
		/// the behavior is undefined otherwise
		std::vector<int> CreatePrincipalQuantumNumberOffsetTable() const
		{
			std::vector offsets(AzimuthalShells.back().AngularMomentum().Value() + 1, 0);

			if (EcpElectronCount == 0)
			{
				return offsets;
			}

			if (!offsets.empty())
			{
				if (EcpElectronCount >= 2)
				{
					offsets[0] = 1;
				}
				if (EcpElectronCount >= 4)
				{
					offsets[0] = 2;
				}
				if (EcpElectronCount >= 12)
				{
					offsets[0] = 3;
				}
				if (EcpElectronCount >= 30)
				{
					offsets[0] = 4;
				}
				if (EcpElectronCount >= 62)
				{
					offsets[0] = 5;
				}
			}
			if (offsets.size() > 1)
			{
				if (EcpElectronCount >= 10)
				{
					offsets[1] = 1;
				}
				if (EcpElectronCount >= 18)
				{
					offsets[1] = 2;
				}
				if (EcpElectronCount >= 36)
				{
					offsets[1] = 3;
				}
				if (EcpElectronCount >= 68)
				{
					offsets[1] = 4;
				}
			}
			else if (offsets.size() > 2)
			{
				if (EcpElectronCount >= 28)
				{
					offsets[2] = 1;
				}
				if (EcpElectronCount >= 60)
				{
					offsets[2] = 2;
				}
				if (EcpElectronCount >= 78)
				{
					offsets[2] = 3;
				}
			}
			else if (offsets.size() > 3)
			{
				if (EcpElectronCount >= 60)
				{
					offsets[3] = 1;
				}
				if (EcpElectronCount >= 92)
				{
					offsets[3] = 2;
				}
			}
			else if (offsets.size() > 4)
			{
				if (EcpElectronCount >= 110)
				{
					offsets[4] = 1;
				}
			}

			return offsets;
		}
	};


	namespace Detail
	{
		template <OwnershipSemantics Semantics>
		class BasisSetImpl
		{
			static constexpr auto AlternativeOwnershipSemantics =
			        static_cast<OwnershipSemantics>(!static_cast<bool>(Semantics));
			using DataType = std::unordered_map<Element, ElementaryBasisSet>;
			using StorageType =
			        std::conditional_t<Semantics == OwnershipSemantics::Value, DataType, std::shared_ptr<DataType>>;
			friend BasisSetImpl<OwnershipSemantics::Reference>;
			friend BasisSetImpl<OwnershipSemantics::Value>;


		public:
			BasisSetImpl() : m_DataStorage(CreateStorage())
			{
				/* NO CODE */
			}

			/// remarks: This conversion constructor will always deep-copy the operand
			explicit BasisSetImpl(const BasisSetImpl<AlternativeOf(Semantics)>& other)
			{
				if constexpr (Semantics == OwnershipSemantics::Value)
				{
					m_DataStorage = other.Data();
				}
				else
				{
					m_DataStorage = std::make_shared<DataType>(other.Data());
				}
			}

			explicit BasisSetImpl(BasisSetImpl<AlternativeOf(Semantics)>&& other)
			{
				static_assert(AlternativeOf(Semantics) == OwnershipSemantics::Value);
				m_DataStorage = std::make_shared<DataType>(std::move(other.Data()));
			}

			template <typename..., OwnershipSemantics SemanticsToo = Semantics>  // poor man's `requires`
			std::enable_if_t<SemanticsToo == OwnershipSemantics::Reference, BasisSetImpl> Clone() const
			{
				return BasisSetImpl{std::make_shared<DataType>(*m_DataStorage)};
			}

			const ElementaryBasisSet& operator[](const Element element) const
			{
				return Data().at(element);
			}

			ElementaryBasisSet& operator[](const Element element)
			{
				return Data().at(element);
			}

			ElementaryBasisSet& AddOrOverwriteEntryOf(const Element element)
			{
				Data()[element] = {};
				return Data()[element];
			}

			ElementaryBasisSet& AddEntryFor(const Element element)
			{
				if (Data().find(element) != Data().end())
				{
					throw std::logic_error(
					        "Duplicated entry from same element detected. This can be an indication of a bug in the "
					        "input file or input parser. If you intend for an overwrite, please use "
					        "OverwriteEntryOf(...) or AddOrOverwriteEntryOf(...) instead.");
				}
				return Data()[element];
			}

			ElementaryBasisSet& OverwriteEntryOf(const Element element)
			{
				Data().at(element) = {};
				return Data()[element];
			}

			bool Has(const Element element) const
			{
				return Data().find(element) != Data().end();
			}

			/// Value comparison will be done for BasisSet, reference comparison will be done for SharedBasisSet.
			/// To compare the value of SharedBasisSet, one should use <c>EqualsTo</c> method instead.
			/// Please note that value comparison is depended on the underlying storage representation.
			template <OwnershipSemantics OtherSemantics>
			constexpr bool operator==(const BasisSetImpl<OtherSemantics>& other) const noexcept
			{
				if constexpr (Semantics != OtherSemantics)
				{
					return false;
				}
				else if constexpr (Semantics == OwnershipSemantics::Reference)
				{
					return m_DataStorage == other.m_DataStorage;
				}
				else
				{
					return EqualsTo(other, 0);
				}
			}

			/// Value comparison will be done for BasisSet, reference comparison will be done for SharedBasisSet.
			/// To compare the value of SharedBasisSet, one should use <c>NotEqualsTo</c> method instead
			/// Please note that value comparison is depended on the underlying storage representation.
			template <OwnershipSemantics OtherSemantics>
			constexpr bool operator!=(const BasisSetImpl<OtherSemantics>& other) const noexcept
			{
				return !(*this == other);
			}

			template <OwnershipSemantics OtherSemantics>
			bool EqualsTo(const BasisSetImpl<OtherSemantics>& other, const Scalar tolerance = 1e-15) const noexcept
			{
				if constexpr (OtherSemantics == Semantics)
				{
					if (this == &other)
					{
						return true;
					}
				}

				if (Data().size() != other.Data().size())
				{
					return false;
				}

				for (const auto& [element, elementaryBasis] : Data())
				{
					const auto iteratorToOtherElementaryBasis = other.Data().find(element);
					if (iteratorToOtherElementaryBasis == other.Data().end())
					{
						return false;
					}

					if (elementaryBasis.NotEqualsTo(iteratorToOtherElementaryBasis->second, tolerance))
					{
						return false;
					}
				}

				return true;
			}

			template <OwnershipSemantics OtherSemantics>
			bool NotEqualsTo(const BasisSetImpl<OtherSemantics>& other, const Scalar tolerance = 1e-15) const noexcept
			{
				return !EqualsTo(other, tolerance);
			}

			bool IsInStandardRepresentation() const noexcept
			{
				if (Data().empty())
				{
					return true;
				}

				for (const auto& [_, elementaryBasis] : Data())
				{
					if (!elementaryBasis.IsInStandardRepresentation()
					    || (elementaryBasis.AzimuthalShells.empty()
					        && elementaryBasis.SemiLocalEcpChannels.empty()))
					{
						return false;
					}
				}

				return true;
			}

			void StandardizeRepresentation()
			{
				std::vector<Element> nullElements;

				for (auto& [element, elementaryBasis] : Data())
				{
					const auto& angularMomentumBlocks = elementaryBasis.AzimuthalShells;
					const auto& ecpProjectors = elementaryBasis.SemiLocalEcpChannels;
					if ((angularMomentumBlocks.empty()
					     || ranges::all_of(angularMomentumBlocks, &AzimuthalShell::IsEmpty))
					    && (ecpProjectors.empty() || ranges::all_of(ecpProjectors, &SemiLocalEcpChannel::IsEmpty)))
					{
						nullElements.emplace_back(element);
					}

					elementaryBasis.StandardizeRepresentation();
				}

				for (const auto element : nullElements)
				{
					Data().erase(element);
				}
			}

			BasisSetImpl ToStandardizedRepresentation() const
			{
				if constexpr (Semantics == OwnershipSemantics::Reference)
				{
					BasisSetImpl result = Clone();
					result.StandardizeRepresentation();
					return result;
				}
				else
				{
					BasisSetImpl result = *this;
					result.StandardizeRepresentation();
					return result;
				}
			}

			auto begin() noexcept
			{
				return Data().begin();
			}

			auto end() noexcept
			{
				return Data().end();
			}

			auto begin() const noexcept
			{
				return Data().begin();
			}

			auto end() const noexcept
			{
				return Data().end();
			}

			auto cbegin() const noexcept
			{
				return Data().cbegin();
			}

			auto cend() const noexcept
			{
				return Data().cend();
			}


		private:
			explicit BasisSetImpl(StorageType&& storage) : m_DataStorage(std::move(storage))
			{
				/* NO CODE */
			}

			static StorageType CreateStorage()
			{
				if constexpr (Semantics == OwnershipSemantics::Reference)
				{
					return std::make_shared<DataType>();
				}
				else
				{
					return {};
				}
			}

			const DataType& Data() const noexcept
			{
				if constexpr (Semantics == OwnershipSemantics::Reference)
				{
					return *m_DataStorage;
				}
				else
				{
					return m_DataStorage;
				}
			}

			DataType& Data() noexcept
			{
				if constexpr (Semantics == OwnershipSemantics::Reference)
				{
					return *m_DataStorage;
				}
				else
				{
					return m_DataStorage;
				}
			}


		private:
			StorageType m_DataStorage{};
		};
	}  // namespace Detail

	using BasisSet = Detail::BasisSetImpl<OwnershipSemantics::Value>;
	using SharedBasisSet = Detail::BasisSetImpl<OwnershipSemantics::Reference>;

	namespace Detail
	{
		template <OwnershipSemantics Semantics>
		class BasisSetLibraryImpl
		{
			friend Builder<BasisSetLibraryImpl>;
			using DataStorage = std::unordered_map<std::string, BasisSetImpl<Semantics>>;

		public:
			BasisSetLibraryImpl Clone() const
			{
				static_assert(Semantics == OwnershipSemantics::Reference);

				BasisSetLibraryImpl result;

				for (const auto& [name, basisData] : *this)
				{
					result.m_Data[name] = basisData.Clone();
				}

				return result;
			}

			const BasisSetImpl<Semantics>& operator[](const std::string& basisSetName) const
			{
				return m_Data.at(basisSetName);
			}

			// One cannot mutate basis sets inside an already built library, even if the library was mutable.
			// This is to ensure stable iterators and references into each basis sets of the library
			const BasisSetImpl<Semantics>& operator[](const std::string& basisSetName)
			{
				return m_Data.at(basisSetName);
			}

			bool Has(const std::string& basisSetName) const noexcept
			{
				// we need to support C++17
				// ReSharper disable once CppUseAssociativeContains
				return m_Data.find(basisSetName) != m_Data.end();
			}


		private:
			DataStorage m_Data = {};
		};
	}  // namespace Detail

	using BasisSetLibrary = Detail::BasisSetLibraryImpl<OwnershipSemantics::Value>;
	using SharedBasisSetLibrary = Detail::BasisSetLibraryImpl<OwnershipSemantics::Reference>;
}  // namespace SecChem::BasisSet::Gaussian

template <SecChem::OwnershipSemantics Semantics>
class SecChem::Builder<SecChem::BasisSet::Gaussian::Detail::BasisSetLibraryImpl<Semantics>>
{
	using Product = BasisSet::Gaussian::Detail::BasisSetLibraryImpl<Semantics>;

public:
	Builder() : m_ProductPtr(std::make_unique<Product>())
	{
		/* NO CODE */
	}

	Product Build()
	{
		if (m_ProductPtr == nullptr)
		{
			throw std::logic_error("Builder already consumed");
		}

		auto h = std::move(*m_ProductPtr);
		m_ProductPtr.reset();
		return h;
	}

	Builder& AddBasisSet(std::string name, BasisSet::Gaussian::Detail::BasisSetImpl<Semantics> definition)
	{
		if (m_ProductPtr == nullptr)
		{
			throw std::logic_error("Builder already consumed");
		}

		if (m_ProductPtr->m_Data.find(name) != m_ProductPtr->m_Data.end())
		{
			throw std::logic_error("A basis set with the same name already exists");
		}

		m_ProductPtr->m_Data.emplace(std::move(name), std::move(definition));
		return *this;
	}

private:
	std::unique_ptr<Product> m_ProductPtr = nullptr;
};
