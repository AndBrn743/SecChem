#pragma once
#if !defined(SECCHEM_GAUSSIAN_BASIS_SET_INTERNAL)
#error Do not include internal header files directly
#endif

namespace SecChem::BasisSet::Gaussian
{
	namespace Detail
	{
		template <OwnershipSemantics Semantics>
		class BasisSetImpl
		{
			static constexpr auto AlternativeOwnershipSemantics =
			        static_cast<OwnershipSemantics>(!static_cast<bool>(Semantics));
			using DataType = std::unordered_map<Element, std::vector<AngularMomentumBlock>>;
			using StorageType =
			        std::conditional_t<Semantics == OwnershipSemantics::Value, DataType, std::shared_ptr<DataType>>;
			friend BasisSetImpl<OwnershipSemantics::Reference>;
			friend BasisSetImpl<OwnershipSemantics::Value>;

			static constexpr auto IsInStandardStorageOrder =
			        [](const AngularMomentumBlock& lhs, const AngularMomentumBlock& rhs)
			{
				if (lhs.AngularMomentum() != rhs.AngularMomentum())
				{
					return lhs.AngularMomentum() < rhs.AngularMomentum();
				}

				//----------------------

				if (lhs.HasOrbital() != rhs.HasOrbital())
				{
					return lhs.HasOrbital();
				}

				if (lhs.HasOrbital())
				{
					if (lhs.ExponentSet().size() != rhs.ExponentSet().size())
					{
						return lhs.ExponentSet().size() > rhs.ExponentSet().size();
					}

					return lhs.ExponentSet()[0] > rhs.ExponentSet()[0];
				}

				//----------------------

				if (lhs.HasSemiLocalEcp() != rhs.HasSemiLocalEcp())
				{
					return lhs.HasSemiLocalEcp();
				}

				if (lhs.HasSemiLocalEcp())
				{
					if (lhs.SemiLocalEcp().TermCount() != rhs.SemiLocalEcp().TermCount())
					{
						return lhs.SemiLocalEcp().TermCount() > rhs.SemiLocalEcp().TermCount();
					}

					return lhs.SemiLocalEcp().GaussianExponent(0) > rhs.SemiLocalEcp().GaussianExponent(0);
				}

				return false;
			};

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

			const std::vector<AngularMomentumBlock>& operator[](const Element element) const
			{
				return Data().at(element);
			}

			std::vector<AngularMomentumBlock>& operator[](const Element element)
			{
				return Data().at(element);
			}

			std::vector<AngularMomentumBlock>& AddOrOverwriteEntryOf(const Element element)
			{
				Data()[element] = {};
				return Data()[element];
			}

			std::vector<AngularMomentumBlock>& AddEntryFor(const Element element)
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

			std::vector<AngularMomentumBlock>& OverwriteEntryOf(const Element element)
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

				for (const auto& [element, angularMomentumBlocks] : Data())
				{
					const auto iteratorToOtherAngularMomentumBlocks = other.Data().find(element);
					if (iteratorToOtherAngularMomentumBlocks == other.Data().end())
					{
						return false;
					}

					const auto& [_, otherAngularMomentumBlocks] = *iteratorToOtherAngularMomentumBlocks;
					if (angularMomentumBlocks.size() != otherAngularMomentumBlocks.size())
					{
						return false;
					}

					for (auto it0 = angularMomentumBlocks.cbegin(), it1 = otherAngularMomentumBlocks.cbegin();
					     it0 != angularMomentumBlocks.cend();
					     ++it0, ++it1)
					{
						if (it0->NotEqualsTo(*it1, tolerance))
						{
							return false;
						}
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

				for (const auto& [_, angularMomentumBlocks] : Data())
				{
					if (angularMomentumBlocks.empty()
					    || !std::is_sorted(
					            angularMomentumBlocks.cbegin(), angularMomentumBlocks.cend(), IsInStandardStorageOrder)
					    || !AngularMomentumBlockConcatSets(angularMomentumBlocks).empty())
					{
						return false;
					}
				}

				return true;
			}

			void StandardizeRepresentation()
			{
				std::vector<Element> nullElements;

				for (auto& [element, angularMomentumBlocks] : Data())
				{
					if (angularMomentumBlocks.empty()
					    || std::none_of(angularMomentumBlocks.cbegin(),
					                    angularMomentumBlocks.cend(),
					                    [](const AngularMomentumBlock& amb)
					                    { return amb.HasOrbital() || amb.HasSemiLocalEcp(); }))
					{
						nullElements.emplace_back(element);
					}

					if (angularMomentumBlocks.size() <= 1)
					{
						continue;
					}

					std::sort(angularMomentumBlocks.begin(), angularMomentumBlocks.end(), IsInStandardStorageOrder);

					const auto concatSets = AngularMomentumBlockConcatSets(angularMomentumBlocks);
					if (concatSets.empty())
					{
						continue;
					}

					const auto concatenatedAmbCount =
					        angularMomentumBlocks.size()
					        - std::accumulate(
					                concatSets.cbegin(),
					                concatSets.cend(),
					                Eigen::Index{0},
					                [](const Eigen::Index acc, const auto& iteratorPair)
					                { return acc + std::distance(iteratorPair.first, iteratorPair.second) - 1; });

					std::vector<AngularMomentumBlock> concatenatedAngularMomentumBlocks;
					concatenatedAngularMomentumBlocks.reserve(concatenatedAmbCount);

					auto ambIterator = angularMomentumBlocks.cbegin();
					auto concatSetIterator = concatSets.cbegin();
					while (ambIterator != angularMomentumBlocks.cend())
					{
						if (ambIterator == concatSetIterator->first)
						{
							concatenatedAngularMomentumBlocks.emplace_back(
							        AngularMomentumBlock::Concat(concatSetIterator->first, concatSetIterator->second));
							ambIterator = concatSetIterator->second;
							++concatSetIterator;
						}
						else
						{
							concatenatedAngularMomentumBlocks.emplace_back(*ambIterator);
							++ambIterator;
						}
					}

					angularMomentumBlocks = concatenatedAngularMomentumBlocks;
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

			static auto AngularMomentumBlockConcatSets(
			        const std::vector<AngularMomentumBlock>& sortedAngularMomentumBlocks)
			{
				using Iterator = std::vector<AngularMomentumBlock>::const_iterator;
				std::vector<std::pair<Iterator, Iterator>> concatSets;
				auto it0 = sortedAngularMomentumBlocks.begin();
				while (it0 != sortedAngularMomentumBlocks.end())
				{
					auto it1 = std::find_if_not(std::next(it0),
					                            sortedAngularMomentumBlocks.end(),
					                            [l0 = it0->AngularMomentum()](const AngularMomentumBlock& block)
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
		return *std::exchange(m_ProductPtr, std::make_unique<Product>()).release();
	}

	Builder& AddBasisSet(std::string name, BasisSet::Gaussian::Detail::BasisSetImpl<Semantics> definition)
	{
		if (m_ProductPtr->m_Data.find(name) != m_ProductPtr->m_Data.end())
		{
			throw std::logic_error("A basis set with same name already exist in the basis set library");
		}

		m_ProductPtr->m_Data[std::move(name)] = std::move(definition);
		return *this;
	}

private:
	std::unique_ptr<Product> m_ProductPtr = nullptr;
};
