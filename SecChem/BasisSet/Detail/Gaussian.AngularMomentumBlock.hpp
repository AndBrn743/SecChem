#pragma once
#if !defined(SECCHEM_GAUSSIAN_BASIS_SET_INTERNAL)
#error Do not include internal header files directly
#endif

namespace SecChem::BasisSet::Gaussian
{
	template <typename Derived>
	class AbstractAngularMomentumBlock
	    : public SecUtility::IEquatableWithTolerance<AbstractAngularMomentumBlock<Derived>>
	{
		friend SecUtility::IEquatableWithTolerance<AbstractAngularMomentumBlock>;
		friend Derived;

	public:
		const auto& ExponentSet() const noexcept
		{
			return static_cast<const Derived&>(*this).ExponentSet_Impl();
		}

		const auto& ContractionSets() const noexcept
		{
			return static_cast<const Derived&>(*this).ContractionSets_Impl();
		}

		Eigen::Index PrimitiveShellCount() const noexcept
		{
			return ExponentSet().size();
		}

		Eigen::Index ContractedShellCount() const noexcept
		{
			return ContractionSets().cols();
		}

		Eigen::Index PrimitiveCartesianOrbitalCount() const noexcept
		{
			return PrimitiveShellCount() * m_AzimuthalQuantumNumber.CartesianMagneticQuantumNumberCount();
		}

		Eigen::Index PrimitiveSphericalOrbitalCount() const noexcept
		{
			return PrimitiveShellCount() * m_AzimuthalQuantumNumber.MagneticQuantumNumberCount();
		}

		Eigen::Index ContractedCartesianOrbitalCount() const noexcept
		{
			return ContractedShellCount() * m_AzimuthalQuantumNumber.CartesianMagneticQuantumNumberCount();
		}

		Eigen::Index ContractedSphericalOrbitalCount() const noexcept
		{
			return ContractedShellCount() * m_AzimuthalQuantumNumber.MagneticQuantumNumberCount();
		}

		const AzimuthalQuantumNumber& AngularMomentum() const noexcept
		{
			return m_AzimuthalQuantumNumber;
		}

	private:
		const auto& ExponentSet_Impl() const noexcept = delete;

		const auto& ContractionSets_Impl() const noexcept = delete;

		bool EqualsTo_Impl(const AbstractAngularMomentumBlock& other, const Scalar tolerance) const noexcept
		{
			return PrimitiveShellCount() == other.PrimitiveShellCount()
			       && ContractedShellCount() == other.ContractedShellCount()
			       && (ExponentSet() - other.ExponentSet()).cwiseAbs().maxCoeff() <= tolerance
			       && (ContractionSets() - other.ContractionSets()).cwiseAbs().maxCoeff() <= tolerance;
		}

		explicit constexpr AbstractAngularMomentumBlock(const AzimuthalQuantumNumber l) noexcept
		    : m_AzimuthalQuantumNumber(l)
		{
			/* NO CODE */
		}

		constexpr AbstractAngularMomentumBlock(const AbstractAngularMomentumBlock&) noexcept = default;
		constexpr AbstractAngularMomentumBlock(AbstractAngularMomentumBlock&&) noexcept = default;
		constexpr AbstractAngularMomentumBlock& operator=(const AbstractAngularMomentumBlock&) noexcept = default;
		constexpr AbstractAngularMomentumBlock& operator=(AbstractAngularMomentumBlock&&) noexcept = default;

#if __cplusplus >= 202002L
		constexpr
#endif
		        ~AbstractAngularMomentumBlock() noexcept = default;

		AzimuthalQuantumNumber m_AzimuthalQuantumNumber;
	};

	class AngularMomentumBlockSegmentView : public AbstractAngularMomentumBlock<AngularMomentumBlockSegmentView>
	{
		using Base = AbstractAngularMomentumBlock;
		friend Base;
		friend AngularMomentumBlock;
		using Scalar = SecChem::Scalar;

		AngularMomentumBlockSegmentView(const AzimuthalQuantumNumber azimuthalQuantumNumber,
		                                const Eigen::Map<const Eigen::VectorX<Scalar>>& exponentSet,
		                                const Eigen::Block<const Eigen::MatrixX<Scalar>>& contractionSets) noexcept
		    : Base(azimuthalQuantumNumber), m_ExponentSetView(exponentSet), m_ContractionSetsView(contractionSets)
		{
			/* NO CODE */
		}

		const auto& ExponentSet_Impl() const noexcept
		{
			return m_ExponentSetView;
		}

		const auto& ContractionSets_Impl() const noexcept
		{
			return m_ContractionSetsView;
		}

		Eigen::Map<const Eigen::VectorX<SecChem::Scalar>> m_ExponentSetView;
		Eigen::Block<const Eigen::MatrixX<Scalar>> m_ContractionSetsView;
	};

	class AngularMomentumBlock : public SecUtility::IEquatableWithTolerance<AngularMomentumBlock>,
	                             public AbstractAngularMomentumBlock<AngularMomentumBlock>
	{
		friend IEquatableWithTolerance<AngularMomentumBlock>;
		using Base = AbstractAngularMomentumBlock;
		friend Base;
		using SegmentationTable = std::vector<std::pair<Eigen::Index, Eigen::Index>>;

	public:
		AngularMomentumBlock(const AzimuthalQuantumNumber angularMomentum,
		                     ContractedRadialOrbitalSet contractedRadialOrbitalSet)
		    : Base(angularMomentum), m_NullableContractedRadialOrbitalSet(std::move(contractedRadialOrbitalSet)),
		      m_SegmentationTable{{0, 0},
		                          {m_NullableContractedRadialOrbitalSet.value().PrimitiveShellCount(),
		                           m_NullableContractedRadialOrbitalSet.value().ContractedShellCount()}}
		{
			/* NO CODE */
		}

		AngularMomentumBlock(const AzimuthalQuantumNumber angularMomentum,
		                     ContractedRadialOrbitalSet contractedRadialOrbitalSet,
		                     SemiLocalEcp ecp)
		    : Base(angularMomentum), m_NullableContractedRadialOrbitalSet(std::move(contractedRadialOrbitalSet)),
		      m_SegmentationTable{{0, 0},
		                          {m_NullableContractedRadialOrbitalSet.value().PrimitiveShellCount(),
		                           m_NullableContractedRadialOrbitalSet.value().ContractedShellCount()}},
		      m_NullableSemiLocalEcp(std::move(ecp))
		{
			/* NO CODE */
		}

		AngularMomentumBlock(const AzimuthalQuantumNumber angularMomentum, std::nullopt_t, SemiLocalEcp ecp)
		    : Base(angularMomentum), m_NullableSemiLocalEcp(std::move(ecp))
		{
			/* NO CODE */
		}

		AngularMomentumBlock& AddOrOverrideSemiLocalEcp(SemiLocalEcp ecp)
		{
			m_NullableSemiLocalEcp = std::move(ecp);
			return *this;
		}

		AngularMomentumBlock& AddOrOverrideContractedRadialOrbitalSet(
		        ContractedRadialOrbitalSet contractedRadialOrbitalSet)
		{
			m_NullableContractedRadialOrbitalSet = contractedRadialOrbitalSet;
			m_SegmentationTable = {{0, 0},
			                       {m_NullableContractedRadialOrbitalSet.value().PrimitiveShellCount(),
			                        m_NullableContractedRadialOrbitalSet.value().ContractedShellCount()}};
			return *this;
		}

		bool HasOrbital() const noexcept
		{
			return m_NullableContractedRadialOrbitalSet.has_value();
		}

		bool HasSemiLocalEcp() const noexcept
		{
			return m_NullableSemiLocalEcp.has_value();
		}

		// GCC require elaborated type specifier
		// ReSharper disable once CppRedundantElaboratedTypeSpecifier
		const class SemiLocalEcp& SemiLocalEcp() const
		{
			if (!m_NullableSemiLocalEcp.has_value())
			{
				throw std::logic_error("AngularMomentumBlock::SemiLocalEcp() called without a SemiLocalEcp present. "
				                       "Use HasSemiLocalEcp() to check availability.");
			}

			return m_NullableSemiLocalEcp.value();
		}

		// GCC require elaborated type specifier
		// ReSharper disable once CppRedundantElaboratedTypeSpecifier
		const class ContractedRadialOrbitalSet& ContractedRadialOrbitalSet() const noexcept
		{
			assert(m_NullableContractedRadialOrbitalSet.has_value());
			return m_NullableContractedRadialOrbitalSet.value();
		}

		template <typename InputIterator, typename Getter>
		static AngularMomentumBlock Concat(InputIterator begin, const InputIterator end, Getter get)
		{
			if (begin == end)
			{
				throw std::runtime_error("AngularMomentumBlock: input range is empty");
			}

			if (std::distance(begin, end) == 1)
			{
				return get(*begin);
			}

			const auto azimuthalQuantumNumber = get(*begin).m_AzimuthalQuantumNumber;

			if (std::any_of(begin,
			                end,
			                [get, azimuthalQuantumNumber](const auto& item)
			                { return get(item).m_AzimuthalQuantumNumber != azimuthalQuantumNumber; }))
			{
				throw std::runtime_error("AngularMomentumBlock: cannot concat blocks with different angular momenta");
			}

			// clang-format off
			auto crs = ContractedRadialOrbitalSet::ConcatNullable(
					begin, end, [get](const auto& block) -> const auto& { return get(block).m_NullableContractedRadialOrbitalSet; });
			auto ecp = SemiLocalEcp::ConcatNullable(
					begin, end, [get](const auto& block) -> const auto& { return get(block).m_NullableSemiLocalEcp; });
			auto segTable = ConcatSegmentationTable(
					begin, end, [get](const auto& block) -> const auto& { return get(block).m_SegmentationTable; });
			// clang-format on

			return {azimuthalQuantumNumber, std::move(crs), std::move(segTable), std::move(ecp)};
		}

		template <typename InputIterator>
		static AngularMomentumBlock Concat(InputIterator begin, const InputIterator end)
		{
			return Concat(begin, end, [](auto&& item) -> decltype(auto) { return std::forward<decltype(item)>(item); });
		}

#if __cplusplus >= 202002L
		constexpr Eigen::Index SegmentCount() const noexcept
#else
		Eigen::Index SegmentCount() const noexcept
#endif
		{
			assert(HasOrbital());
			return static_cast<Eigen::Index>(m_SegmentationTable.size() - 1);
		}

		AngularMomentumBlockSegmentView Segment(const Eigen::Index index) const noexcept
		{
			assert(HasOrbital());
			assert(index >= 0 && index < SegmentCount());

			const auto [p0, c0] = m_SegmentationTable[index];
			const auto [p1, c1] = m_SegmentationTable[index + 1];

			return AngularMomentumBlockSegmentView{AngularMomentum(),
			                                       {ExponentSet().data() + p0, p1 - p0},
			                                       ContractionSets().block(p0, c0, p1 - p0, c1 - c0)};
		}

		using IEquatableWithTolerance<AngularMomentumBlock>::EqualsTo;
		using IEquatableWithTolerance<AngularMomentumBlock>::NotEqualsTo;
		using IEquatableWithTolerance<AngularMomentumBlock>::operator==;
		using IEquatableWithTolerance<AngularMomentumBlock>::operator!=;

	private:
		AngularMomentumBlock(const AzimuthalQuantumNumber angularMomentum,
		                     std::optional<Gaussian::ContractedRadialOrbitalSet> nullableContractedRadialOrbitalSet,
		                     std::vector<std::pair<Eigen::Index, Eigen::Index>> segmentationTable,
		                     std::optional<Gaussian::SemiLocalEcp> nullableSemiLocalEcp)
		    : Base(angularMomentum),
		      m_NullableContractedRadialOrbitalSet(std::move(nullableContractedRadialOrbitalSet)),
		      m_SegmentationTable(std::move(segmentationTable)), m_NullableSemiLocalEcp(std::move(nullableSemiLocalEcp))
		{
			/* NO CODE */
		}

		bool EqualsTo_Impl(const AngularMomentumBlock& other, const SecChem::Scalar tolerance) const noexcept
		{
			return m_NullableSemiLocalEcp.has_value() == other.m_NullableSemiLocalEcp.has_value()
			       && m_NullableContractedRadialOrbitalSet.has_value()
			                  == other.m_NullableContractedRadialOrbitalSet.has_value()
			       && (!m_NullableContractedRadialOrbitalSet.has_value()
			           || static_cast<const Base&>(*this).EqualsTo(other, tolerance))
			       && (!m_NullableSemiLocalEcp.has_value() || SemiLocalEcp().EqualsTo(other.SemiLocalEcp(), tolerance));
		}

		const Eigen::VectorXd& ExponentSet_Impl() const noexcept
		{
			return ContractedRadialOrbitalSet().ExponentSet();
		}

		const Eigen::MatrixXd& ContractionSets_Impl() const noexcept
		{
			return ContractedRadialOrbitalSet().ContractionSets();
		}

		template <typename InputIterator, typename Getter>
		static SegmentationTable ConcatSegmentationTable(const InputIterator begin, const InputIterator end, Getter get)
		{
			static_assert(std::is_lvalue_reference_v<decltype(get(*begin))>);
			const auto segmentCount = std::accumulate(begin,
			                                          end,
			                                          std::size_t{1},  // start from 1, not 0
			                                          [get](const std::size_t acc, const auto& segmentationTable)
			                                          {
				                                          if (get(segmentationTable).empty())
				                                          {
					                                          return acc;
				                                          }
				                                          return acc + get(segmentationTable).size() - 1;
			                                          });
			if (segmentCount == 1)
			{
				return {};
			}

			SegmentationTable result;
			result.reserve(segmentCount);

			std::pair<Eigen::Index, Eigen::Index> offsetPair = {0, 0};
			result.emplace_back(offsetPair);
			for (auto it = begin; it != end; ++it)
			{
				const SegmentationTable& segmentationTable = get(*it);

				if (segmentationTable.empty())
				{
					continue;
				}

				std::transform(std::next(segmentationTable.begin()),
				               segmentationTable.end(),
				               std::back_inserter(result),
				               [offsetPair](const std::pair<Eigen::Index, Eigen::Index>& seg)
				               { return std::pair{offsetPair.first + seg.first, offsetPair.second + seg.second}; });
				offsetPair = result.back();
			}

			return result;
		}

		std::optional<Gaussian::ContractedRadialOrbitalSet> m_NullableContractedRadialOrbitalSet = std::nullopt;
		SegmentationTable m_SegmentationTable = {};
		std::optional<Gaussian::SemiLocalEcp> m_NullableSemiLocalEcp = std::nullopt;
	};
}  // namespace SecChem::BasisSet::Gaussian
