//
// Created by Andy on 11/23/2025.
//

#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <filesystem>
#include <memory>
#include <numeric>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

#include "../Utility/IEquatableWithTolerance.hpp"

#include "../AzimuthalQuantumNumber.hpp"
#include "../Element.hpp"
#include "../System.hpp"

namespace SecChem::BasisSet::Gaussian
{
	class ContractedRadialOrbitalSet;

	class SemiLocalEcp;

	template <typename>
	class AbstractAngularMomentumBlock;

	class AngularMomentumBlock;

	namespace Detail
	{
		template <OwnershipSemantics Semantics>
		class BasisSetImpl;

		template <OwnershipSemantics Semantics>
		class BasisSetLibraryImpl;
	}  // namespace Detail

	using BasisSet = Detail::BasisSetImpl<OwnershipSemantics::Value>;

	using SharedBasisSet = Detail::BasisSetImpl<OwnershipSemantics::Reference>;

	using BasisSetLibrary = Detail::BasisSetLibraryImpl<OwnershipSemantics::Value>;

	using SharedBasisSetLibrary = Detail::BasisSetLibraryImpl<OwnershipSemantics::Reference>;
}  // namespace SecChem::BasisSet::Gaussian

template <>
struct SecUtility::Traits<SecChem::BasisSet::Gaussian::ContractedRadialOrbitalSet>
{
	static constexpr auto DefaultEqualityComparisonTolerance = 1e-15;
};

template <>
struct SecUtility::Traits<SecChem::BasisSet::Gaussian::SemiLocalEcp>
{
	static constexpr auto DefaultEqualityComparisonTolerance = 1e-15;
};

template <typename Derived>
struct SecUtility::Traits<SecChem::BasisSet::Gaussian::AbstractAngularMomentumBlock<Derived>>
{
	static constexpr auto DefaultEqualityComparisonTolerance = 1e-15;
};

template <>
struct SecUtility::Traits<SecChem::BasisSet::Gaussian::AngularMomentumBlock>
{
	static constexpr auto DefaultEqualityComparisonTolerance = 1e-15;
};

namespace SecChem::BasisSet::Gaussian
{
	class ContractedRadialOrbitalSet : public SecUtility::IEquatableWithTolerance<ContractedRadialOrbitalSet>
	{
		friend IEquatableWithTolerance;

		struct ContractionSetViewDescription
		{
			Eigen::Index Offset;
			Eigen::Index Size;
		};

	public:
		static constexpr Scalar ZeroTolerance = 1e-15;

		class ContractionSetView
		{
		public:
			auto Exponents() const noexcept
			{
				return m_Exponents;
			}

			auto ContractionCoefficients() const noexcept
			{
				return m_Contractions;
			}

			Scalar Exponent(const Eigen::Index index) const noexcept
			{
				return m_Exponents[index];
			}

			Scalar ContractionCoefficient(const Eigen::Index index) const noexcept
			{
				return m_Contractions[index];
			}

		private:
			friend class ContractedRadialOrbitalSet;
			friend std::allocator<ContractionSetView>;
			ContractionSetView(const Eigen::Map<const Eigen::VectorX<Scalar>> exponents,     // NOLINT
			                   const Eigen::Map<const Eigen::VectorX<Scalar>> contractions)  // NOLINT
			    : m_Exponents(exponents), m_Contractions(contractions)
			{
				assert(m_Exponents.size() == m_Contractions.size());
			}

			const Eigen::Map<const Eigen::VectorX<Scalar>> m_Exponents;
			const Eigen::Map<const Eigen::VectorX<Scalar>> m_Contractions;
		};

		ContractedRadialOrbitalSet(Eigen::VectorXd exponentSet, Eigen::MatrixXd contractionSets)
		    : m_ExponentSet(std::move(exponentSet)), m_ContractionSets(std::move(contractionSets)),
		      m_ContractionSetViewDescriptions(BuildContractionSetViewDescriptions(m_ExponentSet, m_ContractionSets))
		{
			ValidateInput();
		}

		const Eigen::VectorXd& ExponentSet() const noexcept
		{
			return m_ExponentSet;
		}

		const Eigen::MatrixXd& ContractionSets() const noexcept
		{
			return m_ContractionSets;
		}

		Eigen::Index PrimitiveShellCount() const noexcept
		{
			return m_ExponentSet.size();
		}

		Eigen::Index ContractedShellCount() const noexcept
		{
			return m_ContractionSets.cols();
		}

		ContractionSetView ContractionSet(const Eigen::Index index) const noexcept
		{
			const auto [offset, size] = m_ContractionSetViewDescriptions[index];
			return {Eigen::Map<const Eigen::VectorX<Scalar>>{ExponentSet().data() + offset, size},
			        Eigen::Map<const Eigen::VectorX<Scalar>>{ContractionSets().col(index).data() + offset, size}};
		}

		auto AccumulationCoefficientSet(const Eigen::Index index) const noexcept
		{
			return m_ContractionSets.row(index);
		}

		template <typename InputIterator, typename Getter>
		static ContractedRadialOrbitalSet Concat(InputIterator begin, const InputIterator end, Getter get)
		{
			if (begin == end)
			{
				throw std::runtime_error("ContractedRadialOrbitalSet: input range is empty");
			}

			if (std::distance(begin, end) == 1)
			{
				return get(*begin);
			}

			using IndexPair = std::pair<Eigen::Index, Eigen::Index>;
			const auto dimensions =
			        std::accumulate(begin,
			                        end,
			                        IndexPair{0, 0},
			                        [get](const IndexPair& acc, const auto& block)
			                        {
				                        return IndexPair{acc.first + get(block).m_ContractionSets.rows(),
				                                         acc.second + get(block).m_ContractionSets.cols()};
			                        });

			Eigen::VectorX<Scalar> exponentSet = Eigen::VectorX<Scalar>::Zero(dimensions.first);
			Eigen::MatrixX<Scalar> contractionSets = Eigen::MatrixX<Scalar>::Zero(dimensions.first, dimensions.second);

			for (IndexPair offsets = {0, 0}; begin != end; ++begin)
			{
				exponentSet.segment(offsets.first, get(*begin).m_ExponentSet.size()) = get(*begin).m_ExponentSet;
				contractionSets.block(offsets.first,
				                      offsets.second,
				                      get(*begin).m_ContractionSets.rows(),
				                      get(*begin).m_ContractionSets.cols()) = get(*begin).m_ContractionSets;

				offsets.first += get(*begin).m_ContractionSets.rows();
				offsets.second += get(*begin).m_ContractionSets.cols();
			}

			return {std::move(exponentSet), std::move(contractionSets)};
		}

		template <typename InputIterator>
		static ContractedRadialOrbitalSet Concat(InputIterator begin, const InputIterator end)
		{
			return Concat(begin, end, [](auto&& item) -> decltype(auto) { return std::forward<decltype(item)>(item); });
		}

	private:
		static std::vector<ContractionSetViewDescription> BuildContractionSetViewDescriptions(
		        const Eigen::VectorXd& exponentSet, const Eigen::MatrixXd& contractionSets)
		{
			std::vector<ContractionSetViewDescription> descriptions;
			descriptions.reserve(exponentSet.size());

			for (Eigen::Index i = 0; i < contractionSets.cols(); i++)
			{
				const auto contractionSet = contractionSets.col(i);
				const auto head = std::find_if(contractionSet.begin(),
				                               contractionSet.end(),
				                               [](const Scalar c) { return std::abs(c) >= ZeroTolerance; });
				if (head == contractionSet.end())
				{
					throw std::runtime_error("Contraction coefficients from a contraction set is all zero");
				}

				const auto tail = std::find_if(std::reverse_iterator(contractionSet.end()),
				                               std::reverse_iterator{contractionSet.begin()},
				                               [](const Scalar c) { return std::abs(c) >= ZeroTolerance; });
				descriptions.push_back({std::distance(contractionSet.begin(), head), std::distance(head, tail.base())});
			}

			return descriptions;
		}

		// this method shall be called only from ctors
		void ValidateInput() const
		{
			if (m_ExponentSet.size() == 0)
			{
				throw std::invalid_argument("ContractedRadialOrbitalSet: empty exponent set");
			}

			if (m_ContractionSets.size() == 0)
			{
				throw std::invalid_argument("ContractedRadialOrbitalSet: empty coefficient set");
			}

			if (m_ExponentSet.size() != m_ContractionSets.rows())
			{
				throw std::invalid_argument("ContractedRadialOrbitalSet: size mismatch");
			}
		}

		Eigen::VectorX<Scalar> m_ExponentSet;
		Eigen::MatrixX<Scalar> m_ContractionSets;
		std::vector<ContractionSetViewDescription> m_ContractionSetViewDescriptions;

		bool EqualsTo_Impl(const ContractedRadialOrbitalSet& other,
		                   const Scalar tolerance = ZeroTolerance) const noexcept
		{
			return m_ExponentSet.size() == other.m_ExponentSet.size()
			       && m_ContractionSets.rows() == other.m_ContractionSets.rows()
			       && m_ContractionSets.cols() == other.m_ContractionSets.cols()
			       && (m_ExponentSet - other.m_ExponentSet).cwiseAbs().maxCoeff() <= tolerance
			       && (m_ContractionSets - other.m_ContractionSets).cwiseAbs().maxCoeff() <= tolerance;
		}
	};

	class SemiLocalEcp : public SecUtility::IEquatableWithTolerance<SemiLocalEcp>
	{
		friend IEquatableWithTolerance;

	public:
		static constexpr Scalar ZeroTolerance = 1e-15;

		template <typename CoefficientSet, typename RExponentSet, typename GaussianExponentSet>
		SemiLocalEcp(const Eigen::Index electronCount,
		             const CoefficientSet& coefficientSet,
		             const RExponentSet& rExponentSet,
		             const GaussianExponentSet& gaussianExponentSet)
		    : m_EcpElectronCount(electronCount),
		      m_Data(CreateDataSet(coefficientSet, rExponentSet, gaussianExponentSet))
		{
			if (electronCount < 0)
			{
				throw std::runtime_error("ECP electron count can't be negative");
			}

			static_assert(std::is_base_of_v<Eigen::EigenBase<CoefficientSet>, CoefficientSet>);
			static_assert(std::is_base_of_v<Eigen::EigenBase<RExponentSet>, RExponentSet>);
			static_assert(std::is_base_of_v<Eigen::EigenBase<GaussianExponentSet>, GaussianExponentSet>);
		}

		Eigen::Index ElectronCount() const noexcept
		{
			return m_EcpElectronCount;
		}

		Eigen::Index TermCount() const noexcept
		{
			return m_Data.rows();
		}

		auto Coefficients() const noexcept
		{
			return m_Data.col(0);
		}

		auto RExponents() const noexcept
		{
			return m_Data.col(1);
		}

		auto GaussianExponents() const noexcept
		{
			return m_Data.col(2);
		}

		const auto& Coefficient(const Eigen::Index index) const noexcept
		{
			return m_Data(index, 0);
		}

		const auto& RExponent(const Eigen::Index index) const noexcept
		{
			return m_Data(index, 1);
		}

		const auto& GaussianExponent(const Eigen::Index index) const noexcept
		{
			return m_Data(index, 2);
		}

		template <typename InputIterator, typename Getter>
		static SemiLocalEcp Concat(InputIterator begin, const InputIterator end, const Getter get)
		{
			if (begin == end)
			{
				throw std::runtime_error("SemiLocalEcp: input range is empty");
			}

			if (std::distance(begin, end) == 1)
			{
				return get(*begin);
			}

			const auto totalEcpElectronCount = std::accumulate(begin,
			                                                   end,
			                                                   Eigen::Index{0},
			                                                   [get](const Eigen::Index acc, const auto& ecp)
			                                                   { return acc + get(ecp).m_EcpElectronCount; });

			Eigen::Matrix<Scalar, Eigen::Dynamic, 3> data(std::accumulate(begin,
			                                                              end,
			                                                              Eigen::Index{0},
			                                                              [get](const Eigen::Index acc, const auto& ecp)
			                                                              { return acc + get(ecp).m_Data.rows(); }),
			                                              3);

			for (Eigen::Index offset = 0; begin != end; ++begin)
			{
				data.middleRows(offset, get(*begin).m_Data.rows()) = get(*begin).m_Data;
				offset += get(*begin).m_Data.rows();
			}

			return SemiLocalEcp{totalEcpElectronCount, std::move(data)};
		}

		template <typename InputIterator>
		static SemiLocalEcp Concat(InputIterator begin, const InputIterator end)
		{
			return Concat(begin, end, [](auto&& item) -> decltype(auto) { return std::forward<decltype(item)>(item); });
		}

		template <typename InputIterator, typename Getter>
		static std::optional<SemiLocalEcp> ConcatNullable(InputIterator begin,
		                                                  const InputIterator end,
		                                                  const Getter get)
		{
			if (begin == end)
			{
				return std::nullopt;
			}

			if (std::distance(begin, end) == 1)
			{
				return get(*begin);
			}

			const auto totalEcpElectronCount = std::accumulate(begin,
			                                                   end,
			                                                   Eigen::Index{0},
			                                                   [get](const Eigen::Index acc, const auto& ecp)
			                                                   {
				                                                   if (!get(ecp).has_value())
				                                                   {
					                                                   return acc;
				                                                   }
				                                                   return acc + get(ecp).value().m_EcpElectronCount;
			                                                   });

			const auto totalTermCount = std::accumulate(begin,
			                                            end,
			                                            Eigen::Index{0},
			                                            [get](const Eigen::Index acc, const auto& ecp)
			                                            {
				                                            if (!get(ecp).has_value())
				                                            {
					                                            return acc;
				                                            }
				                                            return acc + get(ecp).value().m_Data.rows();
			                                            });

			if (totalTermCount == 0)
			{
				return std::nullopt;
			}

			Eigen::Matrix<Scalar, Eigen::Dynamic, 3> data(totalTermCount, 3);

			for (Eigen::Index offset = 0; begin != end; ++begin)
			{
				if (get(*begin).has_value())
				{
					const auto& ecpData = get(*begin).value().m_Data;
					data.middleRows(offset, ecpData.rows()) = ecpData;
					offset += ecpData.rows();
				}
			}

			return SemiLocalEcp{totalEcpElectronCount, std::move(data)};
		}

		template <typename InputIterator>
		static std::optional<SemiLocalEcp> ConcatNullable(InputIterator begin, const InputIterator end)
		{
			return ConcatNullable(
			        begin, end, [](auto&& item) -> decltype(auto) { return std::forward<decltype(item)>(item); });
		}

	private:
		explicit SemiLocalEcp(const Eigen::Index electronCount, Eigen::Matrix<Scalar, Eigen::Dynamic, 3> data)
		    : m_EcpElectronCount(electronCount), m_Data(std::move(data))
		{
			/* NO CODE */
		}

		template <typename CoefficientSet, typename RExponentSet, typename GaussianExponentSet>
		static Eigen::Matrix<Scalar, Eigen::Dynamic, 3> CreateDataSet(const CoefficientSet& coefficientSet,
		                                                              const RExponentSet& rExponentSet,
		                                                              const GaussianExponentSet& gaussianExponentSet)
		{
			ValidateInput(coefficientSet, rExponentSet, gaussianExponentSet);
			Eigen::Matrix<Scalar, Eigen::Dynamic, 3> data(coefficientSet.size(), 3);
			data.col(0) = coefficientSet;
			data.col(1) = rExponentSet;
			data.col(2) = gaussianExponentSet;
			return data;
		}

		template <typename C, typename R, typename G>
		static void ValidateInput(const C& c, const R& r, const G& g)
		{
			if (c.size() == 0)
			{
				throw std::invalid_argument("SemiLocalEcp: empty coefficient set");
			}

			if (c.size() != r.size() || c.size() != g.size())
			{
				throw std::invalid_argument("SemiLocalEcp: size mismatch");
			}

			if (!(c.rows() == 1 || c.cols() == 1))
			{
				throw std::invalid_argument("SemiLocalEcp: coefficients must be a vector");
			}

			if (!(r.rows() == 1 || r.cols() == 1))
			{
				throw std::invalid_argument("SemiLocalEcp: r exponents must be a vector");
			}

			if (!(g.rows() == 1 || g.cols() == 1))
			{
				throw std::invalid_argument("SemiLocalEcp: gaussian exponents must be a vector");
			}
		}

		bool EqualsTo_Impl(const SemiLocalEcp& other, const Scalar tolerance = ZeroTolerance) const noexcept
		{
			return m_EcpElectronCount == other.m_EcpElectronCount && m_Data.rows() == other.m_Data.rows()
			       && (m_Data - other.m_Data).cwiseAbs().maxCoeff() <= tolerance;
		}

		Eigen::Index m_EcpElectronCount;
		Eigen::Matrix<Scalar, Eigen::Dynamic, 3> m_Data;
	};

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

		bool EqualsTo(const AbstractAngularMomentumBlock& other, const Scalar tolerance) const noexcept
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

	class AngularMomentumBlock : public SecUtility::IEquatableWithTolerance<AngularMomentumBlock>,
	                             public AbstractAngularMomentumBlock<AngularMomentumBlock>
	{
		friend IEquatableWithTolerance<AngularMomentumBlock>;
		using Base = AbstractAngularMomentumBlock;
		friend Base;

	public:
		AngularMomentumBlock(const AzimuthalQuantumNumber angularMomentum,
		                     ContractedRadialOrbitalSet contractedRadialOrbitalSet)
		    : Base(angularMomentum), m_ContractedRadialOrbitalSet(std::move(contractedRadialOrbitalSet))
		{
			/* NO CODE */
		}

		AngularMomentumBlock(const AzimuthalQuantumNumber angularMomentum,
		                     ContractedRadialOrbitalSet contractedRadialOrbitalSet,
		                     SemiLocalEcp ecp)
		    : Base(angularMomentum), m_ContractedRadialOrbitalSet(std::move(contractedRadialOrbitalSet)),
		      m_NullableSemiLocalEcp(std::move(ecp))
		{
			/* NO CODE */
		}

		AngularMomentumBlock& AllOrOverrideSemiLocalEcp(SemiLocalEcp ecp)
		{
			m_NullableSemiLocalEcp = std::move(ecp);
			return *this;
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
			return m_ContractedRadialOrbitalSet;
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

			auto crs = ContractedRadialOrbitalSet::Concat(
			        begin, end, [get](const auto& block) { return get(block).m_ContractedRadialOrbitalSet; });
			auto ecp = SemiLocalEcp::ConcatNullable(
			        begin, end, [get](const auto& block) { return get(block).m_NullableSemiLocalEcp; });

			return ecp.has_value() ? AngularMomentumBlock{azimuthalQuantumNumber, crs, ecp.value()}
			                       : AngularMomentumBlock{azimuthalQuantumNumber, crs};
		}

		template <typename InputIterator>
		static AngularMomentumBlock Concat(InputIterator begin, const InputIterator end)
		{
			return Concat(begin, end, [](auto&& item) -> decltype(auto) { return std::forward<decltype(item)>(item); });
		}

		using IEquatableWithTolerance<AngularMomentumBlock>::EqualsTo;
		using IEquatableWithTolerance<AngularMomentumBlock>::NotEqualsTo;
		using IEquatableWithTolerance<AngularMomentumBlock>::operator==;
		using IEquatableWithTolerance<AngularMomentumBlock>::operator!=;

	private:
		bool EqualsTo_Impl(const AngularMomentumBlock& other, const Scalar tolerance) const noexcept
		{
			return m_NullableSemiLocalEcp.has_value() == other.m_NullableSemiLocalEcp.has_value()
			       && static_cast<const Base&>(*this).EqualsTo(static_cast<const Base&>(other), tolerance)
			       && (!m_NullableSemiLocalEcp.has_value() || SemiLocalEcp().EqualsTo(other.SemiLocalEcp(), tolerance));
		}

		const Eigen::VectorXd& ExponentSet_Impl() const noexcept
		{
			return m_ContractedRadialOrbitalSet.ExponentSet();
		}

		const Eigen::MatrixXd& ContractionSets_Impl() const noexcept
		{
			return m_ContractedRadialOrbitalSet.ContractionSets();
		}

		Gaussian::ContractedRadialOrbitalSet m_ContractedRadialOrbitalSet;
		std::optional<Gaussian::SemiLocalEcp> m_NullableSemiLocalEcp;
	};

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
					return lhs.AngularMomentum() <= rhs.AngularMomentum();
				}

				if (lhs.ExponentSet().size() != rhs.ExponentSet().size())
				{
					return lhs.ExponentSet().size() <= rhs.ExponentSet().size();
				}

				return lhs.ExponentSet()[0] <= rhs.ExponentSet()[0];
			};

		public:
			BasisSetImpl() : m_DataStorage(CreateStorage())
			{
				/* NO CODE */
			}

			explicit BasisSetImpl(const BasisSetImpl<AlternativeOf(Semantics)>& other)
			{
				if constexpr (Semantics == OwnershipSemantics::Value)
				{
					m_DataStorage = *other.m_DataStorage;
				}
				else
				{
					m_DataStorage = std::make_shared<DataType>(other.m_DataStorage);
				}
			}

			explicit BasisSetImpl(BasisSetImpl<AlternativeOf(Semantics)>&& other)
			{
				static_assert(AlternativeOf(Semantics) == OwnershipSemantics::Value);
				m_DataStorage = std::make_shared<DataType>(std::move(other.m_DataStorage));
			}

			BasisSetImpl Clone() const
			{
				static_assert(Semantics == OwnershipSemantics::Reference);
				return {std::make_shared<DataType>(*m_DataStorage)};
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
					        "OverwriteEntryOf(...), or rarely, AddOrOverwriteEntryOf(...) instead.");
				}
				return Data()[element];
			}

			std::vector<AngularMomentumBlock>& OverwriteEntryOf(const Element element)
			{
				Data().at(element) = {};
				return Data().at(element);
			}

			bool Has(const Element element) const
			{
				return Data().find(element) != Data().end();
			}

			/// Value comparison will be done for BasisSet, reference comparison will be done for SharedBasisSet.
			/// To compare the value of SharedBasisSet, one should use <c>EqualsTo</c> method instead.
			/// Please note that value comparison is depended on the underlying storage representation.
			template <OwnershipSemantics OtherSemantics>
			bool operator==(const BasisSetImpl<OtherSemantics>& other) const noexcept
			{
				return Semantics == OtherSemantics && m_DataStorage == other.m_DataStorage;
			}

			/// Value comparison will be done for BasisSet, reference comparison will be done for SharedBasisSet.
			/// To compare the value of SharedBasisSet, one should use <c>NotEqualsTo</c> method instead
			/// Please note that value comparison is depended on the underlying storage representation.
			template <OwnershipSemantics OtherSemantics>
			bool operator!=(const BasisSetImpl<OtherSemantics>& other) const noexcept
			{
				return !(*this == other);
			}

			template <OwnershipSemantics OtherSemantics>
			bool EqualsTo(const BasisSetImpl<OtherSemantics>& other, const Scalar tolerance = 1e-15) const noexcept
			{
				if (Data().size() != other.Size())
				{
					return false;
				}

				for (const auto& [element, angularMomentumBlocks] : Data())
				{
					const auto iteratorToOtherAngularMomentumBlocks = other.find(element);
					if (iteratorToOtherAngularMomentumBlocks == other.end())
					{
						return false;
					}

					const auto& otherAngularMomentumBlocks = *iteratorToOtherAngularMomentumBlocks;
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
					if (!std::is_sorted(
					            angularMomentumBlocks.cbegin(), angularMomentumBlocks.cend(), IsInStandardStorageOrder))
					{
						return false;
					}
				}

				return true;
			}

			void StandardizeRepresentation()
			{
				for (auto& [_, angularMomentumBlocks] : Data())
				{
					std::sort(angularMomentumBlocks.begin(), angularMomentumBlocks.end(), IsInStandardStorageOrder);
				}
			}

			BasisSetImpl ToStandardizeRepresentation() const
			{
				BasisSetImpl result;
				result.StandardizeRepresentation();
				return result;
			}

			bool IsInStandardAndConcatenatedRepresentation() const noexcept
			{
				if (Data().empty())
				{
					return true;
				}

				for (const auto& [_, angularMomentumBlocks] : Data())
				{
					if (!std::is_sorted(
					            angularMomentumBlocks.cbegin(), angularMomentumBlocks.cend(), IsInStandardStorageOrder)
					    || !AngularMomentumBlockConcatSets(angularMomentumBlocks).empty())
					{
						return false;
					}
				}

				return true;
			}

			void StandardizeAndConcatenatedRepresentation()
			{
				for (auto& [_, angularMomentumBlocks] : Data())
				{
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
			}

			BasisSetImpl ToStandardizeAndConcatenatedRepresentation() const
			{
				BasisSetImpl result;
				result.StandardizeAndConcatenatedRepresentation();
				return result;
			}


		private:
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
					                            { return block.AngularMomentum() != l0; });

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
