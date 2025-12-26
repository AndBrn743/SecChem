//
// Created by Andy on 11/23/2025.
//

#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <filesystem>
#include <memory>
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
			assert(m_ExponentSet.size() > 0);
			assert(m_ExponentSet.size() == m_ContractionSets.rows());
			assert(m_ExponentSet.size() >= m_ContractionSets.cols());
			assert(m_ContractionSets.cols() >= 0);
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

		Eigen::VectorX<Scalar> m_ExponentSet;
		Eigen::MatrixX<Scalar> m_ContractionSets;
		std::vector<ContractionSetViewDescription> m_ContractionSetViewDescriptions;

		bool EqualTo_Impl(const ContractedRadialOrbitalSet& other,
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
		SemiLocalEcp(const CoefficientSet& coefficientSet,
		             const RExponentSet& rExponentSet,
		             const GaussianExponentSet& gaussianExponentSet)
		    : m_Data(CreateDataSet(coefficientSet, rExponentSet, gaussianExponentSet))
		{
			static_assert(std::is_base_of_v<Eigen::EigenBase<CoefficientSet>, CoefficientSet>);
			static_assert(std::is_base_of_v<Eigen::EigenBase<RExponentSet>, RExponentSet>);
			static_assert(std::is_base_of_v<Eigen::EigenBase<GaussianExponentSet>, GaussianExponentSet>);
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


	private:
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

		bool EqualTo_Impl(const SemiLocalEcp& other, const Scalar tolerance = ZeroTolerance) const noexcept
		{
			return m_Data.rows() == other.m_Data.rows() && (m_Data - other.m_Data).cwiseAbs().maxCoeff() <= tolerance;
		}

		Eigen::Matrix<Scalar, Eigen::Dynamic, 3> m_Data;
	};

	class AngularMomentumBlock
	{
	public:
		AngularMomentumBlock(const AzimuthalQuantumNumber angularMomentum,
		                     ContractedRadialOrbitalSet contractedRadialOrbitalSet)
		    : m_AzimuthalQuantumNumber(angularMomentum),
		      m_ContractedRadialOrbitalSet(std::move(contractedRadialOrbitalSet))
		{
			/* NO CODE */
		}

		AngularMomentumBlock(const AzimuthalQuantumNumber angularMomentum,
		                     ContractedRadialOrbitalSet contractedRadialOrbitalSet,
		                     SemiLocalEcp ecp)
		    : m_AzimuthalQuantumNumber(angularMomentum),
		      m_ContractedRadialOrbitalSet(std::move(contractedRadialOrbitalSet)),
		      m_NullableSemiLocalEcp(std::move(ecp))
		{
			/* NO CODE */
		}

		AngularMomentumBlock& OverrideOrAddSemiLocalEcp(SemiLocalEcp ecp)
		{
			m_NullableSemiLocalEcp = std::move(ecp);
			return *this;
		}

		const AzimuthalQuantumNumber& AngularMomentum() const noexcept
		{
			return m_AzimuthalQuantumNumber;
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

		const Eigen::VectorXd& ExponentSet() const noexcept
		{
			return m_ContractedRadialOrbitalSet.ExponentSet();
		}

		const Eigen::MatrixXd& ContractionSets() const noexcept
		{
			return m_ContractedRadialOrbitalSet.ContractionSets();
		}

		Eigen::Index PrimitiveShellCount() const noexcept
		{
			return m_ContractedRadialOrbitalSet.PrimitiveShellCount();
		}

		Eigen::Index ContractedShellCount() const noexcept
		{
			return m_ContractedRadialOrbitalSet.ContractedShellCount();
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

	private:
		AzimuthalQuantumNumber m_AzimuthalQuantumNumber;
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