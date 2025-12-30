#pragma once
#if !defined(SECCHEM_GAUSSIAN_BASIS_SET_INTERNAL)
#error Do not include internal header files directly
#endif

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
}  // namespace SecChem::BasisSet::Gaussian
