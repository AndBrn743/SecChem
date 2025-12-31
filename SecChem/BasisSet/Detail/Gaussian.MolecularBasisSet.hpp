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
		using BasisAssignment = std::pair<std::size_t, const std::vector<BasisSet::Gaussian::AngularMomentumBlock>*>;

	public:
		const BasisSet::Gaussian::SharedBasisSetLibrary& SharedBasisSetLibrary() const noexcept
		{
			return m_Library;
		}

		const SharedMolecule& Molecule() const noexcept
		{
			return m_Molecule;
		}


	private:
		MolecularBasisSet(const BasisSet::Gaussian::SharedBasisSetLibrary& library,
		                  const SharedMolecule& molecule,
		                  std::vector<BasisAssignment> basisAssignments)
		    : m_Library(library), m_Molecule(molecule), m_BasisAssignments(std::move(basisAssignments))
		{
			assert(molecule.AtomCount() == m_BasisAssignments.size());
		}

		BasisSet::Gaussian::SharedBasisSetLibrary m_Library;
		SharedMolecule m_Molecule;
		std::vector<BasisAssignment> m_BasisAssignments;
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

	Builder& SetGlobalDefaultBasisSetTo(const std::string& globalDefaultBasisSetName)
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

	Builder& SetElementaryDefaultBasisSetTo(const std::string& defaultBasisSetName, const Element element)
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
		std::vector<Gaussian::MolecularBasisSet::BasisAssignment> assignments{};

		for (std::string line; std::getline(input, line); /* NO CODE */)
		{
			const auto [atom, basisPtr] = parser.ParseLine(
			        line, std::as_const(atoms), std::as_const(m_Library), std::as_const(m_DefaultBasisSet));
			assert(basisPtr != nullptr);
			assignments.emplace_back(atoms.size(), basisPtr);
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
