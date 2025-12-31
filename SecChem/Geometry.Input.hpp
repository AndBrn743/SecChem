//
// Created by Andy on 12/31/2025.
//

#pragma once

#include "Element.hpp"
#include "Utility/Parser.hpp"
#include <Eigen/Dense>

namespace SecChem::Geometry::Input
{
	struct BasicGeometryLineParsingResult
	{
		// GCC require it
		// ReSharper disable once CppRedundantElaboratedTypeSpecifier
		class Element Element;
		Eigen::Vector3d Position;
	};


	template <typename TokenIterator, typename LengthToBohrRadius>
	auto ParseBasicCartesianCoordinateLine(const TokenIterator tokenBegin, const LengthToBohrRadius toBohrRadius)
	        -> std::enable_if_t<std::is_same_v<decltype(SecUtility::Parse<double>(tokenBegin[0])), double>
	                                    && std::is_same_v<std::decay_t<decltype(toBohrRadius(3.14))>, double>,
	                            BasicGeometryLineParsingResult>
	{
		return {Element::Symbol2Element(tokenBegin[0]),
		        {toBohrRadius(SecUtility::Parse<double>(tokenBegin[1])),
		         toBohrRadius(SecUtility::Parse<double>(tokenBegin[2])),
		         toBohrRadius(SecUtility::Parse<double>(tokenBegin[3]))}};
	}


	template <typename RangeOfTokens, typename LengthToBohrRadius>
	auto ParseBasicCartesianCoordinateLine(const RangeOfTokens tokens, const LengthToBohrRadius toBohrRadius)
	        -> std::enable_if_t<std::is_same_v<decltype(SecUtility::Parse<double>(tokens.begin()[0])), double>
	                                    && std::is_same_v<decltype(SecUtility::Parse<double>(tokens.end()[0])), double>
	                                    && std::is_same_v<std::decay_t<decltype(toBohrRadius(3.14))>, double>,
	                            BasicGeometryLineParsingResult>
	{
		const auto tokenCount = std::distance(tokens.begin(), tokens.end());
		if (tokenCount != 4)
		{
			throw std::runtime_error("Unexpected number of tokens. Expecting 4, got " + std::to_string(tokenCount));
		}

		return ParseBasicCartesianCoordinateLine(tokens.begin(), toBohrRadius);
	}


	template <std::size_t IndexingConvention = 0,
	          typename TokenIterator,
	          typename ReferenceAtomIterator,
	          typename LengthToBohrRadius,
	          typename AngleToRadian>
	auto ParseInternalCoordinateLine(TokenIterator tokenBegin,
	                                 const TokenIterator tokenEnd,
	                                 const ReferenceAtomIterator atomBegin,
	                                 const ReferenceAtomIterator atomEnd,
	                                 const LengthToBohrRadius toBohrRadius,
	                                 const AngleToRadian toRadian)
	        -> std::enable_if_t<std::is_same_v<decltype(SecUtility::Parse<std::int64_t>(*tokenBegin)), std::int64_t>
	                                    && std::is_same_v<decltype(SecUtility::Parse<double>(*tokenBegin)), double>
	                                    && std::is_same_v<std::decay_t<decltype(*atomBegin)>, Atom>
	                                    && std::is_same_v<std::decay_t<decltype(atomBegin[0])>, Atom>
	                                    && std::is_same_v<std::decay_t<decltype(toBohrRadius(3.14))>, double>
	                                    && std::is_same_v<std::decay_t<decltype(toRadian(3.14))>, double>,
	                            BasicGeometryLineParsingResult>
	{
		using SecUtility::Parse;
		static_assert(IndexingConvention == 0 || IndexingConvention == 1, "IndexingConvention must be set to 0 or 1");
		static_assert(std::is_same_v<std::decay_t<decltype(toRadian(3.14))>, double>);

		const auto atomCount = std::distance(atomBegin, atomEnd);
		const auto tokenCount = std::distance(tokenBegin, tokenEnd);
		if (tokenCount == 0)
		{
			throw std::runtime_error("No token was provided [LOC43]");
		}

		const auto element = Element::Symbol2Element(*tokenBegin);
		++tokenBegin;

		if (atomCount == 0)
		{
			if (tokenCount != 1)
			{
				throw std::runtime_error("Unexpected number of tokens. Expecting 1, got " + std::to_string(tokenCount)
				                         + " instead");
			}

			return {element, {}};
		}

		//--------------------------------------------------

		if (tokenCount < 3)
		{
			throw std::runtime_error("Unexpected number of tokens. Expecting at least 3, got "
			                         + std::to_string(tokenCount) + " instead");
		}

		const auto i = Parse<std::uint64_t>(*tokenBegin) - IndexingConvention;
		++tokenBegin;
		const auto length = toBohrRadius(Parse<double>(*tokenBegin));
		++tokenBegin;

		if (length <= 0)
		{
			throw std::runtime_error("Bound length must be positive [LOC44]");
		}
		if (i >= atomCount)
		{
			throw std::out_of_range("Bad atom index [LOC44]");
		}

		if (atomCount == 1)
		{
			if (tokenCount != 3)
			{
				throw std::runtime_error("Unexpected number of tokens. Expecting 3, got " + std::to_string(tokenCount)
				                         + " instead");
			}

			return {element, atomBegin[0].Position() + Eigen::Vector3d{0, 0, length}};
		}

		//--------------------------------------------------

		if (tokenCount < 5)
		{
			throw std::runtime_error("Unexpected number of tokens. Expecting at least 5, got "
			                         + std::to_string(tokenCount) + " instead");
		}

		const auto j = Parse<std::uint64_t>(*tokenBegin) - IndexingConvention;
		++tokenBegin;
		const auto angle = toRadian(Parse<double>(*tokenBegin));
		++tokenBegin;

		if (j >= atomCount)
		{
			throw std::out_of_range("Bad atom index [LOC45]");
		}
		if (i == j)
		{
			throw std::runtime_error("Duplication of reference atom [LOC45]");
		}

		if (atomCount == 2)
		{
			if (tokenCount != 5)
			{
				throw std::runtime_error("Unexpected number of tokens. Expecting 5, got " + std::to_string(tokenCount)
				                         + " instead");
			}

			if (i == 0 /* && j == 1, which is implicit*/)
			{
				return {element,
				        atomBegin[i].Position() + length * Eigen::Vector3d{std::sin(angle), 0, std::cos(angle)}};
			}
			else  // i == 1 && j == 0
			{
				return {element,
				        atomBegin[i].Position() + length * Eigen::Vector3d{std::sin(angle), 0, -std::cos(angle)}};
			}
		}

		//--------------------------------------------------

		if (tokenCount != 7)
		{
			throw std::runtime_error("Unexpected number of tokens. Expecting 7, got " + std::to_string(tokenCount)
			                         + " instead");
		}

		const auto k = Parse<std::uint64_t>(*tokenBegin) - IndexingConvention;
		++tokenBegin;
		const auto dihedral = toRadian(Parse<double>(*tokenBegin));
		++tokenBegin;

		if (k >= atomCount)
		{
			throw std::out_of_range("Bad atom index [LOC49]");
		}
		if (i == k || j == k)
		{
			throw std::runtime_error("Duplication of reference atom [LOC49]");
		}

		const Eigen::Vector3d e1 = (atomBegin[j].Position() - atomBegin[i].Position()).normalized();

		const Eigen::Vector3d v2 = (atomBegin[k].Position() - atomBegin[j].Position()).cross(e1);
		const auto norm = v2.norm();
		if (norm < 1e-9)
		{
			throw std::runtime_error("Z-matrix reference atoms are collinear [LOC50]");
		}

		const Eigen::Vector3d e2 = v2 / norm;
		const Eigen::Vector3d e3 = e1.cross(e2);

		return {element,
		        atomBegin[i].Position()
		                + length
		                          * (std::cos(angle) * e1
		                             + std::sin(angle) * (std::cos(dihedral) * e3 + std::sin(dihedral) * e2))};
	}


	template <std::size_t IndexingConvention = 0,
	          typename RangeOfTokens,
	          typename RangeOfReferenceAtoms,
	          typename LengthToBohrRadius,
	          typename AngleToRadian>
	auto ParseInternalCoordinateLine(const RangeOfTokens& tokens,
	                                 const RangeOfReferenceAtoms& atoms,
	                                 const LengthToBohrRadius toBohrRadius,
	                                 const AngleToRadian toRadian)
	        -> std::enable_if_t<
	                std::is_same_v<decltype(SecUtility::Parse<std::int64_t>(*tokens.begin())), std::int64_t>
	                        && std::is_same_v<decltype(SecUtility::Parse<std::int64_t>(*tokens.end())), std::int64_t>
	                        && std::is_same_v<decltype(SecUtility::Parse<double>(*tokens.begin())), double>
	                        && std::is_same_v<decltype(SecUtility::Parse<double>(*tokens.end())), double>
	                        && std::is_same_v<std::decay_t<decltype(*atoms.begin())>, Atom>
	                        && std::is_same_v<std::decay_t<decltype(*atoms.end())>, Atom>
	                        && std::is_same_v<std::decay_t<decltype(atoms.begin()[0])>, Atom>
	                        && std::is_same_v<std::decay_t<decltype(toBohrRadius(3.14))>, double>
	                        && std::is_same_v<std::decay_t<decltype(toRadian(3.14))>, double>,
	                BasicGeometryLineParsingResult>
	{
		return ParseInternalCoordinateLine<IndexingConvention>(
		        tokens.begin(), tokens.end(), atoms.begin(), atoms.end(), toBohrRadius, toRadian);
	}


	template <std::size_t IndexingConvention = 0,
	          typename TokenIterator,
	          typename RangeOfReferenceAtoms,
	          typename LengthToBohrRadius,
	          typename AngleToRadian>
	auto ParseInternalCoordinateLine(TokenIterator tokenBegin,
	                                 const TokenIterator tokenEnd,
	                                 const RangeOfReferenceAtoms& atoms,
	                                 const LengthToBohrRadius toBohrRadius,
	                                 const AngleToRadian toRadian)
	        -> std::enable_if_t<std::is_same_v<decltype(SecUtility::Parse<std::int64_t>(*tokenBegin)), std::int64_t>
	                                    && std::is_same_v<decltype(SecUtility::Parse<double>(*tokenBegin)), double>
	                                    && std::is_same_v<std::decay_t<decltype(*atoms.begin())>, Atom>
	                                    && std::is_same_v<std::decay_t<decltype(*atoms.end())>, Atom>
	                                    && std::is_same_v<std::decay_t<decltype(atoms.begin()[0])>, Atom>
	                                    && std::is_same_v<std::decay_t<decltype(toBohrRadius(3.14))>, double>
	                                    && std::is_same_v<std::decay_t<decltype(toRadian(3.14))>, double>,
	                            BasicGeometryLineParsingResult>
	{
		return ParseInternalCoordinateLine<IndexingConvention>(
		        tokenBegin, tokenEnd, atoms.begin(), atoms.end(), toBohrRadius, toRadian);
	}


	template <std::size_t IndexingConvention = 0,
	          typename RangeOfTokens,
	          typename ReferenceAtomIterator,
	          typename LengthToBohrRadius,
	          typename AngleToRadian>
	auto ParseInternalCoordinateLine(const RangeOfTokens& tokens,
	                                 const ReferenceAtomIterator atomBegin,
	                                 const ReferenceAtomIterator atomEnd,
	                                 const LengthToBohrRadius toBohrRadius,
	                                 const AngleToRadian toRadian)
	        -> std::enable_if_t<
	                std::is_same_v<decltype(SecUtility::Parse<std::int64_t>(*tokens.begin())), std::int64_t>
	                        && std::is_same_v<decltype(SecUtility::Parse<std::int64_t>(*tokens.end())), std::int64_t>
	                        && std::is_same_v<decltype(SecUtility::Parse<double>(*tokens.begin())), double>
	                        && std::is_same_v<decltype(SecUtility::Parse<double>(*tokens.end())), double>
	                        && std::is_same_v<std::decay_t<decltype(*atomBegin)>, Atom>
	                        && std::is_same_v<std::decay_t<decltype(atomBegin[0])>, Atom>
	                        && std::is_same_v<std::decay_t<decltype(toBohrRadius(3.14))>, double>
	                        && std::is_same_v<std::decay_t<decltype(toRadian(3.14))>, double>,
	                BasicGeometryLineParsingResult>
	{
		return ParseInternalCoordinateLine<IndexingConvention>(
		        tokens.begin(), tokens.end(), atomBegin, atomEnd, toBohrRadius, toRadian);
	}
}  // namespace SecChem::Geometry::Input
