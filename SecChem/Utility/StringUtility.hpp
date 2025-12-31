//
// Created by Andy on 12/25/2025.
//

#pragma once

#include <cstring>
#include <string>
#include <string_view>
#include <vector>

#define CONSTEXPR23

#if __has_include(<SecUtility/IO/StringUtility.hpp>)
#include <SecUtility/IO/StringUtility.hpp>
#else
namespace SecUtility
{
	inline constexpr auto ToLower = [](const unsigned char c) noexcept { return std::tolower(c); };
	inline constexpr auto ToUpper = [](const unsigned char c) noexcept { return std::toupper(c); };

	inline constexpr auto AsciiToLower = [](const unsigned char c) noexcept
	{ return (c >= 'A' && c <= 'Z') ? (c + ('a' - 'A')) : c; };
	inline constexpr auto AsciiToUpper = [](const unsigned char c) noexcept
	{ return (c >= 'a' && c <= 'z') ? (c + ('A' - 'a')) : c; };

	template <typename CharLowerer>
	constexpr int CaseInsensitiveComparison(const std::string_view lhs,
	                                        const std::string_view rhs,
	                                        const CharLowerer lowerCaseOf) noexcept
	{
		const std::size_t minSize = lhs.size() < rhs.size() ? lhs.size() : rhs.size();

		for (std::size_t i = 0; i < minSize; i++)
		{
			const auto la = lowerCaseOf(static_cast<unsigned char>(lhs[i]));
			const auto lb = lowerCaseOf(static_cast<unsigned char>(rhs[i]));

			if (la < lb)
			{
				return -1;
			}
			if (la > lb)
			{
				return 1;
			}
		}

		if (lhs.size() < rhs.size())
		{
			return -1;
		}
		if (lhs.size() > rhs.size())
		{
			return 1;
		}
		return 0;
	}

	template <typename StringLike, typename StringLikeToo>
	int CaseInsensitiveComparison(const StringLike& lhs, const StringLikeToo& rhs) noexcept
	{
		return CaseInsensitiveComparison(
		        static_cast<std::string_view>(lhs), static_cast<std::string_view>(rhs), ToLower);
	}

	template <typename StringLike, typename StringLikeToo>
	constexpr int CaseInsensitiveAsciiComparison(const StringLike& lhs, const StringLikeToo& rhs) noexcept
	{
		return CaseInsensitiveComparison(
		        static_cast<std::string_view>(lhs), static_cast<std::string_view>(rhs), AsciiToLower);
	}

	template <typename StringLike, typename StringLikeToo>
	bool IsCaseInsensitiveStringEqual(const StringLike& lhs, const StringLikeToo& rhs)
	{
		return CaseInsensitiveComparison(lhs, rhs) == 0;
	}

	template <typename StringLike, typename StringLikeToo>
	constexpr bool IsCaseInsensitiveAsciiStringEqual(const StringLike& lhs, const StringLikeToo& rhs)
	{
		return CaseInsensitiveAsciiComparison(lhs, rhs) == 0;
	}

	inline std::vector<std::string> SplitRespectingQuotes(const std::string& line)
	{
		std::vector<std::string> tokens;
		std::string current;
		bool isInQuotes = false;

		for (const char c : line)
		{
			if (c == '"')
			{
				isInQuotes = !isInQuotes;
				current.push_back(c);
			}
			else if (std::isspace(c) && !isInQuotes)
			{
				if (!current.empty())
				{
					tokens.push_back(current);
					current.clear();
				}
			}
			else
			{
				current.push_back(c);
			}
		}

		if (!current.empty())
		{
			tokens.push_back(current);
		}

		return tokens;
	}
}  // namespace SecUtility
#endif
