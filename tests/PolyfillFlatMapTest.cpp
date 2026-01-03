//
// Created by Andy on 1/3/2026.
//

#include <SecChem/Polyfill/flat_map.hpp>
#include <catch2/catch_test_macros.hpp>
using polyfill::flat_map;

TEST_CASE("flat_map - Construction", "[flat_map][construction]")
{
	SECTION("Default construction")
	{
		flat_map<int, std::string> m;
		REQUIRE(m.empty());
		REQUIRE(m.size() == 0);
	}

	SECTION("Initializer list construction")
	{
		flat_map<int, std::string> m = {{1, "one"}, {3, "three"}, {2, "two"}};
		REQUIRE(m.size() == 3);
		REQUIRE(m[1] == "one");
		REQUIRE(m[2] == "two");
		REQUIRE(m[3] == "three");
	}

	SECTION("Range construction")
	{
		std::vector<std::pair<int, std::string>> vec = {{5, "five"}, {1, "one"}, {3, "three"}};
		flat_map<int, std::string> m(vec.begin(), vec.end());
		REQUIRE(m.size() == 3);
		REQUIRE(m[1] == "one");
		REQUIRE(m[3] == "three");
		REQUIRE(m[5] == "five");
	}

	SECTION("Custom comparator")
	{
		flat_map<int, std::string, std::greater<>> m = {{1, "one"}, {2, "two"}, {3, "three"}};
		auto it = m.begin();
		REQUIRE((*it).first == 3);
		++it;
		REQUIRE((*it).first == 2);
		++it;
		REQUIRE((*it).first == 1);
	}
}

TEST_CASE("flat_map - Insertion", "[flat_map][insertion]")
{
	flat_map<int, std::string> m;

	SECTION("Insert lvalue")
	{
		std::pair<int, std::string> p = {1, "one"};
		auto [it, inserted] = m.insert(p);
		REQUIRE(inserted);
		REQUIRE((*it).first == 1);
		REQUIRE((*it).second == "one");
		REQUIRE(m.size() == 1);
	}

	SECTION("Insert rvalue")
	{
		auto [it, inserted] = m.insert({2, "two"});
		REQUIRE(inserted);
		REQUIRE((*it).first == 2);
		REQUIRE((*it).second == "two");
	}

	SECTION("Insert duplicate key")
	{
		m.insert({1, "one"});
		auto [it, inserted] = m.insert({1, "ONE"});
		REQUIRE_FALSE(inserted);
		REQUIRE((*it).second == "one");  // Original value unchanged
		REQUIRE(m.size() == 1);
	}

	SECTION("Emplace")
	{
		auto [it, inserted] = m.emplace(3, "three");
		REQUIRE(inserted);
		REQUIRE((*it).first == 3);
		REQUIRE((*it).second == "three");
	}

	SECTION("Emplace duplicate")
	{
		m.emplace(1, "one");
		auto [it, inserted] = m.emplace(1, "ONE");
		REQUIRE_FALSE(inserted);
		REQUIRE((*it).second == "one");
	}

	SECTION("Range insert")
	{
		std::vector<std::pair<int, std::string>> vec = {{1, "one"}, {2, "two"}, {3, "three"}};
		m.insert(vec.begin(), vec.end());
		REQUIRE(m.size() == 3);
		REQUIRE(m[1] == "one");
		REQUIRE(m[2] == "two");
		REQUIRE(m[3] == "three");
	}

	SECTION("Multiple insertions maintain sorted order")
	{
		m.insert({5, "five"});
		m.insert({1, "one"});
		m.insert({3, "three"});
		m.insert({2, "two"});
		m.insert({4, "four"});

		auto it = m.begin();
		for (int i = 1; i <= 5; ++i)
		{
			REQUIRE((*it).first == i);
			++it;
		}
	}
}

TEST_CASE("flat_map - Element access", "[flat_map][access]")
{
	flat_map<int, std::string> m = {{1, "one"}, {2, "two"}, {3, "three"}};

	SECTION("operator[] existing key")
	{
		REQUIRE(m[1] == "one");
		REQUIRE(m[2] == "two");
		REQUIRE(m[3] == "three");
	}

	SECTION("operator[] non-existing key")
	{
		m[4] = "four";
		REQUIRE(m.size() == 4);
		REQUIRE(m[4] == "four");
	}

	SECTION("operator[] creates default value")
	{
		auto& val = m[5];
		REQUIRE(val == "");
		REQUIRE(m.size() == 4);
	}

	SECTION("at() existing key")
	{
		REQUIRE(m.at(1) == "one");
		REQUIRE(m.at(2) == "two");
	}

	SECTION("at() non-existing key throws")
	{
		REQUIRE_THROWS_AS(m.at(99), std::out_of_range);
	}

	SECTION("at() const version")
	{
		const auto& cm = m;
		REQUIRE(cm.at(1) == "one");
		REQUIRE_THROWS_AS(cm.at(99), std::out_of_range);
	}

	SECTION("Modify through operator[]")
	{
		m[1] = "ONE";
		REQUIRE(m[1] == "ONE");
		REQUIRE(m.size() == 3);  // Size unchanged
	}
}

TEST_CASE("flat_map - Erasure", "[flat_map][erasure]")
{
	SECTION("Erase by key")
	{
		flat_map<int, std::string> m = {{1, "one"}, {2, "two"}, {3, "three"}};
		auto count = m.erase(2);
		REQUIRE(count == 1);
		REQUIRE(m.size() == 2);
		REQUIRE(m.find(2) == m.end());
		REQUIRE(m[1] == "one");
		REQUIRE(m[3] == "three");
	}

	SECTION("Erase non-existing key")
	{
		flat_map<int, std::string> m = {{1, "one"}};
		auto count = m.erase(99);
		REQUIRE(count == 0);
		REQUIRE(m.size() == 1);
	}

	SECTION("Erase by iterator")
	{
		flat_map<int, std::string> m = {{1, "one"}, {2, "two"}, {3, "three"}};
		auto it = m.find(2);
		auto next = m.erase(it);
		REQUIRE(m.size() == 2);
		REQUIRE((*next).first == 3);
		REQUIRE(m.find(2) == m.end());
	}

	SECTION("Erase range")
	{
		flat_map<int, std::string> m = {{1, "one"}, {2, "two"}, {3, "three"}, {4, "four"}, {5, "five"}};
		auto first = m.find(2);
		auto last = m.find(4);
		auto next = m.erase(first, last);
		REQUIRE(m.size() == 3);
		REQUIRE((*next).first == 4);
		REQUIRE(m.find(2) == m.end());
		REQUIRE(m.find(3) == m.end());
	}

	SECTION("Clear")
	{
		flat_map<int, std::string> m = {{1, "one"}, {2, "two"}};
		m.clear();
		REQUIRE(m.empty());
		REQUIRE(m.size() == 0);
	}
}

TEST_CASE("flat_map - Lookup", "[flat_map][lookup]")
{
	flat_map<int, std::string> m = {{1, "one"}, {3, "three"}, {5, "five"}};

	SECTION("find existing key")
	{
		auto it = m.find(3);
		REQUIRE(it != m.end());
		REQUIRE((*it).first == 3);
		REQUIRE((*it).second == "three");
	}

	SECTION("find non-existing key")
	{
		auto it = m.find(2);
		REQUIRE(it == m.end());
	}

	SECTION("count")
	{
		REQUIRE(m.count(1) == 1);
		REQUIRE(m.count(3) == 1);
		REQUIRE(m.count(2) == 0);
		REQUIRE(m.count(99) == 0);
	}

	SECTION("contains")
	{
		REQUIRE(m.contains(1));
		REQUIRE(m.contains(3));
		REQUIRE(m.contains(5));
		REQUIRE_FALSE(m.contains(2));
		REQUIRE_FALSE(m.contains(4));
	}

	SECTION("lower_bound")
	{
		auto it = m.lower_bound(3);
		REQUIRE((*it).first == 3);

		it = m.lower_bound(2);
		REQUIRE((*it).first == 3);

		it = m.lower_bound(0);
		REQUIRE((*it).first == 1);

		it = m.lower_bound(6);
		REQUIRE(it == m.end());
	}

	SECTION("upper_bound")
	{
		auto it = m.upper_bound(3);
		REQUIRE((*it).first == 5);

		it = m.upper_bound(2);
		REQUIRE((*it).first == 3);

		it = m.upper_bound(5);
		REQUIRE(it == m.end());
	}

	SECTION("equal_range")
	{
		auto [first, last] = m.equal_range(3);
		REQUIRE((*first).first == 3);
		REQUIRE((*last).first == 5);
		REQUIRE(std::distance(first, last) == 1);

		auto [first2, last2] = m.equal_range(2);
		REQUIRE((*first2).first == 3);
		REQUIRE(first2 == last2);
	}
}

TEST_CASE("flat_map - Iterators", "[flat_map][iterators]")
{
	flat_map<int, std::string> m = {{1, "one"}, {2, "two"}, {3, "three"}};

	SECTION("Forward iteration")
	{
		auto it = m.begin();
		REQUIRE((*it).first == 1);
		++it;
		REQUIRE((*it).first == 2);
		++it;
		REQUIRE((*it).first == 3);
		++it;
		REQUIRE(it == m.end());
	}

	SECTION("Backward iteration")
	{
		auto it = m.end();
		--it;
		REQUIRE((*it).first == 3);
		--it;
		REQUIRE((*it).first == 2);
		--it;
		REQUIRE((*it).first == 1);
		REQUIRE(it == m.begin());
	}

	SECTION("Reverse iteration")
	{
		auto it = m.rbegin();
		REQUIRE((*it).first == 3);
		++it;
		REQUIRE((*it).first == 2);
		++it;
		REQUIRE((*it).first == 1);
		++it;
		REQUIRE(it == m.rend());
	}

	SECTION("Iterator arithmetic")
	{
		auto it = m.begin();
		it += 2;
		REQUIRE((*it).first == 3);

		it -= 1;
		REQUIRE((*it).first == 2);

		auto it2 = it + 1;
		REQUIRE((*it2).first == 3);

		auto it3 = it2 - 1;
		REQUIRE((*it3).first == 2);
	}

	SECTION("Iterator difference")
	{
		auto first = m.begin();
		auto last = m.end();
		REQUIRE(last - first == 3);
	}

	SECTION("Range-based for loop")
	{
		int expected = 1;
		for (auto [key, value] : m)
		{
			REQUIRE(key == expected);
			++expected;
		}
	}

	SECTION("Const iterators")
	{
		const auto& cm = m;
		auto it = cm.begin();
		REQUIRE((*it).first == 1);

		auto cit = m.cbegin();
		REQUIRE((*cit).first == 1);
	}
}

TEST_CASE("flat_map - Capacity", "[flat_map][capacity]")
{
	SECTION("Empty map")
	{
		flat_map<int, std::string> m;
		REQUIRE(m.empty());
		REQUIRE(m.size() == 0);
	}

	SECTION("Non-empty map")
	{
		flat_map<int, std::string> m = {{1, "one"}, {2, "two"}};
		REQUIRE_FALSE(m.empty());
		REQUIRE(m.size() == 2);
	}

	SECTION("Size changes with operations")
	{
		flat_map<int, std::string> m;
		m[1] = "one";
		REQUIRE(m.size() == 1);
		m[2] = "two";
		REQUIRE(m.size() == 2);
		m.erase(1);
		REQUIRE(m.size() == 1);
		m.clear();
		REQUIRE(m.size() == 0);
	}
}

TEST_CASE("flat_map - Complex types", "[flat_map][complex]")
{
	SECTION("String keys")
	{
		flat_map<std::string, int> m = {{"one", 1}, {"two", 2}, {"three", 3}};
		REQUIRE(m["one"] == 1);
		REQUIRE(m["three"] == 3);
		REQUIRE(m.size() == 3);
	}

	SECTION("String values with move semantics")
	{
		flat_map<int, std::string> m;
		std::string s = "temporary";
		m.insert({1, std::move(s)});
		REQUIRE(m[1] == "temporary");
	}

	SECTION("Non-copyable values")
	{
		struct MoveOnly
		{
			int value;
			MoveOnly(int v) : value(v)
			{
			}
			MoveOnly(const MoveOnly&) = delete;
			MoveOnly& operator=(const MoveOnly&) = delete;
			MoveOnly(MoveOnly&&) = default;
			MoveOnly& operator=(MoveOnly&&) = default;
		};

		flat_map<int, MoveOnly> m;
		m.emplace(1, MoveOnly(42));
		REQUIRE(m.at(1).value == 42);
	}
}

TEST_CASE("flat_map - Edge cases", "[flat_map][edge]")
{
	SECTION("Single element")
	{
		flat_map<int, std::string> m;
		m[1] = "one";
		REQUIRE(m.size() == 1);
		REQUIRE(m.begin() != m.end());
		REQUIRE(m[1] == "one");
	}

	SECTION("Duplicate insertions")
	{
		flat_map<int, std::string> m;
		m.insert({1, "one"});
		m.insert({1, "ONE"});
		m.insert({1, "One"});
		REQUIRE(m.size() == 1);
		REQUIRE(m[1] == "one");
	}

	SECTION("Erase all elements one by one")
	{
		flat_map<int, std::string> m = {{1, "one"}, {2, "two"}, {3, "three"}};
		m.erase(1);
		m.erase(2);
		m.erase(3);
		REQUIRE(m.empty());
	}

	SECTION("Insert after clear")
	{
		flat_map<int, std::string> m = {{1, "one"}, {2, "two"}};
		m.clear();
		m[3] = "three";
		REQUIRE(m.size() == 1);
		REQUIRE(m[3] == "three");
	}
}
