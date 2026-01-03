//
// Created by Andy on 1/3/2026.
//

#pragma once

#include <algorithm>
#include <functional>
#include <iterator>
#include <stdexcept>
#include <utility>
#include <vector>


namespace polyfill
{
	template <typename Key,
	          typename T,
	          typename Compare = std::less<Key>,
	          typename KeyContainer = std::vector<Key>,
	          typename MappedContainer = std::vector<T>>
	class flat_map
	{
	public:
		using key_type = Key;
		using mapped_type = T;
		using value_type = std::pair<Key, T>;
		using size_type = std::size_t;
		using difference_type = std::ptrdiff_t;
		using key_compare = Compare;
		using reference = std::pair<const Key&, T&>;
		using const_reference = std::pair<const Key&, const T&>;

		// Iterator implementation
		template <bool IsConst>
		class iterator_impl
		{
		public:
			using iterator_category = std::random_access_iterator_tag;
			using value_type = std::conditional_t<IsConst, const flat_map::value_type, flat_map::value_type>;
			using difference_type = flat_map::difference_type;
			using pointer = value_type*;
			using reference = std::conditional_t<IsConst, flat_map::const_reference, flat_map::reference>;

		private:
			using key_iter =
			        std::conditional_t<IsConst, typename KeyContainer::const_iterator, typename KeyContainer::iterator>;
			using mapped_iter = std::conditional_t<IsConst,
			                                       typename MappedContainer::const_iterator,
			                                       typename MappedContainer::iterator>;

			key_iter k_it;
			mapped_iter m_it;

		public:
			iterator_impl() = default;
			iterator_impl(key_iter k, mapped_iter m) : k_it(k), m_it(m)
			{
			}

			// Allow conversion from non-const to const iterator
			template <bool WasConst, typename = std::enable_if_t<IsConst && !WasConst>>
			iterator_impl(const iterator_impl<WasConst>& other) : k_it(other.k_it), m_it(other.m_it)
			{
			}

			reference operator*() const
			{
				return reference{*k_it, *m_it};
			}

			iterator_impl& operator++()
			{
				++k_it;
				++m_it;
				return *this;
			}
			iterator_impl operator++(int)
			{
				auto tmp = *this;
				++(*this);
				return tmp;
			}
			iterator_impl& operator--()
			{
				--k_it;
				--m_it;
				return *this;
			}
			iterator_impl operator--(int)
			{
				auto tmp = *this;
				--(*this);
				return tmp;
			}

			iterator_impl& operator+=(difference_type n)
			{
				k_it += n;
				m_it += n;
				return *this;
			}
			iterator_impl& operator-=(difference_type n)
			{
				k_it -= n;
				m_it -= n;
				return *this;
			}

			iterator_impl operator+(difference_type n) const
			{
				return iterator_impl(k_it + n, m_it + n);
			}
			iterator_impl operator-(difference_type n) const
			{
				return iterator_impl(k_it - n, m_it - n);
			}

			difference_type operator-(const iterator_impl& other) const
			{
				return k_it - other.k_it;
			}

			bool operator==(const iterator_impl& other) const
			{
				return k_it == other.k_it;
			}
			bool operator!=(const iterator_impl& other) const
			{
				return k_it != other.k_it;
			}
			bool operator<(const iterator_impl& other) const
			{
				return k_it < other.k_it;
			}
			bool operator>(const iterator_impl& other) const
			{
				return k_it > other.k_it;
			}
			bool operator<=(const iterator_impl& other) const
			{
				return k_it <= other.k_it;
			}
			bool operator>=(const iterator_impl& other) const
			{
				return k_it >= other.k_it;
			}

			friend class flat_map;
			template <bool>
			friend class iterator_impl;
		};

		using iterator = iterator_impl<false>;
		using const_iterator = iterator_impl<true>;
		using reverse_iterator = std::reverse_iterator<iterator>;
		using const_reverse_iterator = std::reverse_iterator<const_iterator>;

	private:
		KeyContainer keys;
		MappedContainer values;
		[[no_unique_address]] Compare comp;

		template <typename K>
		auto lower_bound_impl(const K& k) const
		{
			return std::lower_bound(keys.begin(), keys.end(), k, comp);
		}

		template <typename K>
		auto lower_bound_impl(const K& k)
		{
			return std::lower_bound(keys.begin(), keys.end(), k, comp);
		}

	public:
		// Constructors
		flat_map() = default;
		explicit flat_map(const Compare& c) : comp(c)
		{
		}

		template <typename InputIt>
		flat_map(InputIt first, InputIt last, const Compare& c = Compare()) : comp(c)
		{
			insert(first, last);
		}

		flat_map(std::initializer_list<value_type> init, const Compare& c = Compare()) : comp(c)
		{
			insert(init.begin(), init.end());
		}

		// Capacity
		bool empty() const noexcept
		{
			return keys.empty();
		}
		size_type size() const noexcept
		{
			return keys.size();
		}
		size_type max_size() const noexcept
		{
			return keys.max_size();
		}

		// Iterators
		iterator begin() noexcept
		{
			return iterator(keys.begin(), values.begin());
		}
		const_iterator begin() const noexcept
		{
			return const_iterator(keys.begin(), values.begin());
		}
		const_iterator cbegin() const noexcept
		{
			return begin();
		}

		iterator end() noexcept
		{
			return iterator(keys.end(), values.end());
		}
		const_iterator end() const noexcept
		{
			return const_iterator(keys.end(), values.end());
		}
		const_iterator cend() const noexcept
		{
			return end();
		}

		reverse_iterator rbegin() noexcept
		{
			return reverse_iterator(end());
		}
		const_reverse_iterator rbegin() const noexcept
		{
			return const_reverse_iterator(end());
		}
		const_reverse_iterator crbegin() const noexcept
		{
			return rbegin();
		}

		reverse_iterator rend() noexcept
		{
			return reverse_iterator(begin());
		}
		const_reverse_iterator rend() const noexcept
		{
			return const_reverse_iterator(begin());
		}
		const_reverse_iterator crend() const noexcept
		{
			return rend();
		}

		// Element access
		T& at(const Key& k)
		{
			auto it = lower_bound_impl(k);
			if (it == keys.end() || comp(k, *it))
			{
				throw std::out_of_range("flat_map::at");
			}
			auto idx = std::distance(keys.begin(), it);
			return values[idx];
		}

		const T& at(const Key& k) const
		{
			auto it = lower_bound_impl(k);
			if (it == keys.end() || comp(k, *it))
			{
				throw std::out_of_range("flat_map::at");
			}
			auto idx = std::distance(keys.begin(), it);
			return values[idx];
		}

		T& operator[](const Key& k)
		{
			auto it = lower_bound_impl(k);
			if (it == keys.end() || comp(k, *it))
			{
				auto idx = std::distance(keys.begin(), it);
				keys.insert(it, k);
				return *values.insert(values.begin() + idx, T{});
			}
			auto idx = std::distance(keys.begin(), it);
			return values[idx];
		}

		T& operator[](Key&& k)
		{
			auto it = lower_bound_impl(k);
			if (it == keys.end() || comp(k, *it))
			{
				auto idx = std::distance(keys.begin(), it);
				keys.insert(it, std::move(k));
				return *values.insert(values.begin() + idx, T{});
			}
			auto idx = std::distance(keys.begin(), it);
			return values[idx];
		}

		// Modifiers
		void clear() noexcept
		{
			keys.clear();
			values.clear();
		}

		std::pair<iterator, bool> insert(const value_type& v)
		{
			return emplace(v.first, v.second);
		}

		std::pair<iterator, bool> insert(value_type&& v)
		{
			return emplace(std::move(v.first), std::move(v.second));
		}

		template <typename InputIt>
		void insert(InputIt first, InputIt last)
		{
			for (; first != last; ++first)
			{
				insert(*first);
			}
		}

		template <typename... Args>
		std::pair<iterator, bool> emplace(Args&&... args)
		{
			value_type v(std::forward<Args>(args)...);
			auto it = lower_bound_impl(v.first);

			if (it != keys.end() && !comp(v.first, *it))
			{
				auto idx = std::distance(keys.begin(), it);
				return {iterator(keys.begin() + idx, values.begin() + idx), false};
			}

			auto idx = std::distance(keys.begin(), it);
			keys.insert(it, std::move(v.first));
			values.insert(values.begin() + idx, std::move(v.second));
			return {iterator(keys.begin() + idx, values.begin() + idx), true};
		}

		iterator erase(const_iterator pos)
		{
			auto idx = pos.k_it - keys.begin();
			keys.erase(keys.begin() + idx);
			values.erase(values.begin() + idx);
			return iterator(keys.begin() + idx, values.begin() + idx);
		}

		iterator erase(const_iterator first, const_iterator last)
		{
			auto first_idx = first.k_it - keys.begin();
			auto last_idx = last.k_it - keys.begin();
			keys.erase(keys.begin() + first_idx, keys.begin() + last_idx);
			values.erase(values.begin() + first_idx, values.begin() + last_idx);
			return iterator(keys.begin() + first_idx, values.begin() + first_idx);
		}

		size_type erase(const Key& k)
		{
			auto it = lower_bound_impl(k);
			if (it == keys.end() || comp(k, *it))
			{
				return 0;
			}
			auto idx = std::distance(keys.begin(), it);
			keys.erase(it);
			values.erase(values.begin() + idx);
			return 1;
		}

		// Lookup
		iterator find(const Key& k)
		{
			auto it = lower_bound_impl(k);
			if (it == keys.end() || comp(k, *it))
			{
				return end();
			}
			auto idx = std::distance(keys.begin(), it);
			return iterator(keys.begin() + idx, values.begin() + idx);
		}

		const_iterator find(const Key& k) const
		{
			auto it = lower_bound_impl(k);
			if (it == keys.end() || comp(k, *it))
			{
				return end();
			}
			auto idx = std::distance(keys.begin(), it);
			return const_iterator(keys.begin() + idx, values.begin() + idx);
		}

		size_type count(const Key& k) const
		{
			return find(k) != end() ? 1 : 0;
		}

		bool contains(const Key& k) const
		{
			return find(k) != end();
		}

		iterator lower_bound(const Key& k)
		{
			auto it = lower_bound_impl(k);
			auto idx = std::distance(keys.begin(), it);
			return iterator(keys.begin() + idx, values.begin() + idx);
		}

		const_iterator lower_bound(const Key& k) const
		{
			auto it = lower_bound_impl(k);
			auto idx = std::distance(keys.begin(), it);
			return const_iterator(keys.begin() + idx, values.begin() + idx);
		}

		iterator upper_bound(const Key& k)
		{
			auto it = std::upper_bound(keys.begin(), keys.end(), k, comp);
			auto idx = std::distance(keys.begin(), it);
			return iterator(keys.begin() + idx, values.begin() + idx);
		}

		const_iterator upper_bound(const Key& k) const
		{
			auto it = std::upper_bound(keys.begin(), keys.end(), k, comp);
			auto idx = std::distance(keys.begin(), it);
			return const_iterator(keys.begin() + idx, values.begin() + idx);
		}

		std::pair<iterator, iterator> equal_range(const Key& k)
		{
			return {lower_bound(k), upper_bound(k)};
		}

		std::pair<const_iterator, const_iterator> equal_range(const Key& k) const
		{
			return {lower_bound(k), upper_bound(k)};
		}

		// Observers
		key_compare key_comp() const
		{
			return comp;
		}
	};
}  // namespace polyfill