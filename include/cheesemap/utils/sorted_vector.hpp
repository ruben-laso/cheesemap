#pragma once

#include <algorithm>
#include <utility>
#include <vector>

namespace chs
{
	template<typename T, typename Compare = std::less<T>, typename Allocator = std::allocator<T>>
	class sorted_vector : public std::vector<T, Allocator>
	{
		using base_type = std::vector<T, Allocator>;

		using size_type = typename std::vector<T, Allocator>::size_type;

		size_type max_size_ = std::numeric_limits<size_type>::max();

		Compare comp_{};

		public:
		using base_type::base_type;

		sorted_vector(const std::size_t max_size) : max_size_(max_size) { this->reserve(max_size_); }

		auto full() const { return std::cmp_equal(this->size(), max_size_); }

		auto insert(const T & value) -> typename base_type::iterator
		{
			const auto is_full = not std::cmp_less(this->size(), max_size_);

			// Do not insert the element if the vector is full and cmp(value, back()) is false
			if (is_full and not comp_(value, this->back())) { return this->end(); }

			auto it = std::lower_bound(this->begin(), this->end(), value, comp_);

			// If the element is found and the vector is not full
			if (not is_full) { it = this->emplace(it, value); }
			// If the vector is full, displace elements to the right
			else if (it != this->end())
			{
				// Shift elements to the right (last element is lost)
				std::move_backward(it, this->end() - 1, this->end());
				*it = value;
			}

			return it;
		}

		auto insert(T && value) -> typename base_type::iterator
		{
			const auto is_full = not std::cmp_less(this->size(), max_size_);

			// Do not insert the element if the vector is full and cmp(value, back())
			if (is_full and not comp_(value, this->back())) { return this->end(); }

			const auto it = std::lower_bound(this->begin(), this->end(), value, comp_);

			// If the element is found and the vector is not full
			if (not is_full) { this->emplace(it, std::move(value)); }
			// If the vector is full, displace elements to the right
			else if (it != this->end())
			{
				// Shift elements to the right (last element is lost)
				std::move_backward(it, this->end() - 1, this->end());
				*it = std::move(value);
			}

			return it;
		}
	};
} // namespace chs