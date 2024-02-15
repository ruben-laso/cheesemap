#pragma once

#include <array>
#include <optional>
#include <vector>

#include "cheesemap/utils/Box.hpp"
#include "cheesemap/utils/Cell.hpp"

namespace chs::Slice
{
	template<typename Point_type>
	class Smart
	{
		protected:
		static constexpr auto SPARSE_TO_DENSE_THRESHOLD = 0.8;

		static constexpr auto Dim = 2;

		using resolution_type   = double;
		using dimensions_array  = std::array<resolution_type, Dim>;
		using dimensions_vector = std::vector<resolution_type>;
		using indices_array     = std::array<std::size_t, Dim>;
		using indices_vector    = std::vector<std::size_t>;
		using cell_type         = Cell<Point_type>;

		static constexpr dimensions_array DEFAULT_RESOLUTIONS = []() {
			std::array<resolution_type, Dim> res{};
			res.fill(1);
			return res;
		}();

		// Dimension of each cell
		dimensions_vector resolutions_{ DEFAULT_RESOLUTIONS.begin(), DEFAULT_RESOLUTIONS.end() };

		// Bounding box of the map
		Box box_;

		// Number of cells of the map on each dimension
		indices_array sizes_{ Dim };

		// Cells of the map (using sparse representation)
		std::vector<cell_type>   cells_sparse_;
		arma::SpMat<std::size_t> indices_;

		// Cells of the map (using dense representation)
		std::vector<cell_type> cells_dense_;

		// Use sparse
		bool use_sparse_ = true;

		[[nodiscard]] inline auto idx2box(const std::size_t i, const std::size_t j) const
		{
			Point center{ box_.min()[0] + (static_cast<resolution_type>(i) + 0.5) * resolutions_[0],
				      box_.min()[1] + (static_cast<resolution_type>(j) + 0.5) * resolutions_[1], 0 };
			Point radii{ resolutions_[0] / 2, resolutions_[1] / 2, 0 };

			return Box{ center, radii };
		}

		inline void add_point_dense(Point_type & point)
		{
			const auto [i, j] = coord2indices(point);
			const auto idx    = i * sizes_[1] + j;
			cells_dense_[idx].add_point(&point);
		}

		inline void add_point_sparse(Point_type & point)
		{
			const auto [i, j] = coord2indices(point);

			// If the cell is empty
			if (std::as_const(indices_).at(i, j) == 0)
			{
				// Create a new cell
				cells_sparse_.emplace_back(Cell<Point_type>{ idx2box(i, j) });
				// Store the index of the cell in the sparse matrix
				indices_.at(i, j) = cells_sparse_.size();
			}

			// Insert in (i, j) the index of the cell - 1
			// (-1 because we can't store a 0 in a sparse matrix)
			const auto idx = std::as_const(indices_).at(i, j) - 1;
			cells_sparse_[idx].add_point(&point);

			if (density() > SPARSE_TO_DENSE_THRESHOLD) { sparse2dense(); }
		}

		inline void sparse2dense()
		{
			cells_dense_ = {};

			for (const auto i : ranges::views::indices(sizes_[0]))
			{
				for (const auto j : ranges::views::indices(sizes_[1]))
				{
					const auto idx = std::as_const(indices_).at(i, j);
					if (idx == 0) { cells_dense_.emplace_back(Cell<Point_type>{ idx2box(i, j) }); }
					else { cells_dense_.emplace_back(std::move(cells_sparse_[idx - 1])); }
				}
			}

			use_sparse_ = false;

			indices_.resize(0, 0);
			cells_sparse_ = {};
		}

		inline void dense2sparse()
		{
			cells_sparse_ = {};
			indices_.reset();

			for (const auto i : ranges::views::indices(sizes_[0]))
			{
				for (const auto j : ranges::views::indices(sizes_[1]))
				{
					const auto idx = i * sizes_[1] + j;
					if (not cells_dense_[idx].empty())
					{
						cells_sparse_.emplace_back(cells_dense_[idx]);
						indices_.at(i, j) = cells_sparse_.size();
					}
				}
			}

			use_sparse_ = true;

			cells_dense_ = {};
		}

		public:
		Smart() = delete;

		Smart(const Box & box, const resolution_type res) : Smart(box, dimensions_vector(Dim, res)) {}

		Smart(const Box & box, const dimensions_vector & res) :
		        resolutions_(res),
		        box_(box),
		        sizes_({ static_cast<std::size_t>(std::ceil((box.max()[0] - box.min()[0]) / res[0])),
		                 static_cast<std::size_t>(std::ceil((box.max()[1] - box.min()[1]) / res[1])) }),
		        indices_(sizes_[0], sizes_[1])
		{}

		[[nodiscard]] inline auto resolutions() const -> const auto & { return resolutions_; }

		[[nodiscard]] inline auto resolutions() { return resolutions_; }

		[[nodiscard]] inline auto size() const { return sizes_[0] * sizes_[1]; }

		[[nodiscard]] inline auto density() const
		{
			if (use_sparse_)
			{
				return static_cast<double>(indices_.n_nonzero) / static_cast<double>(indices_.n_elem);
			}
			return 1.0;
		}

		[[nodiscard]] inline auto coord2indices(const Point_type & p) const
		{
			const auto rel  = p - box_.min();
			const auto size = box_.max() - box_.min();
			const auto i    = static_cast<std::size_t>(std::clamp(rel[0], 0.0, size[0]) / resolutions_[0]);
			const auto j    = static_cast<std::size_t>(std::clamp(rel[1], 0.0, size[1]) / resolutions_[1]);
			return std::make_pair(i, j);
		}

		[[nodiscard]] inline auto to_global_idx_as_dense(const std::size_t i, const std::size_t j) const
		{
			return i * sizes_[1] + j;
		}

		inline void add_point(Point_type & point)
		{
			if (use_sparse_) { add_point_sparse(point); }
			else { add_point_dense(point); }
		}

		[[nodiscard]] inline auto at(const std::size_t i, const std::size_t j)
		        -> std::optional<std::reference_wrapper<cell_type>>
		{
			if (use_sparse_)
			{
				const auto idx = std::as_const(indices_).at(i, j);
				if (std::cmp_equal(idx, 0)) { return std::nullopt; }
				return { cells_sparse_[idx - 1] };
			}

			return { cells_dense_[i * sizes_[1] + j] };
		}

		[[nodiscard]] inline auto at(const std::size_t i, const std::size_t j) const
		        -> std::optional<std::reference_wrapper<const cell_type>>
		{
			if (use_sparse_)
			{
				const auto idx = std::as_const(indices_).at(i, j);
				if (std::cmp_equal(idx, 0)) { return std::nullopt; }
				return { cells_sparse_[idx - 1] };
			}

			return { cells_dense_[i * sizes_[1] + j] };
		}
	};
} // namespace chs::Slice