#include "micromechanics/cell_placement.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <compare>
#include <numbers>
#include <stdexcept>

namespace physicore::mechanics::micromechanics {

namespace {
// A lattice site: its squared distance from the origin in lattice units, which is an integer and so
// compares exactly, and its integer lattice coordinates.
struct lattice_site
{
	sindex_t squared_distance;
	std::array<sindex_t, 3> coordinates;

	// by distance first; ties are broken by the coordinates, so the order is fully determined
	auto operator<=>(const lattice_site&) const = default;
};

// Sites of a line with unit spacing, within reach of the origin
std::vector<lattice_site> line_sites(sindex_t reach)
{
	std::vector<lattice_site> sites;
	for (sindex_t i = -reach; i <= reach; ++i)
	{
		sites.push_back({ i * i, { i, 0, 0 } });
	}
	return sites;
}

// Sites i * (1, 0) + j * (1/2, sqrt(3)/2) of the hexagonal lattice; squared distance i^2 + i j + j^2
std::vector<lattice_site> hexagonal_sites(sindex_t reach)
{
	std::vector<lattice_site> sites;
	for (sindex_t j = -2 * reach; j <= 2 * reach; ++j)
	{
		for (sindex_t i = -2 * reach; i <= 2 * reach; ++i)
		{
			const sindex_t squared_distance = i * i + i * j + j * j;
			if (squared_distance <= reach * reach)
			{
				sites.push_back({ squared_distance, { i, j, 0 } });
			}
		}
	}
	return sites;
}

// Sites of the face-centred cubic lattice: integer triples with an even sum, in units of spacing / sqrt(2);
// the squared distance in those units is i^2 + j^2 + k^2, and nearest neighbours are at squared distance 2
std::vector<lattice_site> face_centred_cubic_sites(sindex_t reach)
{
	std::vector<lattice_site> sites;
	for (sindex_t k = -reach; k <= reach; ++k)
	{
		for (sindex_t j = -reach; j <= reach; ++j)
		{
			for (sindex_t i = -reach; i <= reach; ++i)
			{
				const sindex_t squared_distance = i * i + j * j + k * k;
				if ((i + j + k) % 2 == 0 && squared_distance <= reach * reach)
				{
					sites.push_back({ squared_distance, { i, j, k } });
				}
			}
		}
	}
	return sites;
}

// Radius, in units of the spacing, of the ball that holds count sites at the density of the lattice
real_t ball_radius(index_t dims, index_t count)
{
	const auto n = static_cast<real_t>(count);
	if (dims == 1)
	{
		return 0.5 * n;
	}
	if (dims == 2)
	{
		// 2 / sqrt(3) sites per unit area
		return std::sqrt(n * std::numbers::sqrt3 / (2.0 * std::numbers::pi));
	}
	// sqrt(2) sites per unit volume
	return std::cbrt(3.0 * n / (4.0 * std::numbers::pi * std::numbers::sqrt2));
}
} // namespace

std::vector<real_t> packed_ball_positions(index_t dims, index_t count, real_t spacing)
{
	if (dims < 1 || dims > 3)
	{
		throw std::invalid_argument("Elements can be packed in 1, 2 or 3 dimensions");
	}
	if (spacing <= 0.0)
	{
		throw std::invalid_argument("The spacing of packed elements must be positive");
	}

	// Collect all sites within a ball that is certainly large enough, then keep the closest ones.
	// The margin of two spacings covers the difference between a ball and the lattice sites near it.
	const real_t radius = ball_radius(dims, count) + 2.0;

	std::vector<lattice_site> sites;
	if (dims == 1)
	{
		sites = line_sites(static_cast<sindex_t>(std::ceil(radius)));
	}
	else if (dims == 2)
	{
		sites = hexagonal_sites(static_cast<sindex_t>(std::ceil(radius)));
	}
	else
	{
		// lattice units are spacing / sqrt(2)
		sites = face_centred_cubic_sites(static_cast<sindex_t>(std::ceil(radius * std::numbers::sqrt2)));
	}

	if (sites.size() < count)
	{
		throw std::logic_error("Not enough lattice sites were generated");
	}
	std::ranges::sort(sites);

	std::vector<real_t> positions(count * dims);
	for (index_t n = 0; n < count; ++n)
	{
		const auto i = static_cast<real_t>(sites[n].coordinates[0]);
		const auto j = static_cast<real_t>(sites[n].coordinates[1]);
		const auto k = static_cast<real_t>(sites[n].coordinates[2]);

		if (dims == 1)
		{
			positions[n] = spacing * i;
		}
		else if (dims == 2)
		{
			positions[2 * n] = spacing * (i + 0.5 * j);
			positions[2 * n + 1] = spacing * 0.5 * std::numbers::sqrt3 * j;
		}
		else
		{
			const real_t unit = spacing / std::numbers::sqrt2;
			positions[3 * n] = unit * i;
			positions[3 * n + 1] = unit * j;
			positions[3 * n + 2] = unit * k;
		}
	}

	return positions;
}

} // namespace physicore::mechanics::micromechanics
