// Prints the intra- and inter-cellular pair potential and radial force as CSV.
//
// The notebooks in mechanics/micromechanics/docs/notebooks compare their own formulas against this
// output, so the documentation cannot silently drift away from the implementation.
//
// Usage: mechanics.micromechanics.potential_curve [equilibrium_distance steepness intra_well_depth
//                                                  inter_well_depth cutoff_factor [samples]]

#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <span>
#include <string>

#include <micromechanics/cell_type.h>
#include <micromechanics/pair_potential.h>
#include <micromechanics/simulation_parameters.h>

using namespace physicore;
using namespace physicore::mechanics::micromechanics;

int main(int argc, char** argv)
{
	const std::span<char*> args(argv, static_cast<std::size_t>(argc));

	simulation_parameters params;
	cell_type type;
	index_t samples = 301;

	try
	{
		if (args.size() != 1 && args.size() != 6 && args.size() != 7)
		{
			throw std::invalid_argument("expected 0, 5 or 6 arguments");
		}
		if (args.size() >= 6)
		{
			type.equilibrium_distance = std::stod(args[1]);
			params.steepness = std::stod(args[2]);
			type.intra_well_depth = std::stod(args[3]);
			type.inter_well_depth = std::stod(args[4]);
			params.cutoff_factor = std::stod(args[5]);
		}
		if (args.size() == 7)
		{
			samples = std::stoull(args[6]);
		}
		if (samples < 2)
		{
			throw std::invalid_argument("samples must be at least 2");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "potential_curve: " << e.what() << '\n'
				  << "usage: " << args[0]
				  << " [equilibrium_distance steepness intra_well_depth inter_well_depth cutoff_factor [samples]]\n";
		return EXIT_FAILURE;
	}

	const pair_potential intra = intra_potential(type, params);
	const pair_potential inter = inter_potential(type, params);

	// sample slightly past the cutoff so the truncation is visible
	const real_t r_max = 1.1 * intra.cutoff_distance;

	std::cout << std::setprecision(std::numeric_limits<real_t>::max_digits10);
	std::cout << "# equilibrium_distance=" << type.equilibrium_distance << " steepness=" << params.steepness
			  << " intra_well_depth=" << type.intra_well_depth << " inter_well_depth=" << type.inter_well_depth
			  << " cutoff_factor=" << params.cutoff_factor << '\n';
	std::cout << "r,intra_energy,intra_force,inter_energy,inter_force\n";

	for (index_t i = 0; i < samples; ++i)
	{
		const real_t r = r_max * static_cast<real_t>(i) / static_cast<real_t>(samples - 1);
		const real_t r2 = r * r;

		// radial force, positive = repulsive
		std::cout << r << ',' << potential_energy(intra, r2) << ',' << force_coefficient(intra, r2) * r << ','
				  << potential_energy(inter, r2) << ',' << force_coefficient(inter, r2) * r << '\n';
	}

	return EXIT_SUCCESS;
}
