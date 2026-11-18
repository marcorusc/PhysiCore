---
layout: default
title: Micromechanics Implementation
parent: Mechanics Module
nav_order: 1
description: "Guide to the sub-cellular element implementation of the mechanics module: model, API, configuration and kernels"
---

# Micromechanics Implementation
{: .no_toc }

Micromechanics is PhysiCore's sub-cellular element implementation of the mechanics module. A cell is not one particle but a cluster of elements held together by a short-range potential, so cells have a shape, deform against each other, and can grow and divide by gaining and splitting elements.

## Table of Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Model

The model follows Sandersius and Newman, *Phys. Biol.* **5** (2008) 015002.

### Equation of Motion

Every element $i$ moves without inertia:

$$
\eta\,\dot{\mathbf{x}}_i = \sum_{j} \mathbf{F}_{ij} + \mathbf{F}_i^{\text{wall}} + \eta\,v_c\,\mathbf{d}_c + \boldsymbol{\xi}_i
$$

- $\eta$ is the drag of one element, set by the type of its cell.
- $\mathbf{F}_{ij}$ is the force from element $j$, from the pair potential below.
- $\mathbf{F}_i^{\text{wall}}$ is the force of the walls of the domain, if they are switched on.
- $v_c\,\mathbf{d}_c$ is the migration velocity of the element's cell, if the cell is motile.
- $\boldsymbol{\xi}_i$ is thermal noise with $\langle \xi\,\xi \rangle = 2\,kT\,\eta$, if the thermal energy $kT$ is not zero.

Positions are advanced with the forward Euler scheme (Euler-Maruyama with noise).

### Pair Potential

$$
V(r) = u\left(t^2 - 2t\right), \qquad t = \exp\!\left[\rho\left(1 - \frac{r^2}{r_\text{eq}^2}\right)\right]
$$

The minimum $-u$ is at the equilibrium distance $r_\text{eq}$; the stiffness there is $8\rho^2 u / r_\text{eq}^2$. Interactions are cut off at `cutoff_factor` equilibrium distances.

### Three Kinds of Pairs

| Two elements of | Equilibrium distance | Well depth | Attraction multiplied by |
|---|---|---|---|
| the same cell, of type $t$ | $r_t$ | $u_t$ | 1 |
| two cells of the same type $t$ | $r_t$ | $c_t$ | $a_{tt}$ |
| two cells of the types $t$ and $s$ | $(r_t + r_s)/2$ | $\sqrt{c_t c_s}$ | $\sqrt{a_{ts} a_{st}}$ |

The well inside a cell ($u_t$) is deeper than the well between cells ($c_t$), which keeps a cell together. An adhesion affinity $a$, between 0 and 1, weakens only the attractive side of the potential: cells with the affinity 0 do not stick but still cannot overlap. The geometric means are those of PhysiCell.

### What Else the Model Has

| Feature | In short |
|---|---|
| Scaling | `scale_to_elements` turns the radius, stiffness and damping of a cell into element values for any number of elements (three dimensions) |
| Growth | `grow_cell` adds one element next to the cell |
| Division | `divide_cell` splits the elements of a cell by a plane, at any ratio |
| Motility | a persistent random walk of the cell: speed, persistence time, bias |
| Pressure | the virial of the forces from outside a cell, positive under compression |
| Walls | a soft repulsion at the faces of the domain, with a hard stop behind it |
| Reproducibility | all random numbers follow from one seed; results do not depend on the number of threads |

Every parameter, with hints on choosing it, is explained in `mechanics/micromechanics/docs/Parameters.md`. Each feature has a notebook in `mechanics/micromechanics/docs/notebooks/`.

## Structure

**Namespace:** `physicore::mechanics::micromechanics`

**Location:** `mechanics/micromechanics/`

The module follows the layout of BioFVM:

| Part | Where | What |
|---|---|---|
| Public API | `include/micromechanics/` | `environment`, `environment_builder`, `cell_type`, `simulation_parameters`, the cell and sub-agent interfaces, `solver`, `solver_registry`, `serializer` |
| Library | `src/` | data, containers, settings reader, VTK output |
| Kernel | `kernels/openmp_solver/` | the solver backend, registered as `"openmp_solver"` |
| Examples | `examples/` | small programs that print CSV or VTK |
| Documentation | `docs/` | parameter reference and notebooks |

### Data

Cells and their elements ("sub-agents") live in two containers of the shared `generic_agent_and_data_container`, each a structure of arrays:

- **Sub-agents** (`environment::agents`): position, velocity and the index of their cell.
- **Cells** (`environment::cells`): type, motility, and what the solver derives from the sub-agents (position, velocity, number of sub-agents, pressure).

A sub-agent refers to its cell by the cell's index. Containers remove an item by moving the last one into its place, so cells and sub-agents must be created and removed through `environment`, which keeps those indices valid.

### Solver

A backend implements `solver`:

```cpp
class solver
{
public:
	virtual void initialize(environment& e) = 0;
	virtual void solve(environment& e, index_t iterations) = 0;
};
```

`environment` is a `timestep_executor`: `run_single_timestep()` calls `solver->solve(*this, 1)`.

The OpenMP kernel does, in every iteration:

1. choose the motility direction of every motile cell;
2. sort the sub-agents into a uniform grid and compute the pair forces and the wall forces;
3. advance the positions, add the noise, and put back any sub-agent that left a walled domain.

After the last iteration it updates the cells: position, velocity, number of sub-agents and pressure.

## Using It

### From a Settings File

A simulation can be described by the PhysiCell settings file. Fields PhysiCell already has are reused (domain, `dt_mechanics`, the virtual wall option, the volume, adhesion affinities and motility of each cell definition, the file of cell positions). Each cell definition gets one new block, `<sub_cellular_elements>`.

```cpp
#include <micromechanics/environment.h>

using namespace physicore::mechanics::micromechanics;

auto env = environment::create_from_config("settings.xml");

for (int step = 0; step < 1000; ++step)
{
	env->run_single_timestep();
}
```

`mechanics/micromechanics/examples/settings.xml` is a complete small file with two cell types.

### In Code

```cpp
#include <array>

#include <micromechanics/cell_parameters.h>
#include <micromechanics/environment_builder.h>

using namespace physicore;
using namespace physicore::mechanics::micromechanics;

// a cell of radius 10 made of 64 elements, from its stiffness and damping
const cell_parameters cell { .radius = 10.0, .stiffness = 5000.0, .damping = 5000.0 };
cell_type epithelial = scale_to_elements(cell, 64);
epithelial.name = "epithelial";

environment_builder builder;
builder.resize(3, { -100, -100, -100 }, { 100, 100, 100 }, { 20, 20, 20 });
builder.set_time_step(0.002);
builder.set_domain_walls(true);
const index_t type = builder.add_cell_type(epithelial);

auto env = builder.build();

const std::array<real_t, 3> centre = { 0.0, 0.0, 0.0 };
env->create_cell_of_type(type, centre);

env->solver->solve(*env, 500);
```

`build()` refuses values that cannot work, among them a timestep beyond the stability limit of forward Euler; `largest_stable_timestep` returns that limit for a cell type.

### Several Cell Types

```cpp
cell_type first;
first.name = "first";
first.elements_per_cell = 19;
first.adhesion_affinities = { 1.0, 0.1 }; // to itself, to the second type

cell_type second = first;
second.name = "second";
second.adhesion_affinities = { 0.1, 1.0 };

builder.add_cell_type(first);
builder.add_cell_type(second);
```

### Growth and Division

```cpp
// one more element for the first cell
env->grow_cell(0);

// divide it across its longest axis; the daughter gets 30% of the elements
const auto axis = env->longest_axis(0);
cell_interface* daughter = env->divide_cell(0, axis, 0.3);
```

These are operations: when a cell grows or divides is decided by the code that drives the simulation.

### Reading the State

```cpp
for (index_t c = 0; c < env->cells->size(); ++c)
{
	cell_interface* a_cell = env->cells->get_agent_at(c);
	const auto position = a_cell->position(); // centre of its sub-agents
	const real_t pressure = a_cell->pressure();
	const index_t elements = a_cell->agent_count();
}

for (index_t i = 0; i < env->agents->size(); ++i)
{
	agent_interface* element = env->agents->get_agent_at(i);
	const auto position = element->position();
	const index_t owner = element->cell_index();
}
```

The cell values are current after `solver->initialize` and after every `solver->solve`.

### Output

```cpp
#include <micromechanics/vtk_serializer.h>

env->serializer = make_vtk_agents_serializer("output");
env->serialize_state(0.0);
```

This writes `output/micromechanics_agents.pvd` for ParaView: one point per sub-agent with `cell_index`, `cell_type` and `velocity`.

### From Another Project

```cmake
find_package(Micromechanics CONFIG REQUIRED)
target_link_libraries(my_simulation PRIVATE PhysiCore::mechanics.micromechanics)
```

## Limits

- **Forward Euler only.** The timestep must stay below a limit set by the stiffest cell type.
- **The scaling from cell to element values is for three dimensions.** In two dimensions the element values are given directly.
- **The neighbour grid is built by one thread.** The force and position loops run in parallel.
- **One kernel**, for the CPU with OpenMP.
- **The domain is a box with whole-number bounds**; walls are the only boundary condition.
- **Cells have one compartment.** A nucleus or other compartments are not modelled.
- **The list of neighbouring cells** in the cell interface is not filled.
- **No coupling yet** to the reactions-diffusion module or to cell phenotype.

## Further Reading

- `mechanics/micromechanics/docs/Parameters.md` - every parameter and how to choose it
- `mechanics/micromechanics/docs/notebooks/` - one notebook per feature, with the checks behind the numbers
- [BioFVM Implementation](BioFVM.md) - the module whose layout this one follows
