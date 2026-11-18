---
layout: default
title: Mechanics Module
parent: Architecture
nav_order: 3
description: "Cell mechanics in PhysiCore: what a mechanics implementation provides and the two that exist"
---

# Mechanics Module
{: .no_toc }

The mechanics module moves cells: it computes the forces between them and with their surroundings and advances their positions. PhysiCore has two implementations, which differ in what a cell is made of.

## Table of Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Overview

Mechanical interactions act on timescales of minutes, slower than the diffusion of substrates and faster than cell division. A mechanics implementation:

- holds the positions of the cells (or of what cells are made of);
- computes the forces on them: repulsion, adhesion, confinement by the domain, motility;
- advances the positions by one mechanics timestep;
- writes the state for visualisation.

**Location:** `mechanics/`

## Physical Principle

At the scale of cells, inertia is negligible against viscous drag. Positions follow from a balance of forces:

$$
\eta\,\frac{d\mathbf{x}_i}{dt} = \sum \mathbf{F}_i
$$

with $\eta$ the drag coefficient. Both implementations advance this equation with an explicit scheme, so the timestep is limited by the stiffest interaction.

## What an Implementation Looks Like

There is no common base class for mechanics. The two implementations share a layout, the one of BioFVM:

| Part | Role |
|---|---|
| `environment` | owns the cells, the parameters and the solver. It is a `physicore::timestep_executor`: `run_single_timestep()` advances the mechanics by one step and `serialize_state()` writes the state |
| `environment_builder` | puts an environment together in code |
| `environment::create_from_config` | creates an environment from a PhysiCell settings file |
| `solver` | the contract of a backend: `initialize(environment&)` and `solve(environment&, iterations)` |
| `solver_registry` | backends register under a name and are looked up by it |
| `kernels/<backend>/` | a backend; both implementations have `openmp_solver` |
| `serializer` | output, with a VTK writer for ParaView |

```cpp
#include <common/timestep_executor.h>

namespace physicore {

class timestep_executor
{
public:
	virtual void run_single_timestep() = 0;
	virtual void serialize_state(real_t current_time) = 0;
};

} // namespace physicore
```

Cells are stored as structures of arrays in the containers of the common module, and reached through proxy objects; see [Common Module](Common.md).

## Implementations

### PhysiCell Mechanics

**Namespace:** `physicore::mechanics::physicell`

**Location:** `mechanics/physicell/`

One agent per cell: a sphere with a radius. It follows the mechanics of PhysiCell and reads the same settings file. From its code it provides:

- cell-cell repulsion and adhesion, with a strength per cell type and adhesion affinities between types;
- repulsion at the edges of the domain (`virtual_wall_at_domain_edge`);
- spring attachments between cells, with attachment and detachment rates;
- motility with a persistence time and a migration bias, including chemotaxis along substrate gradients, for which it takes the reactions-diffusion module as input.

This summary was written from the code of the module and is not its documentation.

### Micromechanics

**Namespace:** `physicore::mechanics::micromechanics`

**Location:** `mechanics/micromechanics/`

Many elements per cell (the sub-cellular element model). A cell is a cluster of elements held together by a short-range potential, so it has a shape, deforms against its neighbours, and grows and divides by gaining and splitting elements.

👉 **[Learn more about Micromechanics](Micromechanics.md)**

### Choosing

| | PhysiCell mechanics | Micromechanics |
|---|---|---|
| A cell is | one sphere | a cluster of elements |
| Cell shape | fixed | emerges, deforms |
| Cost per cell | one agent | tens to hundreds of elements |
| Suited to | many cells, tissue scale | fewer cells, where shape and contact matter |

## Coupling with Other Modules

### Reactions-Diffusion Module

- PhysiCell mechanics takes the reactions-diffusion module as input, for chemotaxis.
- Micromechanics is not coupled to it yet.

### Phenotype Module

The phenotype decides what cells do: grow, divide, die, change how they move or stick. Mechanics provides the operations and reports back, for example the pressure on a cell. See [Phenotype Module](Architecture-Phenotype.md).

## Next Steps

- **[Micromechanics Implementation](Micromechanics.md)** - Model, API and configuration of the sub-cellular element implementation
- **[Common Module](Common.md)** - Agent containers and structure-of-arrays storage
- **[Reactions-Diffusion Module](Architecture-Diffusion.md)** - Substrate transport

---

**See also:**
- [Architecture Overview](Architecture.md)
