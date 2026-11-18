# Parameters of the sub-cellular element model

This page lists every parameter of the micromechanics module: what it means physically, where it is set,
what happens when it is raised or lowered, and how to choose it. The hints come from the tests and
notebooks of this module unless a paper is named; each one says which.

A simulation is described by a PhysiCell settings file (`PhysiCell_settings.xml`) and created with
`environment::create_from_config`. Fields that PhysiCell already has are reused where they mean the same
here. Two blocks are specific to this model: `<sub_cellular_elements>` in every cell definition and an
optional `<micromechanics>` at the top level. `examples/settings.xml` is a complete small example.

All values must be in **one consistent set of units**; the `units` attributes of the file are not
interpreted. The examples use micrometres, minutes and piconewtons, so a stiffness is in pN/µm, a damping
in pN·min/µm and an energy in pN·µm.

## The model in one paragraph

A cell is a cluster of $N$ elements. Two elements interact through the potential

$$V(r) = u\left(t^2 - 2t\right), \qquad t = \exp\!\left[\rho\left(1 - \frac{r^2}{r_\mathrm{eq}^2}\right)\right]$$

which has its minimum $-u$ at the distance $r_\mathrm{eq}$ and the stiffness
$\kappa = 8\rho^2 u / r_\mathrm{eq}^2$ there. Elements of the same cell use a deeper well ($u$) than
elements of different cells ($c$), which is what keeps a cell together. Every element moves with
$\eta\,\dot{\mathbf x} = \mathbf F + \text{noise}$, plus a drift if its cell is motile.

## Per cell type

These are set in each `<cell_definition>`. In the code they are the fields of `cell_type`.

### Size: `<phenotype><volume><total>`

| | |
|---|---|
| Meaning | Volume of the cell. The cell is taken as a sphere of that volume, radius $R = (3V/4\pi)^{1/3}$ |
| Reused from | PhysiCell |
| How to choose | From the measured cell volume or diameter |

- With the default steepness a relaxed cell comes out **7 to 8% smaller in radius** than requested (9.17 to
  9.32 for a requested 10, for 128 to 1024 elements). If the exact size matters, request a slightly larger
  volume. *(cell scaling notebook)*

### Resolution: `<sub_cellular_elements><elements>`

| | |
|---|---|
| Symbol | $N$ |
| Meaning | Number of elements that make up one cell of this type |
| Code | `cell_type::elements_per_cell` |

- **More elements:** finer shape, smaller elements ($r_\mathrm{eq} \propto N^{-1/3}$), cost growing in
  proportion to $N$, and a smaller stable timestep.
- **How to choose:** the smallest number that resolves the shapes you care about. Sandersius and Newman
  recommend about 100 or more for the scaling relations to hold; from 128 up the radius of a cell changed by
  less than 2% in our runs. Cells of 13 to 37 elements, as in several examples, keep together and interact
  sensibly but are too coarse to represent a measured stiffness. *(cell scaling notebook)*
- Stiffness, damping and volume describe the whole cell and do not need to change when $N$ does.

### Stiffness: `<sub_cellular_elements><stiffness>`

| | |
|---|---|
| Symbol | $\kappa_0$ |
| Meaning | How hard the whole cell resists deformation |
| Units | force per length |
| Code | `cell_parameters::stiffness`; it sets `cell_type::intra_well_depth` |

- The stiffness of one pair of elements follows as
  $\kappa = \kappa_0\,N^{-1/3}\left(1 - \lambda N^{-1/3}\right)$ with $\lambda = 0.75$, and the well depth as
  $u = \kappa\,r_\mathrm{eq}^2 / (8\rho^2)$ (Sandersius and Newman 2008, Eq. 11 and 5).
- **Higher:** the cell deforms less under the same load and relaxes faster; it needs a smaller timestep.
  In the sorting runs, stiffer cells rearranged slightly less. *(cell types notebook)*
- **How to choose:** from a measured cell stiffness. The elastic modulus of the cell is of the order
  $\kappa_0 / (2R)$. That relation is taken from the paper as an order of magnitude and has not been checked
  here against a simulated compression test. The paper's example cell uses $\kappa_0 = 5\times10^{-3}$ N/m.

### Damping: `<sub_cellular_elements><damping>`

| | |
|---|---|
| Symbol | $\eta_0$ |
| Meaning | How slowly the whole cell responds: its viscous resistance |
| Units | force × time per length |
| Code | `cell_parameters::damping`; it sets `cell_type::drag_coefficient` $= \eta_0 / N$ |

- **Higher:** everything the cell does takes longer in proportion. Shapes at rest do not change.
- The cell relaxes over a time of the order $\eta_0 / \kappa_0$. **How to choose:** from a measured recovery
  time after a deformation, together with the stiffness.
- It also sets how hard a motile cell pulls: its total traction is $\eta_0 \times$ migration speed.
- With thermal noise, an element diffuses with $D = kT / (\eta_0/N)$ and the centre of a cell with
  $kT/\eta_0$.

### Adhesion between cells of the type: `<sub_cellular_elements><adhesion_ratio>`

| | |
|---|---|
| Meaning | Depth of the well between two cells of this type, as a fraction of the well inside a cell |
| Default | 0.5 |
| Code | `cell_parameters::adhesion_ratio`; it sets `cell_type::inter_well_depth` $c$ |

- It sets both how strongly two cells stick and how hard they are to push into each other.
- **0.5:** two cells adhere and stay distinct. **1:** nothing tells two cells apart and they merge.
  **0:** the cells do not feel each other at all and pass through one another; such a type does not feel
  the walls either. *(two cells notebook)*
- **How to choose:** keep it clearly below 1. To weaken sticking without softening the contact, lower the
  affinity instead.

### Adhesion affinities: `<mechanics><cell_adhesion_affinities>`

| | |
|---|---|
| Symbol | $a_{ts}$, between 0 and 1 |
| Meaning | How much of the attraction to a cell of another type (or of the same type) is used |
| Reused from | PhysiCell, same tags; a type that is not listed has the affinity 1 |
| Code | `cell_type::adhesion_affinities` |

- An affinity weakens **attraction only**. Cells repel in full whatever their affinity, so two types with
  the affinity 0 do not stick but still cannot overlap.
- The value used for a pair of types is the geometric mean $\sqrt{a_{ts}\,a_{st}}$ of the two entries, as in
  PhysiCell. If either entry is 0 the pair does not attract.
- The entry of a type for itself scales how cohesive a tissue of that type is.
- **How to choose:** this is the knob for mixing. With 1 between two types they stay mixed; with 0.1 they
  sort (differential adhesion). In a packed 3D block the fraction of unlike neighbours fell from 0.46 to
  0.19 with an affinity of 0.1 and enough noise, and stayed at 0.51 with equal affinities.
  *(cell types notebook)*

### Motility: `<phenotype><motility>`

Reused from PhysiCell. A motile cell crawls as a persistent random walk; every element of the cell gets the
same extra velocity.

| Tag | Symbol | Meaning |
|---|---|---|
| `<speed>` | $v$ | speed of a cell that nothing holds back |
| `<persistence_time>` | $\tau$ | mean duration of a straight run |
| `<migration_bias>` | $b$ | 0: random directions; 1: always the bias direction |
| `<options><enabled>` | | whether cells of this type crawl |

- **Speed:** what cell tracking gives over short intervals. A crowded cell moves more slowly, since it only
  pushes with a finite traction.
- **Persistence time:** read it off where the mean squared displacement of tracked cells bends from $t^2$ to
  $t$. Over long times only $v^2\tau$ matters: the effective diffusion coefficient is $v^2\tau/d$ in $d$
  dimensions. Keep the timestep well below $\tau$.
- **Bias:** the drift along the bias direction is $v\,m(b)$ with a nonlinear $m$; in 2D a bias of 0.5 already
  gives 64% of the full speed.
- The bias direction itself (`migration_bias_direction` of a cell) is not in the file. It is set by the code
  that drives the simulation, for instance from a chemical gradient. PhysiCell's chemotaxis tags are not read.

*(motility notebook)*

### Giving the elements directly: `<sub_cellular_elements><element>`

Instead of stiffness and damping of the whole cell, the values of one element can be given. This is needed
in two dimensions, where the scaling above does not hold, and is convenient for model studies in arbitrary
units.

| Tag | Symbol | Meaning |
|---|---|---|
| `<equilibrium_distance>` | $r_\mathrm{eq}$ | distance at which two elements are at rest; an element is a sphere of radius $r_\mathrm{eq}/2$ |
| `<intra_well_depth>` | $u$ | depth of the well between two elements of one cell |
| `<inter_well_depth>` | $c$ | depth of the well between two cells of the type |
| `<drag_coefficient>` | $\eta$ | drag of one element |

With this form `<volume><total>` is not used: the size of a cell follows from `<elements>` and
$r_\mathrm{eq}$.

## For the whole simulation

### Steepness: `<micromechanics><steepness>`

| | |
|---|---|
| Symbol | $\rho$ |
| Default | 2, the value of Sandersius and Newman (2008) and of SEM2 |
| Code | `simulation_parameters::steepness` |

- It sets the range of the potential. **Higher:** shorter ranged, so elements feel fewer than their nearest
  neighbours' neighbours.
- At a given cell stiffness the well depth goes as $1/\rho^2$, so a steeper potential means **shallower
  wells and bonds that thermal noise breaks more easily**.
- With 4 a relaxed cell kept its requested radius exactly, but growing cells became lumpy.
  *(cell scaling and growth notebooks)*
- **How to choose:** leave it at 2 unless there is a reason. It is shared by all cell types.

### Cutoff: `<micromechanics><cutoff_factor>`

| | |
|---|---|
| Default | 2.5, the value of SEM2 |
| Meaning | Two elements do not interact beyond this many equilibrium distances |

- **Lower:** faster, and less accurate, as weak attraction at a distance is dropped. An element that gets
  further than this from everything is lost: nothing pulls it back.

### Thermal energy: `<micromechanics><thermal_energy>`

| | |
|---|---|
| Symbol | $kT$ |
| Default | 0: no noise, a deterministic simulation |
| Units | energy, the same as the well depths |
| Code | `simulation_parameters::thermal_energy` |

- Agitation of the elements, standing for thermal and active fluctuations of the cytoskeleton. It is one
  energy for all cell types; an element then diffuses with $kT/\eta$.
- **What matters is $kT$ against the well depths** $u$ and $c$ of the elements, which for cells given by
  stiffness and damping are computed, not typed in; the `scaled_cell` example prints them.
- **How to choose:**
  - below about a tenth of $c$: elements jiggle, nothing rearranges;
  - $kT \approx 0.4$ to $0.6\,c$ (0.2 to 0.3 for $c = 0.5$, $u = 1$): cells flow past each other and sort,
    with only a few elements lost over 200 time units;
  - a single bond of depth $u$ already breaks at $kT \approx 0.16\,u$, and elements that leave do not come
    back, so check the number of strayed elements for your own parameters.
- *(thermal noise and cell types notebooks)*

### Seed: `<micromechanics><random_seed>`

The same seed gives the same simulation bit for bit, whatever the number of threads. It feeds the thermal
noise, the turns of motile cells and the placement of new elements when a cell grows.

### Walls: `<options><virtual_wall_at_domain_edge>`

| | |
|---|---|
| Reused from | PhysiCell |
| Default | off if the tag is absent |
| Code | `simulation_parameters::domain_walls` |

- With walls, the faces of the `<domain>` push elements back in. A wall is as hard as a contact with a cell
  of the element's own type and never attracts. *(domain walls notebook)*
- Without walls the size of the domain plays no role.
- The domain has whole-number bounds, and cells must start inside it.

### Timestep: `<overall><dt_mechanics>`

Reused from PhysiCell. This is a numerical parameter, not a physical one: results must not depend on it.

- Forward Euler is stable only below

  $$\Delta t_\mathrm{max} = \frac{\eta}{4\kappa} \;\text{(3D)}, \qquad \frac{\eta}{3\kappa}\;\text{(2D)}$$

  with $\eta$ the drag and $\kappa = 8\rho^2 \max(u, c)/r_\mathrm{eq}^2$ the stiffness of the stiffest pair of
  the type. Beyond it the elements of a cell oscillate with growing amplitude and the cell flies apart.
- `environment_builder` and `create_from_config` **refuse a timestep at or above that limit** and name the
  cell type and the limit. `largest_stable_timestep` returns it.
- **How to choose:** a fifth of the limit or less. Stiffer cells, more elements per cell and smaller cells all
  lower the limit.

## A way to calibrate

1. **One cell.** Set volume, stiffness and damping from measurements, choose $N$, and relax a single cell.
   Check its radius and how long it takes to settle. *(cell scaling notebook, `scaled_cell` example)*
2. **Two cells.** Bring two cells of one type into contact and choose the adhesion ratio so that they stick
   and stay distinct. *(two cells notebook)*
3. **Timestep.** Take a fifth of the limit the builder reports, and check that halving it does not change the
   result.
4. **Several types.** Set the affinities between types; they do not disturb steps 1 to 3.
5. **Motion.** Add motility from tracking data, or thermal energy, and check that cells neither freeze nor
   shed elements.

## What is not read

The settings file may contain everything PhysiCell needs; this module ignores what it does not use: the
cycle, death and volume dynamics of the phenotype, secretion, PhysiCell's own adhesion and repulsion
strengths and spring attachments, chemotaxis, and the microenvironment. Some of these will be used once
cell behaviour is coupled to the mechanics.

## References

- Sandersius S A and Newman T J 2008 Modeling cell rheology with the subcellular element model
  *Phys. Biol.* **5** 015002
- Newman T J 2005 Modeling multicellular systems using subcellular elements *Math. Biosci. Eng.* **2** 613
- Chattaraj S, Torre M, Kalcher C, Stukowski A, Morganti S, Reali A and Pasqualini F S 2023 SEM2:
  introducing mechanics in cell and tissue modeling using coarse-grained homogeneous particle dynamics
  *APL Bioeng.*
- Ghaffarizadeh A, Heiland R, Friedman S H, Mumenthaler S M and Macklin P 2018 PhysiCell: an open source
  physics-based cell simulator for 3-D multicellular systems *PLoS Comput. Biol.*
