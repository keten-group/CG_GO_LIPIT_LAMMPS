### Overview

**This model was developed by the Keten Research Group at Northwestern University**

This repository contains all files required to generate structures and run LIPIT-like ballistic impact simulations on multilayered graphene oxide (GO) structures represented by a coarse-grained molecular model with a 4-to-1 mapping relative to atomistic graphene. The model is intended to be run using [LAMMPS](https://www.lammps.org/#gsc.tab=0) molecular dynamics software. **If one or more of these files contibutes to your publication, please cite _Meng et al. 2017_ (at minimum). Note that model parameters have been adjusted in subsequent publications.** Select publications are listed at the end of this document. Figures within this document were created using [OVITO](https://www.ovito.org/): A. Stukowski, Modelling Simul. Mater. Sci. Eng. 18, 015012 (2010).

Additionally, this repository includes a spreadsheet (GO_thickness_data.xlsx) containing calculated specific penetration energies ($E_p^*$) across a range of GO flake sizes, film thicknesses, and impact velocities. This data accompanies publication (1) in the list below, which is currently under review.

<img src='https://drive.google.com/uc?id=1qxH2H8M85IE2fGxyyo3RAiRTflXIp_57' width="1024" height="768">

### About the model

The model maps 4 atoms to 1 CG bead and assigns a bead type based on oxidation type. It has been shown to accurately reproduce in-plane mechanical properties and interlayer adhesion strength of GO at varying degrees of oxidation (see publication (4) in the list at the end of this file). 

**Masses**
- Bead types
  - Type 1: 48.0 g/mol (unoxidized bead)
  - Type 2: 65.0 g/mol (hydroxyl-oxidized bead)
  - Type 3: 64.0 g/mol (epoxide-oxidized bead)
  - Type 4: 96.0 g/mol (projectile bead) 
- NOTE: Type 4 is not necessary for equilibration, but needs to be added to the data file prior to LIPIT

**Bonds**
- Bond style: [hybrid](https://docs.lammps.org/bond_hybrid.html) of [Morse](https://docs.lammps.org/bond_morse.html) and [table](https://docs.lammps.org/bond_table.html)
- Parameters
  - All beads are type 1: Morse, $D$ = 443.07 kcal/mol, $\alpha$ = 1.154 &#197;<sup>-1</sup>, $r_0$ = 2.86 &#197;      
  - At least one bead is type 2, but none are type 3: table, file _bond2-GO.table_
  - At least one bead is type 3: table, file _bond3-GO.table_
                         
**Angles**
- Angle style: [harmonic](https://docs.lammps.org/angle_harmonic.html)
- Parameters
  - All beads are type 1: $K$ = 456.61 kcal/mol-&#197;<sup>2</sup>, $\theta_0$ = 120&deg;
  - At least one bead is type 2, but none are type 3: $K$ = 259.47 kcal/mol-&#197;<sup>2</sup>, $\theta_0$ = 120&deg;
  - At least one bead is type 3: $K$ = 189.93 kcal/mol-&#197;<sup>2</sup>, $\theta_0$ = 120&deg;

 **Pairs**
 - Pair style: [hybrid]() of [lj/cut](https://docs.lammps.org/pair_lj.html) and [lj/charmm/coul/charmm](https://docs.lammps.org/pair_charmm.html)
 - Parameters
   - Type 1 - Type 1: lj/cut, $\varepsilon$ = 0.204 kcal/mol, $\sigma$ = 7.48 &#197;
   - Type 2 - Type 2: lj/cut, $\varepsilon$ = 1.024 kcal/mol, $\sigma$ = 7.48 &#197;
   - Type 3 - Type 3: lj/cut, $\varepsilon$ = 0.638 kcal/mol, $\sigma$ = 7.48 &#197;
   - Type 4 - Type 4: lj/charmm/coul/charmm, $\varepsilon$ = 0.813 kcal/mol, $\sigma$ = 3.46 A
   - Cross interactions are described by Lorentz-Berthelot mixing rules using the LAMMPS command  [pair_modify mix arithmetic](https://docs.lammps.org/pair_modify.html)
 - NOTE: Type 4 is not necessary for equilibration, but needs to be used for the projectile in the _lipit.in_ input script

### About the LAMMPS protocols

Measurements acquired via these LIPIT-like simulations, like those of LIPIT itself, are quite sensitive to geometric choices, particularly relative scales of the projectile and the film. Fore more information, we reccomend reviewing this publication: Z. Meng and S. Keten, “Unraveling the Effect of Material Properties and Geometrical Factors on Ballistic Penetration Energy of Nanoscale Thin Films,” Journal of Applied Mechanics, vol. 85, no. 12, Sept. 2018, doi: 10.1115/1.4041041.

Detailed descriptions of the LAMMPS protocols can be found in publications (1) and (2) in the list at the end of this file. The _*.table_ files must be in the same directory as the input script, or otherwise the path within the input script should be changed. 

Files: 
- (1) equil.in
- (2) lipit.in
- (3) bond2-GO.table
- (4) bond3-GO.table

### About the MATLAB structure generation script

File: CG_graphene_oxide_generation.m

This MATLAB script uses user inputs of "nx" and "ny" to describe the number of units in the zigzag and armchair directions, respectively, for each GO flake. Other user specifications include number of layers (corresponding to thickness), number of flakes per layer in x and y, percents of beads functionalized by hydroxyl and epoxide groups, and fraction of flake overlap between layers. 

The below image demonstrates how the number of units is counted in each direction. Comments within the script contain equations for estimating the number of units needed for a given flake edge length.

<img src='https://drive.google.com/uc?id=15AQClNzgdR7renqEzswIjGi6LBqTADnq' width="620" height="600">

### Select publications (reverse chronological):
1) UNDER REVIEW: H. L. White, W. Chen, N. M. Pugno, and S. Keten, "Rate-dependent size effects govern the inverse thickness dependence of specific penetration energy in nanoscale thin films"
2) H. L. White, A. Giuntoli, M. Fermen-Coker, and S. Keten, “Tailoring flake size and chemistry to improve impact resistance of graphene  oxide thin films,” Carbon, vol. 215, p. 118382, Nov. 2023,  doi: 10.1016/j.carbon.2023.118382.
3) T. Li, Z. Meng, and S. Keten, “Interfacial mechanics and viscoelastic properties of patchy graphene oxide reinforced nanocomposites,” Carbon, vol. 158, pp. 303–313, Mar. 2020, doi: 10.1016/j.carbon.2019.10.039.
4) Z. Meng et al., “A coarse-grained model for the mechanical behavior of graphene oxide,” Carbon, vol. 117, pp. 476–487, June 2017, doi: 10.1016/j.carbon.2017.02.061.


For questions not answered by the provided materials, please contact the Keten Lab (https://keten-group.northwestern.edu/).
