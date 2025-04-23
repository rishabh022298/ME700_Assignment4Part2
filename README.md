# ME700_Assignment4Part2


Problem related to **thermal expansion with linear elasticity** was implemented. The most of the code structure was borrowed from the [hyperelasticity example](https://github.com/Lejeune-Lab-Graduate-Course-Materials/fenicsX/blob/main/hyperelasticity_beam.py).

## Instructions for running the script

### Installing FEniCSx on the SCC

```bash
module load miniconda
mamba create -n fenicsx-env
mamba activate fenicsx-env
mamba install -c conda-forge fenics-dolfinx mpich pyvista
pip install imageio
pip install gmsh
pip install PyYAML
```

### Running the scripts on VSCode Server
Launch VSCode Server.\
Open terminal in VSCode.\
Run the following commands one by one in the termainal to clone the repository (after moving to the desired directory):

```bash
git clone https://github.com/rishabh022298/ME700_Assignment4Part2.git
```
Change the folder:
```bash
cd ME700_Assignment4Part2
```
Activate fenicsx-env:
```bash
conda activate fenicsx-env
```
**Note:** If you are running into some memory related error then try relaunching the server after closing the session.

### File names

#### Part 1: Analytical vs Numerical
This can be found [here](https://github.com/rishabh022298/ME700_Assignment4Part2/blob/main/p1_analytical_vs_num.py)
```bash
python p1_analytical_vs_num.py
```

#### Part 2: h and p refinement
Scipt for h-refinement can be found [here](https://github.com/rishabh022298/ME700_Assignment4Part2/blob/main/p2_h_refinement.py)
```bash
python p2_h_refinement.py
```
Script for p-refinement can be found [here](https://github.com/rishabh022298/ME700_Assignment4Part2/blob/main/p2_p_refinement.py)
```bash
python p2_p_refinement.py
```

#### Part 3: FEA Code Failure
Script for first example using boundary free system can be found [here](https://github.com/rishabh022298/ME700_Assignment4Part2/blob/main/p3_free_boundaries.py)
```bash
python p3_free_boundaries.py
```

Script for second example using poor mesh can be found [here](https://github.com/rishabh022298/ME700_Assignment4Part2/blob/main/p3_poor_mesh.py)
```bash
python p3_poor_mesh.py
```

![Alt Text](figures/P1_analytical_vs_numerical.png)
