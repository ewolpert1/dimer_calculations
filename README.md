# Dimer calculations

A toolkit to help produce and calculate the energy of dimers of molecules. This code was originally built with cages in mind but can be used for more general cases. The creation of dimers relies on [stk](https://stk.readthedocs.io/en/stable/), which comes with the pip install

# Installation

The code can be installed following these steps:

1. Create a `conda` environement:
```
conda create --name your-env-name python=3.11
```

2. Activate the environment:
```
conda activate your-env-name
```

3. Clone the package

```
git clone -b refactoring --single-branch https://git@github.com/ewolpert1/dimer_calculations.git
```

4. Install the package

```
cd dimer_calculations
pip install .
```

# Usage

There are currently three different optimisers set up to work with the code; `OPLS` through [Macromodel](https://www.schrodinger.com/platform/products/macromodel/), [GULP](https://gulp.curtin.edu.au/), and [XTB](https://xtb-docs.readthedocs.io/en/latest/optimization.html).

The code works by first defining the axes that the molecules are displaced along. There are four options of how to define the axes:
1. Using [Pywindow](https://github.com/marcinmiklitz/pywindow)
2. Using Smarts string, where the axis is the centroid of the Smarts string molecule
3. Using Smiles string, where the axis is the centroid of the Smiles string molecule
4. Using the midpoint, where the axis is defined by the midpoint of the vectors provided.

Using these axes you can generate the dimers as a dictionary in which you can than optimise the structures.

Example scripts for using both `OPLS` and `XTB` are provided in the example_scripts folder labelled G17.py and WaaF-1.py respectively.

