# Dimer calculations

`dimer_calculations` is a Python toolkit for generating molecular dimers and
preparing them for energy calculations. It was developed for porous organic
and metal-organic cages, but its axis-based dimer generation can be applied to
other molecules.

The package can:

- define molecular axes from SMARTS, SMILES, molecular fragments, or
  geometric midpoints;
- generate displaced, rotated, and optionally laterally shifted dimers;
- reject configurations with atomic overlap or cage catenation; and
- optimise accepted dimers with xTB, GULP, MacroModel, or OpenMM.

## Installation

Python 3.11 or newer is required. We recommend a clean Conda environment:

```bash
conda create --name dimer-calculations python=3.11
conda activate dimer-calculations
python -m pip install .
```

To run the example notebook, install the notebook dependencies as well:

```bash
python -m pip install ".[notebook]"
jupyter lab dimer_calculations.ipynb
```

Once version 1.0.0 has been released, it can also be installed directly from
GitHub:

```bash
python -m pip install \
  "git+https://github.com/ewolpert1/dimer_calculations.git@v1.0.0"
```

## External optimisation programs

Dimer generation and screening do not require an external optimiser. Energy
optimisation requires at least one of the following separately installed
programs:

- [xTB](https://xtb-docs.readthedocs.io/);
- [GULP](https://gulp.curtin.edu.au/); or
- [MacroModel](https://www.schrodinger.com/platform/products/macromodel/).

The package locates `xtb` and `gulp` automatically when they are on `PATH`.
Otherwise, set the appropriate environment variable before starting Python:

```bash
export XTB_PATH=/path/to/xtb
export GULP_PATH=/path/to/gulp
export SCHRODINGER_PATH=/path/to/schrodinger
```

The resolved values are available as `XTB_PATH`, `GULP_PATH`, and
`SCHRODINGER_PATH` from `dimer_calculations.config`. Optimiser functions also
accept executable paths directly.

## OpenMM

`optimise_dimer_openmm` optimises a dimer with an
[OpenFF](https://openforcefield.org/) force field in
[OpenMM](https://openmm.org/). It runs inside Python, so it needs no
executable path, but the OpenFF packages are installed from conda-forge:

```bash
conda install -c conda-forge openmm openff-toolkit openff-interchange \
  openff-nagl openff-nagl-models
```

The default force field, Sage 2.3.0 (`openff_unconstrained-2.3.0.offxml`),
assigns partial charges with NAGL in seconds. Like the other optimisers, it
writes the optimised dimer to `{output_dir}_opt.mol`, and
`{output_dir}/openmm_opt.out` records the energy before and after
optimisation.

Atoms in `fixed_atom_set` keep their positions, which holds the two cages at
the intended offset. Atom ids count from 0, as `Cage.fix_atom_set` returns
them:

```python
from dimer_calculations import cage, optimiser_functions

fixed_atom_set = cage.Cage.fix_atom_set(dimer, "NCCN", metal_atom=None)
optimiser_functions.optimise_dimer_openmm(
    dimer=dimer,
    output_dir="OpenMM_shell_0_slide_0_rot_0",
    fixed_atom_set=fixed_atom_set,
)
```

This has been tested with OpenMM 8.2, OpenFF Toolkit 0.18, OpenFF Interchange
0.5 and OpenFF NAGL 0.5.

## Basic workflow

1. Load an `stk.Molecule` and place its centroid at the origin.
2. Use `dimer_calculations.axes` to define the molecular axes of interest.
3. Pass two molecules and a pair of axes to
   `dimer_calculations.dimer_generator.DimerGenerator`.
4. Screen the generated structures for overlap and, for cages, catenation.
5. Pass accepted structures to an optimiser in
   `dimer_calculations.optimiser_functions`.

The complete worked example is provided in
[`dimer_calculations.ipynb`](dimer_calculations.ipynb). The example molecular
structures are in [`cages/`](cages/).

## Development and tests

```bash
python -m pip install ".[dev]"
python -m pytest
```

The smoke test generates a dimer without invoking an external optimiser. The
OpenMM tests are skipped unless OpenMM and OpenFF are installed.

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). When using an
archived release, cite the DOI for that specific software version.

## Licence

This project is released under the [MIT Licence](LICENSE).
