# Dimer calculations

`dimer_calculations` is a Python toolkit for generating molecular dimers and
preparing them for energy calculations. It was developed for porous organic
and metal-organic cages, but its axis-based dimer generation can be applied to
other molecules.

The package can:

- define molecular axes from SMARTS, SMILES, molecular fragments, or
  geometric midpoints;
- generate displaced, rotated, and optionally laterally shifted dimers;
- reject configurations with atomic overlap or cage catenation;
- optimise accepted dimers with xTB, GULP, or MacroModel; and
- score dimers with an OpenFF force field in OpenMM, giving interaction and
  binding energies.

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

## OpenMM energies

`dimer_calculations.openmm_energies` scores dimers with an
[OpenFF](https://openforcefield.org/) force field in
[OpenMM](https://openmm.org/), in vacuum, without an external program. The
default force field is Sage 2.3.0 (`openff_unconstrained-2.3.0.offxml`), whose
NAGL partial charges take milliseconds per molecule.

The OpenFF stack is installed from conda-forge rather than pip:

```bash
conda install -c conda-forge openmm openff-toolkit openff-interchange \
  openff-nagl openff-nagl-models
```

It has been tested with OpenMM 8.2, OpenFF Toolkit 0.18, OpenFF Interchange 0.5
and OpenFF NAGL 0.5.

```python
from dimer_calculations.openmm_energies import (
    OpenMMCalculator,
    dimer_energies,
    monomer_energies,
)

calculator = OpenMMCalculator()
reference = monomer_energies(molecule, calculator)

for entry in list_of_dimers:  # from DimerGenerator.generate
    energies = dimer_energies(entry["Dimer"], reference, reference, calculator)
    print(entry["Displacement shell"], entry["Rotation"], energies.e_bind)
```

For a dimer of two different molecules, pass each molecule's own
`monomer_energies`. The monomers must be the molecules the dimers were built
from, so that each cage in a dimer is a rigid copy of its monomer.
`energies.minimised` holds the minimised dimer, and `energies.as_dict()`
returns every energy under the column names below. All energies are in kJ/mol;
a negative interaction or binding energy means the two cages attract.

| Term | Definition | Meaning |
|---|---|---|
| `E_sp` | | The dimer exactly as built, without geometry optimisation. |
| `E_min` | | The dimer after minimisation. |
| `E_int_sp` | E_sp − E(cage 1 as given) − E(cage 2 as given) | The interaction between the two cages without optimisation. |
| `E_int_min` | E_min − E(cage 1 in the minimised dimer) − E(cage 2 in the minimised dimer) | The same interaction, but after minimisation. |
| `E_bind` | E_min − E_min(cage 1 alone) − E_min(cage 2 alone) | The binding energy: the energy released when two relaxed cages form the relaxed dimer. It is the best single number for comparing motifs. |
| `E_def` | E_bind − E_int_min | The deformation energy: how much more internal energy the cages have in the minimised dimer than when relaxed alone. |

Do not compare `E_sp` with `E_min` directly: both include each cage's own
internal energy, and the cages relax during minimisation.

## Basic workflow

1. Load an `stk.Molecule` and place its centroid at the origin.
2. Use `dimer_calculations.axes` to define the molecular axes of interest.
3. Pass two molecules and a pair of axes to
   `dimer_calculations.dimer_generator.DimerGenerator`.
4. Screen the generated structures for overlap and, for cages, catenation.
5. Pass accepted structures to an optimiser in
   `dimer_calculations.optimiser_functions`, or score them with
   `dimer_calculations.openmm_energies`.

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
