# SecChem

[![CI](https://github.com/AndBrn743/SecChem/actions/workflows/linux.yml/badge.svg?branch=master)](https://github.com/AndBrn743/SecChem/actions)
[![CI](https://github.com/AndBrn743/SecChem/actions/workflows/windows.yml/badge.svg?branch=master)](https://github.com/AndBrn743/SecChem/actions)
[![CI](https://github.com/AndBrn743/SecChem/actions/workflows/macos.yml/badge.svg?branch=master)](https://github.com/AndBrn743/SecChem/actions)
[![codecov](https://codecov.io/gh/AndBrn743/SecChem/branch/master/graph/badge.svg)](https://codecov.io/gh/AndBrn743/SecChem)

**SecChem** is a C++17 chemistry *infrastructure* library that provides core data structures and parsers for
representing **molecules** and **Gaussian-type basis sets**, intended to be reused by quantum chemistry (QC)
software.

This repository was split out of the author’s HF/DFT codebase (**SecScf**) to allow focused refactoring,
testing, and experimentation with cleaner abstractions and better software structure.


## Project goals

SecChem focuses on **representation**, not computation.

The primary goals are:

- Strongly typed representations of:
    - Chemical elements
    - Atoms and molecules
    - Gaussian basis sets and ECPs
- Parsing of:
    - Cartesian coordinate input
    - Z-matrix (internal coordinate) input
    - Basis set definitions (example: BSE JSON)
- Explicit data models suitable for **binding to QC backends** (downstream should implement binding themself)
- Clean, readable, pedantic, standard-conforming C++17

This library is intended to be used as **infrastructure** by electronic-structure codes, not as an end-user
program.


## Design philosophy

**Correctness > performance.**

More specifically:

- Correctness, clarity, and explicitness take priority over micro-optimizations
- Readability and maintainability are valued over cleverness
- “Performance tricks” (when needed) should be encapsulated behind clean interfaces
- No hidden global state
- Strong value/reference semantics where appropriate
- Deterministic, predictable behavior

Performance matters, but *only after the model is correct*.


## Non-goals

SecChem deliberately does **not** implement:

- Hartree-Fock, DFT, or post-HF methods
- Wavefunction analysis
- Integral evaluation
- Geometry optimization
- Molecular topology inference or graph algorithms
- Visualization
- Binding for other languages and downstream libraries

Those responsibilities belong to downstream libraries or applications.


## Status

- **Maturity**: Experimental
- **API & ABI stability**: Not guaranteed
- **Primary use**: Personal research infrastructure
- **Contributions**: Not actively solicited, but suggestions and contributions are welcome

The API is expected to evolve as the design is refined.


## Relationship to SecScf

SecChem was originally part of **SecScf**, a Hartree-Fock / DFT codebase.

The original implementation grew organically and was not ideal in terms of structure or testability.
This repository exists to:

- Refactor core chemistry representations
- Improve test coverage
- Experiment with better abstractions
- Make the code public and inspectable

If others find it useful, that’s a bonus.


## Build requirements

- **C++ standard**: C++17 (minimum)
- **Build system**: CMake
- **Required dependency**:
    - A fully C++17 compliant STL (sadly, some people do use ancient STL)
    - [Eigen3](https://gitlab.com/libeigen/eigen/) v5 (fetched by CMake if not found)
    - [range-v3](https://github.com/ericniebler/range-v3) v0.12.0 (fetched by CMake if not found)
- **Test-only dependencies** (fetched by CMake if needed):
    - [Catch2](https://github.com/catchorg/Catch2) v3
    - [nlohmann/json](https://github.com/nlohmann/json) v3.12.0

The Basis Set Exchange (BSE) JSON parser currently included should be considered an **example**, not a stable
or mandatory part of the core library.


## Examples and documentation

This project intentionally does **not** provide user-facing examples in the README.

- The library is in an early stage
- APIs may change
- The test suite serves as the primary form of documentation and usage examples

If you want to understand how the library is meant to be used, start with the tests.


## License & copyright notice

This library is licensed under MIT license.


## Data and attribution

<!--
### Physical constants and unit conversions

This library uses standard **physical constants and unit conversions** where needed (e.g., Bohr radius, Hartree to eV, atomic mass unit, etc.). The numerical values are taken from **CODATA recommended values** and other publicly available references.

- **CODATA 2018 (or latest available)** recommended values are used for:
    - Bohr radius (a₀)
    - Planck constant (h, ħ)
    - Elementary charge (e)
    - Atomic mass unit (u)
    - Electron mass (mₑ)
    - Other derived constants (Hartree energy, etc.)

- All constants are embedded in the library as `constexpr` values to ensure **compile-time safety, reproducibility, and unit correctness**.

- Users who require the **most up-to-date or experimentally precise constants** for high-accuracy calculations should consult the official **CODATA releases** (https://physics.nist.gov/cuu/Constants/) and adjust values accordingly.

- Values are intended for **internal computations, testing, and examples**. They are **not intended as a substitute for official physical reference tables** in precision quantum chemistry simulations.
-->

### Scientific data sources and attribution

This library includes **reference atomic and nuclear data** used to construct default `Element` and `Atom` properties 
(e.g. atomic mass, nuclear radius, and derived radii). These data are compiled from publicly available scientific 
sources and peer-reviewed literature, and are intended **for representation and testing purposes**, not as an 
authoritative physical database.

Specifically:

- **Atomic masses** are taken from publicly available reference tables, including:
    - PubChem periodic table (https://pubchem.ncbi.nlm.nih.gov/ptable/atomic-mass/)
    - Wikipedia (https://en.wikipedia.org/wiki/Neutron) (this library treat neutron as the zeroth element)

- **Nuclear charge radii** for elements with atomic number *Z < 110* are based on:
    - Visscher and Dyall, *Atomic Data and Nuclear Data Tables* **67**, 207 (1997), values reported in atomic units.

- For *Z ≥ 110*, nuclear radii are computed using the empirical relations:
    - Mass number estimation:  
      `A(Z) = 0.004467·Z² + 2.163·Z − 1.168`
    - Nuclear radius (in fm):  
      `r = 0.57 + 0.836·A^{1/3}`

  as described in:
    - D. Andrae, *Physics Reports* **336**, 414 (2000)
    - D. Andrae, *Nuclear charge density distributions in quantum chemistry*, in  
      *Relativistic Electronic Structure Theory, Part 1*, P. Schwerdtfeger (Ed.),  
      Theoretical and Computational Chemistry, Vol. 11, Elsevier (2002)

- **Atomic radii** are derived from published empirical formulas and fitted parameters, including:
    - DOI: 10.1016/j.theochem.2008.06.020
    - DOI: 10.3390/i3020087
    - DOI: 10.1002/qua.22202

In several cases, numerical values may be **rounded, reformatted, or recomputed from published formulas** for 
consistency, simplicity, or testing purposes. Such derived or modified values **should not be treated as primary 
reference data**.

Users requiring authoritative physical constants or high-precision nuclear models for production calculations should 
consult the **original publications** and **specialized databases**.

### Basis set data and attribution

Some tests and examples in this repository use basis set data obtained from the Basis Set Exchange (BSE), including 
well-known basis sets such as def2-SVP and ANO-R0. These data are used solely for testing and validation purposes, 
to exercise parsing, representation, and internal consistency checks.

In some cases, the basis set data may be truncated, reordered, or numerically modified (e.g. rounded coefficients or 
reduced shells) to keep the test cases small and readable. Such modified data must not be considered equivalent to the 
original published basis sets and are not suitable for production quantum chemistry calculations.

The original basis sets and their bibliographic references are curated and distributed by the Basis Set Exchange:

https://www.basissetexchange.org

Users who rely on the original, unmodified basis set data in scientific work should cite the appropriate original 
publications, as indicated by the Basis Set Exchange.


## Final note

This repository exists because *software structure matters*, even (and especially) in scientific code.
