# Actin rings Monte Carlo simulation package

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cumberworth.github.io/ActinRingsMC.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cumberworth.github.io/ActinRingsMC.jl/dev)
[![Build Status](https://github.com/cumberworth/ActinRingsMC.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/cumberworth/ActinRingsMC.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/cumberworth/ActinRingsMC.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/cumberworth/ActinRingsMC.jl)

A Julia simulation package with functions for running adaptive umbrella sampling Monte Carlo simulations of crosslinked filamentous (e.g. actin) rings.

[Repository](https://github.com/cumberworth/ActinRingsMC.jl)

[Documentation](http://www.cumberworth.org/ActinRingsMC.jl)

This package implements the simulation protocol outlined in the paper X.
The input and output files for the simulations in the paper that use this package can be found here.

## Installation

Installation instructions for Julia can be found at the [official website](https://github.com/cumberworth/ActinRingsMC.jl).

The package can be installed by starting the Julia REPL, typing `]` to enter package mode, and running
```
add ActinRingsMC
```
to install from the General x, or by running
```
add https://github.com/cumberworth/ActinRingsMC.jl
```
to install directly from the development repository.

## Running a simulation

In `examples` directory of the repository is a script for umbrella sampling.
It has both an initial run and a continuation run.
The commands can be run in the Julia REPL, or, when in the same directory as the script, it can be run with
```
julia run-umbrella-sampling.jl
```
To achieve good sampling, one should experiment with the number of steps per iteration and the number of iterations.

## Analysis and visualization

The simulations output two data file types per simulation, or in the case of an umbrella sampling run, per iteration.
The file type with an `.ops` extension contains order parameters saved at steps determined by the write interval and can be read in as a dataframe for analysis and plotting.
Currently the order parameters include the energy, lattice height, and ring radius.
The file type with a `.vtf` extension is able to be read by the molecular visualization program, [VMD](https://www.ks.uiuc.edu/Research/vmd/), such that the configurations of the simulations can be spatially visualized.
When viewing in VMD, it is recommended to use "Chain" for "Coloring Method" and "VDW" for "Drawing Method".
The simulations also output a file containing all the system and simulation parameters used in the run to a `.parms` file.

Umbrella sampling runs output additional file types.
The output of these file is based on every step, not only those determined by the write frequency.
The top row gives the lattice height (or bin, although binning has not been tested).
Each subsequent row contains data from an iteration.
`.counts` gives the number of times each lattice height was visited, `.freqs` gives a normalized version of this.
`.biases` includes the biases that were used to for that iteration (in contrast to those that would be calculated from its data).

A related python package, [actinrings](https://github.com/cumberworth/actinrings), includes code for analyzing and plotting the output from these simulations, including Landau free energies and ring constriction forces; see its documentation for details.

# Links

[Julia programming language]()

[Associated paper]()

[Replication package data]()

[Replication analysis package]()

[VMD]()
