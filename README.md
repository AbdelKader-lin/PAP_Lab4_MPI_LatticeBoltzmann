# Parallel Lattice Boltzmann Fluid Simulation with MPI

A C/MPI project exploring distributed-memory parallelization of the **Lattice Boltzmann Method (LBM)** for fluid simulation. The implementation progresses through several exercises and communication strategies, with tools for correctness checking, visualization and performance benchmarking.

## What this project demonstrates

- distributed-memory parallel programming with **MPI**
- domain decomposition for a 2D numerical simulation
- communication between neighboring MPI processes
- handling of ghost cells / halo regions
- correctness testing of communication patterns
- performance benchmarking and weak-scaling experiments
- scientific computing in C

## Project structure

The repository contains a sequence of implementations from `exercise_0.c` onward. These exercises progressively introduce and refine the parallel MPI version of the simulation.

Supporting tools include:

- `check_comm` - isolates communication patterns to make halo/ghost-cell exchanges easier to inspect
- `config.txt` - simulation configuration
- `bench_correction.sh` - benchmarking support
- `gen_animate_gif.sh` - generates an animated visualization of simulation output using gnuplot
- `Makefile` - compilation and build configuration

## Requirements

You need:

- a C compiler such as GCC
- an MPI implementation such as OpenMPI or MPICH
- `make`
- optionally `gnuplot` for rendering simulation results

## Build

```bash
make
```

## Run

Run the simulation with MPI and select an exercise using `--exercise` or `-e`:

```bash
mpirun -np 8 ./lbm --exercise 1
```

Equivalent short form:

```bash
mpirun -np 8 ./lbm -e 1
```

Exercise 0 is the sequential baseline and should be run with a single process.

A different configuration file can be supplied as an argument:

```bash
mpirun -np 8 ./lbm -e 1 config-other.txt
```

## Communication debugging

The `check_comm` executable reproduces the MPI communication independently of the complete simulation, which makes it easier to inspect distributed subdomains and ghost cells.

```bash
mpirun -np 8 ./check_comm -e 1
```

## Benchmarking

Output can be disabled when measuring execution time:

```bash
mpirun -np 8 ./lbm -e 1 --no-out
```

The weak-scaling option increases the mesh size with the number of processes:

```bash
mpirun -np 8 ./lbm -e 6 --scale 8 --no-out
```

This allows the implementation to be evaluated as the amount of parallel work grows with the available MPI processes.

## Visualization

Simulation output can be converted to an animated GIF using the included gnuplot-based script:

```bash
./gen_animate_gif.sh resultat.raw output.gif
```

## Context

Parallel programming / high-performance computing lab focused on applying MPI to a numerical fluid simulation.

## License

The original teaching code is distributed under the BSD license.
