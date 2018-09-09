# Experimental tools

This repository contains experimental tools for GRChombo (see below).
These tools are "experimental" for various reasons: some of them are not tested as
rigourously as the main code, some are tested but not maintained and may have to be
updated to work with the latest version of GRChombo, etc. None of the code in this
repository has been through code review. Use it at your own risk and check that it
does what you want it to do.

Currently the repository contains:

* A multigrid solver for the initial conditions

  * For the moment only the Hamiltonian Constraint is solved, but the intention is to make this more general. 
    See the [wiki](https://github.com/GRChombo/GRChombo/wiki) for further details of the current setup.

* Files for intermediate checkpointing, i.e. writing checkpoints even though the
  levels are not synchronised in time. This is useful for very deep mesh hierarchies.

* A tool for comparing two checkpoint files. This is a great way of testing
  GRChombo (it can even pick up errors in the checkpointing process if one stops and
  restarts before comparing checkpoint files).

# GRChombo
GRChombo is a new open-source code for numerical general relativity simulations. 
It is developed and maintained by a collaboration of numerical relativists with a 
wide range of research interests, from early universe cosmology to astrophysics 
and mathematical general relativity, and has been used in many papers since its
first release in 2015.

GRChombo is written entirely in C++14, using hybrid MPI/OpenMP parallelism and 
vector intrinsics to achieve good performance on the latest architectures.
Furthermore, it makes use of the Chombo library for adaptive mesh refinement
to allow automatic increasing and decreasing of the grid resolution in regions
of arbitrary shape and topology.

Please visit www.grchombo.org for the full list of developers and their
institutions, a list of publications using GRChombo, and some videos.

## Getting started
Detailed installation instructions and usage examples are available in
our [wiki](https://github.com/GRChombo/GRChombo/wiki).

## Contributing
We welcome feedback, bug reports, and contributions. Please consult the [wiki](https://github.com/GRChombo/GRChombo/wiki)
for our coding style and testing policy before filing a pull request.

## License
GRChombo is licensed under the BSD 3-Clause License. Please see LICENSE for details.
