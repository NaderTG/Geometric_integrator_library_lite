# Geometric_integrator_library_lite
This code is reduced version of a bigger library and the main reason is just testing and prototyping. This project is a C++ implementation of geometric numerical integrators, including stochastic integrators. The aim is to develop this code into a fully fledged parallel C++/MPI library.

There are two components to the project: ordinary differential equations and partial differential equations. The former is based on discrete mechanics and the latter is based on discrete multisymplectic geometry.

For finite dimensional systems, we plan on adding the following type of integrators: -Variational integrators (deterministic and stochastic), -and integrators based on generating functions (symplectic and Poisson).

For infinite dimenisonal systems, we rely heavily on fibre bundle, jet bundle and multisymplectic geometry. One advantage of considering jet bundles is that the mathematical abstraction lends itself to code abstraction. Restricting ourselves to fibre bundle and jet bundle description, we can use the code to numerically solve non-Hamiltonian PDEs, such as the heat equation or Poisson equation. We also provide an example implementation of multigrid method for elliptic pdes.

TODO LIST:

    Finish the geometry part
    Add abstract PDEtype class
    Add abstract class for Lagrangian and Hamiltonian

