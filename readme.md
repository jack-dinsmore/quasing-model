# Quising model

This repo contains a Rust implementation of the [Wolff algorithm](https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.62.361) in order to find phase transition temperatures and possibly critical exponents of the Ising model for quasicrystalline lattices.

It is Jack Dinsmore's final project for Stanford's Condensed Matter Theory course in 2023.

It also contains a mostly finished implementation of the [continuous Wolff algorithm](https://journals.aps.org/pre/pdf/10.1103/PhysRevE.66.066110) for the transverse Ising model.

## Directories

* **figs/** contains the figures for the paper
* **python/** contains code to generate the figures and do basic analysis
* **src/** contains the Rust MCMC code
* **paper/** contains the final report.