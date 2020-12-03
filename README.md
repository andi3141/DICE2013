# DICE2013
This is the effort of a simple, straight forward implementation of the DICE model in Julia. The implementation is based on the [DICE manual](http://www.econ.yale.edu/~nordhaus/homepage/homepage/documents/DICE_Manual_100413r1.pdf), which unluckily is incomplete, and this [publication](https://arxiv.org/pdf/1812.08066.pdf).

The package uses JuMP, Ipopt, and MUMPS for optimization, and Plots for plotting.

There are three basis scenarios, "baseline", "optimal", and "2degree" available, which can be called via

`runScenario(scenario)`.

