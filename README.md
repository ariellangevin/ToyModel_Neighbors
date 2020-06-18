# MATLAB Toy Model of Impact of Neighbors Antibiotic Export

This is a toy model describing how export of antibiotics by certain cells can produce harmful effects to their neighbors. These scripts simulate cell growth with different neighbors in the presence of antibiotics, using an agent-based model with Moore neighborhood architecture to represent the spatial interactions between cells and the environment. Each cell's growth and antibiotic concentration are represented by ordinary differential equations. Please see associated manuscript below for more details.

### Related Manuscript:

Wen, X., Langevin, A.M. & Dunlop, M.J. Antibiotic export by efflux pumps affects growth of neighboring bacteria. *Sci Rep* 8, 15120 (2018). https://doi.org/10.1038/s41598-018-33275-4.

## File Guide:

*neighborHard.m*: Hard coded script to enter inputs, run ODE function, as well as plot simulation results.

### Function:

*neighbor_func.m*: This updates the map of cells with their own biomass and antibioitic concentration over time, derived by ODEs. Mathematical equations for model shown in Methods section of manuscript.
