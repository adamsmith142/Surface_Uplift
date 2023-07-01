# SUFR (Surface Uplift From Relict topography)
Code to calculate surface uplift based on the method of Fox 2019. Topotoolbox is required.

This code uses topotoolbox to read in a DEM, extract a stream network, and perform an inversion of a relict river network to infer surface uplift and channel steepness index parameters across a grid.

The code writes matrices of the x grid, y grid, surface uplift and channel steepness that can be used to plot maps. 

Cite as:
Smith, A.G.G., Fox, M., Schwanghart, W., and Carter, A., 2022, Comparing methods for calculating channel steepness index: Earth-Science Reviews, v. 227, p. 103970, doi:10.1016/j.earscirev.2022.103970.

Fox, M., 2019, A linear inverse method to reconstruct paleo-topography: Geomorphology, v. 337, p. 151–164, doi:10.1016/j.geomorph.2019.03.034.


