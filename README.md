# PlutoFiles
This github repo contains code to calculate the initial jet conditions of special relativistic hydrodynamic PLUTO simulations of Long Gamma Ray Bursts in the PLUTO_PY_FILES directory. This directory also contains code to process stellar envelope profiles from Woosley & Heger (2006) in a way that allows them to be easily imported into PLUTO 4.3 simulations as the initial condition.

The PLUTO_SIM_FILES directory contains the code that is used in the PLUTO simulation. The init.c file can be used to set the appropriate initial conditions while the TagCells.cpp file outlines the algorithm that allows PLUTO to adaptively refine the frid around the jet as the jet propagates through the domain.
