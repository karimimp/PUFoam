# PUFoam

**PUFoam** is an [OpenFOAM](http://openfoam.org/) solver to simulate expansion
of polyurethane foams.

It includes **Quadrature Method of Moments** (QMOM) to solve a
**Population Balance Equation** (PBE) determining the bubble size distribution
inside PU foams.

Kinetics of the reactions including gelling, blowing and evaporation of
the physical blowing agents (n-pentane and R11) are incorporated into the solver.

Finally, the kinetics and PBE have been coupled to describe the time
evolution of the foaming/filling process.

In order to compile the solver and the required libraries execute the following:
```
user@machine> ./Allwmake
```


##### *Common compilation errors:*
The common compilation errors and the solution are listed here:

***Error:***

```
fatal error: sys/cdefs.h: No such file or directory
```

***Solution:***

```
sudo apt-get install libc6-dev-i386
sudo apt-get install libc6-dev
```

***Error:***

```
cannot find crti.o: No such file or directory
```

***Solution:***

```
LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LIBRARY_PATH
export LIBRARY_PATH
```

### References:
1. M. Karimi, H. Droghetti, D. L. Marchisio, `PUFoam`: a novel open-source CFD solver for the simulation of expanding and reacting polyurethane foams, Computer Physics Communications, Vol. 217, pp. 138--148, 2017.
2. M. Karimi, D. L. Marchisio, A Baseline Model for the Simulation of
Polyurethane Foams via the Population Balance Equation, Macromolecular Theory
and Simulations 24 (2015) 291–300.
3. M. Karimi, H. Droghetti, D. L. Marchisio, Multiscale Modeling of Expand-ing
Polyurethane Foams via Computational Fluid Dynamics and Popula-tion Balance
Equation, Macromolecular Symposia 360 (2016) 108–122.
4. P. Ferkl, M. Karimi, D. L. Marchisio, J. Kosek, Multi-scale modelling of
expanding polyurethane foams: coupling macro- and bubble-scales, Chemical
Engineering Science (2016).
5. MoDeNa-EUProject, Modelling of morphology development of micro- and
nanostructures, 2015. URL: https://github.com/MoDeNa-EUProject/MoDeNa
6. D. L. Marchisio, R. O. Fox, Computational Models for Polydisperse Particulate
and Multiphase Systems, Cambridge University Press, 2013.