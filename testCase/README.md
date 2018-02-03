# How to run the test case:
This test case demonstrates a 2-dimensional example of the foam expansion inside
 a cup. User shall first provide the required input variables through different
 dictionaries in `constant` directory. It includes the kinetics, rheology,
 simulation mode and PBE properties for the PU foaming process. Further, the
 simulation characteristics such as time step, discretization schemes and etc
 should be defined in `system` directory. Finally to run the simulation the
 `run.sh` script can be executed which sets the boundary and initial conditions,
 run the solver and save the output in a `log` file.