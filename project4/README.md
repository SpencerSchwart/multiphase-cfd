# Instructions #
## C++ Program ##
For project 4, only a C++ program is available. To run the C++ program's solution to project 4, do
```bash
./build.sh
```
then
```bash
./run.sh
```
in project4/cpp. You can also specify which method to calculate interfacial normal and/or which method to solve the VOF advection by, for example,
```bash
YOUNGS=1 BASILISK_VOF=1 ./build.sh
```
where YOUNGS=1 will use Youngs-Finite-Difference method and BASILISK_VOF=1 corresponds to using the VOF advection approach as done by Basilisk (http://basilisk.fr/src/vof.h). By default, the code will use the Mixed-Youngs-Centered (MYC) method for normal calculation and the basic Eulerian Implicit (EI) method for VOF advection.

One can also specify parameters like Reynolds #, level of refinement, simulation duration, timestep size, and the droplet's radius in run.sh.
