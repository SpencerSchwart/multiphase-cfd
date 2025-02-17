# Instructions #
## Basilisk  Program ##
To run the Basilisk solution to project 1, simply do
```bash
make vortex.tst
```
in project1/basilisk/. One can edit the viscosity and wether to use the upwind or central method inside of vortex.c.

## C++ Program ##
To run the C++ program's solution to project 1, do
```bash
./build.sh
```
then
```bash
./vortex
```
in project1/cpp. You can also specify viscosity and/or which method for the advection term to use by
```bash
UPWIND=1 MU=0.5 ./build.sh
```
where UPWIND=0 corresponds to the central method.
