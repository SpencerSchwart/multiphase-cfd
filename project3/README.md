# Instructions #
## Basilisk  Program ##
To run the Basilisk solution to project 3, simply do
```bash
make lid.tst
```
in project1/basilisk/. One can edit the viscosity and wether to use the upwind or central method inside of lid.c.

## C++ Program ##
To run the C++ program's solution to project 2, do
```bash
./build.sh
```
then
```bash
./lid
```
in project2/cpp. You can also specify viscosity and/or which method for the advection term to use by
```bash
UPWIND=1 RE=400 ./build.sh
```
where UPWIND=0 corresponds to the central method.
