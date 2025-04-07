# Instructions #

Only a C++ program (no Basilisk program) for this project is made.

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
UPWIND=1 RE=1000 ./build.sh
```
where UPWIND=0 corresponds to the central method. UPWIND=0 and RE=400 by default.
