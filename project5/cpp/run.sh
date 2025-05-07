#!/bin/bash

#parameters
re=400
level=8
dtmax=0.00001
radius=0.078539816 # droplet radius (0.025*pi)
rho1=1.2           # oil density
rho2=1             # water density
mu1=2.52e-1        # oil dynamics viscosity
mu2=1.26e-2        # water dynamic viscosity
u0=5.0475
tend=1
sigma=1

./lid $level $tend $re $dtmax $radius $rho1 $rho2 $mu1 $mu2


