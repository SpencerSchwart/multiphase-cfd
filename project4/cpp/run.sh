#!/bin/bash

#parameters
re=400
level=8
tend=10
dtmax=0.0005
radius=0.078539816 # droplet radius (0.025*pi)

./lid $level $tend $re $dtmax $radius


