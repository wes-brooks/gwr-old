#!/bin/bash
for i in {0..1799}
do
   Rscript code/simulation/simulations.r 0 $i
done