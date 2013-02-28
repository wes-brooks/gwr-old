#! /bin/sh
for i in 410 411 412 415 419 420 421 422
do
Rscript code/simulation/simulations-i386.r 43 $i
done