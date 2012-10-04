#i /bin/sh

python generate_seeds.py
condor_submit beauty.condor

