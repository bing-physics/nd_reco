# NOTE:
This smearing is for beam monitoring only, for standard SAND/STT smearing, check out this one: https://github.com/bingguo1/FastReco

# nd_reco

* Input : Edep-Sim root file

* Output: ROOT tree file

# how to run

1. cmake .
2. make
3.
```
./smear_for_beamMonitoring -g 870 -v both -i inputEdepSimFile.root -o out.root
```
or check the detailed usage by
```
./smear_for_beamMonitoring -h
```
# how to read output

> root readEdep_smeared.C
