#!/bin/bash

foamCleanTutorials

blockMesh 
checkMesh
mirrorMesh -overwrite -dict mirrorMeshDict_2
mirrorMesh -overwrite
renumberMesh -overwrite -noFunctionObjects



decomposePar
mpirun -np 4 reoFoam.C -parallel | tee log.solver

