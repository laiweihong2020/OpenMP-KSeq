#!/bin/bash
#SBATCH --partition=physical
#SBATCH --constraint=physg5
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --job-name=diag

module load gcc/10.1.0
g++ -Wall -fopenmp -o kseqalignDiag.out kseqalignDiag.cpp
cat mseq-big13-example.dat | ./kseqalignDiag.out