#!/bin/sh
#SBATCH -J Palindromes
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --threads-per-core 1
#SBATCH --cpus-per-task=20
#SBATCH --cpu_bind=none
#SBATCH --mail-user=franklin.delehelle@irit.fr
#SBATCH --mail-type=ALL

BINARY=target/release/palindromes
# Using a symbolic link to the DNA
Y_FILE=$HOME/Y.fasta

WORK_DIR=$(mktemp -d -p /tmpdir/franklin/)

cp $BINARY $WORK_DIR
cp $Y_FILE $WORK_DIR

echo "Working in $WORK_DIR"
cd $WORK_DIR
./palindromes
