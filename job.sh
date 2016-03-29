#!/bin/sh
#SBATCH -J Palindromes
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --threads-per-core 1
#SBATCH --cpus-per-task=20
#SBATCH --cpu_bind=none
#--time=1000:00:00
#SBATCH --mail-user=franklin.delehelle@irit.fr
#SBATCH --mail-type=ALL

BINARY=target/release/palindromes
# Using symlink to DNA
Y_FILE=$HOME/Y.fasta

DATE=$(date +"%d_%m_%Y-%H_%M")
WORK_DIR="/tmpdir/franklin/$DATE"
mkdir $WORK_DIR

cp $BINARY $WORK_DIR
cp $Y_FILE $WORK_DIR

echo "Working in $WORK_DIR"
cd $WORK_DIR
/usr/bin/time -v ./palindromes
