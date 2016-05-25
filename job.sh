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

###################
# Global settings #
###################
NAME="Yhuman_Ychimp"
FASTA_FILE=$HOME/chimp.fasta


BINARY=target/release/asgart
DATE=$(date +"%d_%m_%Y-%H_%M")
WORK_DIR="/tmpdir/franklin/${DATE}_${NAME}"
mkdir $WORK_DIR

cp $BINARY $WORK_DIR

echo "Working in $WORK_DIR"
cd $WORK_DIR

for kmer_size in 500 1500
do
	for max_gap_size in 1000 1500 2000
	do
		export RUST_BACKTRACE=1
		./asgart ~/human.fasta ~/chimp.fasta $kmer_size $max_gap_size
	done
done

#cargo run --release -- ~/chimp.fasta ~/chimp.fasta 500 1500 --reverse --translate
