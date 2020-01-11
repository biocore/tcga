#!/bin/bash -l

### TORQUE stuff here ####

### To send email when the job is completed:
#PBS -m ae
#PBS -M adswafford@ucsd.edu

# Tell Sun grid engine that this is a Bash script
#PBS -S /bin/bash
# Write errors to this file - make sure the directory exists
#PBS -e /home/adswafford/cluster/logs/$PBS_JOBID.elog

# Log output to this file - make sure the directory exists
#PBS -o /home/adswafford/cluster/logs/$PBS_JOBID.olog

#PBS -l walltime=120:00:00
#PBS -l nodes=1:ppn=30

# maximum amount of memory used by any single process
#PBS -l mem=250gb
# ALTERNATIVELY: -l mem=32gb # amount of physical memory used by the job

### Run in the queue named
### #PBS -q med4gb

# name of the job
#PBS -N prep_lung

# BASH stuff here
### Switch to the working directory; by default TORQUE launches processes
### from your home directory.
cd $PBS_O_WORKDIR
echo Working directory is $PBS_O_WORKDIR

# Calculate the number of processors allocated to this run.
NPROCS=`wc -l < $PBS_NODEFILE`

# Calculate the number of nodes allocated.
NNODES=`uniq $PBS_NODEFILE | wc -l`

### Display the job context
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Using ${NPROCS} processors across ${NNODES} nodes

source ~/.bashrc
set -e
cpus=$PBS_NUM_PPN

export TMPDIR=/panfs/panfs1.ucsd.edu/panscratch/$USER
tmp=$(mktemp -d --tmpdir)
export TMPDIR=$tmp
trap "rm -r $tmp; unset TMPDIR" EXIT

cd $tmp

aligner=${aligner:-bowtie2}
persample=/projects/tcga-data/human_filter/per_sample/combined/tcga_per_sample.list
db=/projects/genome_annotation/profiling/dbs/wol/shogun
bowtie2_env=bowtie2
shogun_env=shogun-1.0.6
utree_dir=/projects/genome_annotation/profiling/tools/shogun/UTree/2.0RF/linux64
burst_dir=/projects/genome_annotation/profiling/tools/shogun/BURST/0.99.7f
bowtie2_dir=/projects/genome_annotation/profiling/dbs/wol/shogun

export PATH=$utree_dir:$burst_dir:$bowtie2_dir:$PATH
declare -A a2ext=( [bowtie2]=sam [utree]=tsv [burst]=b6 )

#do shogun
echo 'starting shogun'
conda activate $shogun_env

while read tumor; do
    echo $tumor
    for file in /projects/tcga-data/$tumor/data/Bacterial*_Classified*.fasta
    do
    id=${file##*/}
    id=${id%.*}
    echo 'shogun-ing ' $id
    if [ ! -d /projects/tcga-data/shogun/gotus/$tumor ]
    then
        mkdir /projects/tcga-data/shogun/gotus/$tumor
    fi

    echo Directory is `pwd`
    lext=$aligner.${a2ext[$aligner]}
    if [ ! -f /projects/tcga-data/shogun/gotus/$tumor/$id.$lext.xz ]
       then
       shogun align -t $cpus -d $db -a $aligner -i $file -o .
       ls
       mv alignment.$lext $id.$lext
       shogun assign_taxonomy -d $db -a $aligner -i $id.$lext -o $id.profile.tsv
       for level in phylum genus species
       do
           shogun redistribute -d $db -l $level -i $id.profile.tsv -o $id.redist.$level.txt
       done
       shogun normalize -i $id.profile.tsv -o $id.profile.norm.tsv
       xz -T$cpus -9 $id.$lext
       mv $id.* /projects/tcga-data/shogun/gotus/$tumor/
    fi
    done
    echo $tumor 'done' >  /projects/tcga-data/shogun/gotus/$tumor/$tumor.done
done </projects/tcga-data/shogun/20191015_lung.list

conda deactivate
