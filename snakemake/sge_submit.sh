#!/bin/bash
#Submit to the cluster, give it a unique name
#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=1.9G,h_rt=20:00:00,tmem=1.9G
#$ -pe smp 2

# join stdout and stderr output
#$ -j y
#$ -R y


snakemake -s bullseye.smk \
--jobscript sge_qsub.sh \
--cluster-config sge_cluster.yaml \
--cluster-sync "qsub -l tmem={cluster.tmem},h_vmem={cluster.h_vmem},h_rt={cluster.h_rt}" \
-j 5 \
--nolock \
--rerun-triggers mtime \
--rerun-incomplete \
--latency-wait 100 
