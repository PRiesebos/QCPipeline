# Quality Control pipeline for genotype data using Nextflow.

Pipeline made in colleboration with Orfeas Gkourlias for my UMCG genetics department internship.

## This pipeline uses the following nextflow.config:
```
nextflow.enable.dsl=2

process.executor = "slurm"
process.container = "" // path to .sif container from: apptainer pull docker://ogkourlias/pub-rna-nocache:latest

// Adjust the params to your usecase.
params {
    inputDir = "" // path to the root folder of the output from the pub-rna pipeline. Link: https://github.com/ogkourlias/pub-rna
    kingTableFilter = 0.04419417382
    refPath = '' // path to the plink reference .bed file from 1000 genomes build 38
    refPop = '' // path to reference file sample population list. Link: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
    refAFs = '' // path to plink .afreq file of 1000 genomes
    annotationGTF = '' // path to gencode.v44.primary_assembly.annotation.gtf
    beagleJarDir = '' // path to version 4.1 of beagle.jar. Link: https://faculty.washington.edu/browning/beagle/b4_1.html
    mapFile = '' // path to genetic map file: https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/. Build 38, concat all chromosome files into one.
}

singularity {
    enabled = true
    autoMounts = true
    runOptions = "--bind" // bind to a root folder dir e.g. /scratch/
    cacheDir = "" // path to cache directory
}
```
