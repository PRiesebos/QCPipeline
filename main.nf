"""
File:         main.nf
Created:      2024/04/11
Last Changed: 2024/04/25
Author:       Peter Riesebos & Orfeas Gkourlias
"""

nextflow.enable.dsl=2

process concatCHRFiles {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path vcfFilesPath

    output:
    path "merged_output.vcf.gz"
    path "sorted_merged_output.vcf.gz", emit: sortedVCF

    script:
    """
    mkdir -p ${params.inputDir}/qc_logs
    mkdir -p ${params.inputDir}/figures
    bcftools concat ${vcfFilesPath}/chr*.vcf.gz -Oz -o merged_output.vcf.gz
    bcftools sort merged_output.vcf.gz -Oz -o sorted_merged_output.vcf.gz
    """
}

process splitMultiAllelicVariants {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path vcfFilesPath

    output:
    path "split_sorted_merged_output.vcf.gz"

    script:
    """
    bcftools norm -m -any ${vcfFilesPath} -Oz -o split_sorted_merged_output.vcf.gz
    """
}

// process filterMultiAllelicVariants {
//     errorStrategy 'retry'
//     maxRetries 1

//     time '4h'
//     memory '8 GB'
//     cpus 1

//     input:
//     path vcfFilesPath

//     output:
//     path "no_multi_allelic.vcf.gz"

//     script:
//     """
//     bcftools view --max-alleles 2 ${vcfFilesPath} -Oz -o no_multi_allelic.vcf.gz
//     """
// }

process fixGTAnnot {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path vcfFilesPath

    output:
    path "no_multi_allelic_dotrevive.vcf.gz", emit: unfilteredVCF

    script:
    """
    # 1. fix GT field annotation
    python3 ${projectDir}/bin/dotrevive.py -i ${vcfFilesPath} -o no_multi_allelic_dotrevive.vcf.gz
    """
}


process filterVariants {
    storeDir "${params.inputDir}"
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path vcfFile

    output:
    path "no_multi_allelic_dotrevive-filtered.vcf.gz", emit: filteredVCF
    path "no_multi_allelic_dotrevive-filtered.log.gz"

    script:
    """
    # 1. Define command arguments
    commandArguments="--input ${vcfFile} \
    --output no_multi_allelic \
    --call_rate 0.5 \
    --filtered_depth 5 \
    --genotype_quality 10 \
    --minor_allele_frequency 0.01 \
    --no_snv_vqsr_check \
    --no_indel_vqsr_check \
    --remove_non_pass_snv \
    --remove_non_pass_indel \
    --replace_poor_quality_genotypes \
    --output ${params.inputDir}/no_multi_allelic_dotrevive"

    # 2. Run the custom VCF filter script
    python3 ${projectDir}/bin/custom_vcf_filter.py \${commandArguments} \
    | tee ${params.inputDir}/qc_logs/custom_vcf_filter.log
    """
}

process getMetrics {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path unfilteredVCF
    path filteredVCF

    // output:
    // path "qc_logs/variant_count.txt", emit: variantCount

    script:
    """
    zcat ${unfilteredVCF} | grep -v '^#' | wc -l | xargs -I {} echo -e "variant count:\t{}" > ${params.inputDir}/qc_logs/variant_count.txt
    zcat ${filteredVCF} | grep -v '^#' | wc -l | xargs -I {} echo -e "variant count filtered:\t{}" >> ${params.inputDir}/qc_logs/variant_count.txt
    """
}

process convertToPlinkFormat {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path vcfFile

    output:
    path "*.bed", emit: bedFile
    path "*.bim", emit: bimFile
    path "*.fam", emit: famFile

    script:
    """
    plink2 --vcf ${vcfFile} --make-bed
    """
}

process calculateMissingness {
    storeDir "${params.inputDir}"
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path bedFile
    path bimFile
    path famFile
    path vcfFile

    output:
    path "PRE_FILTER.smiss"
    path "POST_FILTER.smiss", emit: missingFile

    script:
    """
    plink2 --vcf ${vcfFile} --missing --out PRE_FILTER
    plink2 --bfile plink2 --missing --out POST_FILTER
    """
}

// if less than 50 samples use --bad-freqs. Needs an alternative solution?
process createHetFile {
    storeDir "${params.inputDir}"
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path bedFile
    path bimFile
    path famFile

    output:
    path "*.het"

    script:
    """
    plink2 --bfile plink2 --het cols=hom,het,nobs,f --bad-freqs --out heterozygosity_output
    """
}

// Create a txt file with samples where missingness => 50%
process findMissingSamples {
    storeDir "${params.inputDir}"
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path smiss

    output:
    path 'qc_logs/filtered_samples.txt', emit: filteredSamplesFile

    script:
    """
    python ${projectDir}/bin/filterMissingness.py ${smiss} ${params.inputDir}/qc_logs/filtered_samples.txt --threshold 0.5
    """
}

// Filter these samples. Output -> filtered plink files
process filterMissingSamples {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path filteredSamplesFile
    path bedFile
    path bimFile
    path famFile

    output:
    path 'data_keep.bed', emit: bedFile
    path 'data_keep.bim', emit: bimFile
    path 'data_keep.fam', emit: famFile

    script:
    """
    plink2 --bfile plink2 --keep ${filteredSamplesFile} --make-bed --out data_keep
    """
}

process findHetSamples {
    storeDir "${params.inputDir}"
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path het_file

    output:
    path "figures/*.png"
    path "qc_logs/*Failed.txt"
    path "qc_logs/*FailedSamplesOnly.txt", emit: failedHetSamples

    script:
    """
    python3 ${projectDir}/bin/heterozygosityCheck.py ${het_file} ${params.inputDir}
    """
}

process filterHetSamples {
    storeDir "${params.inputDir}"
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path failedHetSamples
    path bedFile
    path bimFile
    path famFile

    output:
    path 'output.bed', emit: bedFile
    path 'output.bim', emit: bimFile
    path 'output.fam', emit: famFile

    script:
    """
    plink2 --bfile data_keep \
       --remove ${failedHetSamples} \
       --make-bed \
       --out output
    """
}

process createMetricsFile {
    storeDir "${params.inputDir}"
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path bedFile
    path bimFile
    path famFile

    output:
    path "qc_logs/metrics_matrix.tsv", emit: metrics_matrix

    script:
    """
    python3 ${projectDir}/bin/CombineQCFiles.py ${params.inputDir} ${params.inputDir}/qc_logs/metrics_matrix.tsv
    """
}


process filterLowAltFreq {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 4

    input:
    path bed
    path bim
    path fam

    output:
    path "${bed.SimpleName}.over10.bed"
    path "${bed.SimpleName}.over10.bim"
    path "${bed.SimpleName}.over10.fam"
    
    script:
    """
    plink2 --sample-counts 'cols=homref,het,homalt,missing' --bfile ${bed.SimpleName} --out ${bed.SimpleName}
    nonalt.py -i ${bed.SimpleName}.scount -o ${bed.SimpleName}-keep.txt
    plink2 --keep ${bed.SimpleName}-keep.txt --bfile ${bed.SimpleName} --make-bed --out ${bed.SimpleName}.over10
    """
}

process filterRelated {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 4

    input:
    path bed
    path bim
    path fam

    output:
    path 'RelatednessCheck.bed', emit: bedFile
    path 'RelatednessCheck.bim', emit: bimFile
    path 'RelatednessCheck.fam', emit: famFile
    path 'related.kin0'
    path 'RelatednessPassedSamples.txt'

    script:
    """
    # 1. Do relatedness check
    plink2 --bfile ${bed.SimpleName}.over10 \
        --make-king-table \
        --king-table-filter ${params.kingTableFilter} \
        --out related \
        --threads 4 \

    # 2. Create sample list of non-related samples
    find_related_samples.R --kin_file related.kin0 --target_bed ${bed}

    # 3. Remove samples that are not on the list created above
    plink2 --bed ${bed} \
        --bim ${bim} \
        --fam ${fam} \
        --keep RelatednessPassedSamples.txt \
        --make-bed \
        --out RelatednessCheck
    """
}


process popProject {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '24 GB'
    cpus 4

    input:
    path bed
    path bim
    path fam
    
    output:
    path '*.png'
    path '*.pdf'
    path '1000G_PC_projections.txt'
    path 'PopAssignResults.txt'
    
    script:
    """
    project_samples_to_superpop.R --ref_bed ${params.refPath} \
        --target_bed ${bed} --ref_pop ${params.refPop}
    """
}

workflow {
    concatCHRFiles(params.inputDir)
    splitMultiAllelicVariants(concatCHRFiles.output.sortedVCF)
    // filterMultiAllelicVariants(splitMultiAllelicVariants.output)
    fixGTAnnot(splitMultiAllelicVariants.output)
    filterVariants(fixGTAnnot.output)
    getMetrics(fixGTAnnot.output.unfilteredVCF, filterVariants.output.filteredVCF)
    convertToPlinkFormat(filterVariants.output.filteredVCF)
    calculateMissingness(convertToPlinkFormat.output, fixGTAnnot.output)
    createHetFile(convertToPlinkFormat.output)
    findMissingSamples(calculateMissingness.output.missingFile)
    filterMissingSamples(findMissingSamples.output.filteredSamplesFile, convertToPlinkFormat.output)
    findHetSamples(createHetFile.output)
    filterHetSamples(findHetSamples.output.failedHetSamples, filterMissingSamples.output)
    createMetricsFile(filterHetSamples.output)
    filterLowAltFreq(filterHetSamples.output)
    filterRelated(filterLowAltFreq.output)
    popProject(filterRelated.output.bedFile, filterRelated.output.bimFile, filterRelated.output.famFile)
}