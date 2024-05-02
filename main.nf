"""
File:         main.nf
Created:      2024/04/11
Last Changed: 2024/05/02
Author:       Peter Riesebos & Orfeas Gkourlias
"""

nextflow.enable.dsl=2

process concatCHRFiles {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '4 GB'
    cpus 1

    input:
    path vcfsPath

    output:
    path "sorted_merged_output.vcf.gz", emit: sortedVcf

    script:
    """
    mkdir -p ${params.inputDir}/figures
    mkdir -p ${params.inputDir}/qc_logs

    # hardcoded concat step so sorting isn't needed
    bcftools concat \
    ${vcfsPath}/chr1.vcf.gz ${vcfsPath}/chr2.vcf.gz ${vcfsPath}/chr3.vcf.gz ${vcfsPath}/chr4.vcf.gz \
    ${vcfsPath}/chr5.vcf.gz ${vcfsPath}/chr6.vcf.gz ${vcfsPath}/chr7.vcf.gz ${vcfsPath}/chr8.vcf.gz \
    ${vcfsPath}/chr9.vcf.gz ${vcfsPath}/chr10.vcf.gz ${vcfsPath}/chr11.vcf.gz ${vcfsPath}/chr12.vcf.gz \
    ${vcfsPath}/chr13.vcf.gz ${vcfsPath}/chr14.vcf.gz ${vcfsPath}/chr15.vcf.gz ${vcfsPath}/chr16.vcf.gz \
    ${vcfsPath}/chr17.vcf.gz ${vcfsPath}/chr18.vcf.gz ${vcfsPath}/chr19.vcf.gz ${vcfsPath}/chr20.vcf.gz \
    ${vcfsPath}/chr21.vcf.gz ${vcfsPath}/chr22.vcf.gz \
    -Oz -o sorted_merged_output.vcf.gz

    # Original steps: 
    # bcftools concat ${vcfsPath}/chr*.vcf.gz -Oz -o ${vcfsPath.SimpleName}.concat.vcf.gz
    # bcftools sort ${vcfsPath.SimpleName}.concat.vcf.gz -Oz -o ${vcfsPath.SimpleName}.sorted.concat.vcf.gz
    """
}


process fixGTAnnot {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path vcf

    output:
    path "${vcf.SimpleName}.gtfixed.sorted.concat.vcf.gz", emit: gtFixedVcf

    script:
    """
    # 1. fix GT field annotation
     dotrevive.py -i ${vcf} -o ${vcf.SimpleName}.gtfixed.sorted.concat.vcf.gz
    """
}


process splitMultiAllelicVariants {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path vcf

    output:
    path "${vcf.SimpleName}.split.gtfixed.sorted.concat.vcf.gz", emit: splitVcf

    script:
    """
    bcftools norm -m -any ${vcf} -Oz -o ${vcf.SimpleName}.split.gtfixed.sorted.concat.vcf.gz
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


process createGraphs {
    storeDir "${params.inputDir}"
    errorStrategy 'retry'
    maxRetries 0

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path vcf

    output:
    path "*stats.tsv.gz"
    path "*regions.tsv.gz"
    path "*rates.png"
    path "*freqs.png"

    publishDir "${params.inputDir}/figures", mode: 'move', pattern: '*.png'

    script:
    """
    # 1 run vcf_stats.py
    vcf_stats.py -i ${vcf} -o output_stats.tsv.gz
    
    # 2 run regions.py
    get_region.py -all -i output_stats.tsv.gz -g ${params.annotationGTF} -o output_stats_regions.tsv.gz

    # 3 run graph_summary.
    graph_summary.R output_stats_regions.tsv.gz
    """
}


process filterVariants {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path vcf

    output:
    path "${vcf.SimpleName}.filtered.vcf.gz", emit: filteredVcf
    path "${vcf.SimpleName}.filtered.log.gz", emit: filteredLog

    script:
    """
    # 1. Define command arguments
    commandArguments="--input ${vcf} \
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
    --output ${vcf.SimpleName}"

    # 2. Run the custom VCF filter script
    custom_vcf_filter.py \${commandArguments} \
    | tee ${params.inputDir}/qc_logs/custom_vcf_filter.log

    mv ${vcf.SimpleName}-filtered.vcf.gz ${vcf.SimpleName}.filtered.vcf.gz
    mv ${vcf.SimpleName}-filtered.log.gz ${vcf.SimpleName}.filtered.log.gz
    """
}


process getVariantCount {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path unfilteredVcf
    path filteredVcf

    output:
    path "variant_count.txt", emit: variantCount

    publishDir "${params.inputDir}/qc_logs", mode: 'move', pattern: '*.txt'

    script:
    """
    zcat ${unfilteredVcf} | grep -v '^#' | wc -l | xargs -I {} echo -e "variant count:\t{}" > variant_count.txt
    zcat ${filteredVcf} | grep -v '^#' | wc -l | xargs -I {} echo -e "variant count filtered:\t{}" >> variant_count.txt
    """
}


process getSampleCount {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path vcf
    path fam

    output:
    path "sample_count.txt", emit: sampleCount

    publishDir "${params.inputDir}/qc_logs", mode: 'move', pattern: '*.txt'

    script:
    """
    zcat ${vcf} | grep -m1 "^#CHROM" | cut -f 10- | tr "\t" "\n"  | wc -l | xargs -I {} echo -e "Sample count:\t{}" > sample_count.txt
    less ${fam} | grep -v '^#' | wc -l | xargs -I {} echo -e "Sample count filtered:\t{}" >> sample_count.txt
    """
}


process convertToPlinkFormat {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path vcf

    output:
    path "*.bed", emit: bed
    path "*.bim", emit: bim
    path "*.fam", emit: fam

    script:
    """
    plink2 --vcf ${vcf} --make-bed
    """
}


process convertToPlinkFormatAlt {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path vcf

    output:
    path "*.bed", emit: bed
    path "*.bim", emit: bim
    path "*.fam", emit: fam

    script:
    """
    plink2 --vcf ${vcf} --make-bed --out beagle_imputed
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
    path bed
    path bim
    path fam
    path vcf

    output:
    path "${vcf.SimpleName}.preFilter.smiss"
    path "${vcf.SimpleName}.postFilter.smiss", emit: missing

    script:
    """
    plink2 --vcf ${vcf} --missing --out ${vcf.SimpleName}.preFilter
    plink2 --bfile ${bed.SimpleName} --missing --out ${vcf.SimpleName}.postFilter
    """
}


process createHetFile {
    storeDir "${params.inputDir}"
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path bed
    path bim
    path fam

    output:
    path "*.het"

    script:
    """
    plink2 --bfile ${bed.SimpleName} --het cols=hom,het,nobs,f --bad-freqs --out heterozygosity_output
    """
}


process findMissingSamples {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path smiss

    output:
    path "${smiss.SimpleName}-filtered-samples.txt", emit: filteredSamples

    publishDir "${params.inputDir}/qc_logs", mode: 'copy', pattern: '*.txt'

    script:
    """
    filterMissingness.py ${smiss} ${smiss.SimpleName}-filtered-samples.txt --threshold 0.5
    # sleep 10s
    """
}


process filterMissingSamples {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path filteredSamples
    path bed
    path bim
    path fam

    output:
    path "${bed.SimpleName}.keep.bed", emit: bed
    path "${bed.SimpleName}.keep.bim", emit: bim
    path "${bed.SimpleName}.keep.fam", emit: fam

    script:
    """
    plink2 --bfile ${bed.SimpleName} --keep ${filteredSamples} --make-bed --out ${bed.SimpleName}.keep
    """
}


process findHetSamples {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path het_file

    output:
    path "*.png"
    path "*Failed.txt"
    path "*FailedSamplesOnly.txt", emit: failedHetSamples

    publishDir "${params.inputDir}/figures", mode: 'move', pattern: '*.png'
    publishDir "${params.inputDir}/qc_logs", mode: 'copy', pattern: '*.txt'

    script:
    """
    heterozygosityCheck.py ${het_file} ./
    """
}


process filterHetSamples {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path failedHetSamples
    path bed
    path bim
    path fam

    output:
    path "${bed.SimpleName}.bed", emit: bed
    path "${bed.SimpleName}.bim", emit: bim
    path "${bed.SimpleName}.fam", emit: fam

    script:
    """
    plink2 --bfile ${bed.SimpleName}.keep \
       --remove ${failedHetSamples} \
       --make-bed \
       --out ${bed.SimpleName}
    """
}


process createMetricsFile {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path bed
    path bim
    path fam

    output:
    path "metrics_matrix.tsv", emit: metrics_matrix

    publishDir "${params.inputDir}/qc_logs", mode: 'move', pattern: '*.tsv'

    script:
    """
    CombineQCFiles.py ${params.inputDir} metrics_matrix.tsv
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
    path "${bed.SimpleName}.over10.bed", emit: lowAltFreqBedFile
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
    path 'RelatednessCheck.bed', emit: bed
    path 'RelatednessCheck.bim', emit: bim
    path 'RelatednessCheck.fam', emit: fam
    path 'related.kin0'
    path 'RelatednessPassedSamples.txt'

    publishDir "${params.inputDir}/qc_logs", mode: 'move', pattern: '*.txt'

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

process convertBackToVCF {
    errorStrategy 'retry'
    maxRetries 1

    time '6h'
    memory '8 GB'
    cpus 1

    input:
    path bed
    path bim
    path fam

    output:
    path "final.vcf.gz", emit: finalVcf

    script:
    """
    plink2 --bed ${bed} --bim ${bim} --fam ${fam} --export vcf-4.2 bgz --out final
    """
}

process runBeagle {
    errorStrategy 'retry'
    maxRetries 1

    time '8h'
    memory '12 GB'
    cpus 4

    input:
    path vcfFile

    output:
    path "beagle_imputed.vcf.gz", emit: beagleVcf
    path "beagle_imputed.log"

    publishDir "${params.inputDir}/qc_logs", mode: 'copy', pattern: '*.log'
    publishDir "${params.inputDir}", mode: 'copy', pattern: '*.gz'


    script:
    """
    java -Xmx20g -jar ${params.beagleJarDir} \
    gtgl=${vcfFile} \
    out=beagle_imputed \
    map=${params.mapFile} \
    gprobs=true
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
    path "*.png"
    path "*.pdf"
    path "1000G_PC_projections.txt"
    path "PopAssignResults.txt"
    
    publishDir "${params.inputDir}/figures", mode: 'move', pattern: '*.{pdf,png}'
    publishDir "${params.inputDir}/qc_logs", mode: 'move', pattern: '*.txt'
    
    script:
    """
    project_samples_to_superpop.R --ref_bed ${params.refPath} \
        --target_bed ${bed} --ref_pop ${params.refPop}
    """
}


process popProjectAfterBeagle {
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
    path "*.png"
    path "*.pdf"
    path "1000G_PC_projections.txt"
    path "PopAssignResults.txt"
    
    publishDir "${params.inputDir}/after_beagle_pop", mode: 'move', pattern: '*.{pdf,png}'
    publishDir "${params.inputDir}/after_beagle_pop", mode: 'move', pattern: '*.txt'
    
    script:
    """
    mkdir -p ${params.inputDir}/after_beagle_pop
    project_samples_to_superpop.R --ref_bed ${params.refPath} \
        --target_bed ${bed} --ref_pop ${params.refPop}
    """
}


process targetPCA {
    storeDir "${params.inputDir}/final"

    time '6h'
    memory '16 GB'
    cpus 1

    input:
    path bed
    path bim
    path fam

    output:
    path "${bed.SimpleName}.toImputation.bed", emit: bed
    path "${bed.SimpleName}.toImputation.bim", emit: bim
    path "${bed.SimpleName}.toImputation.fam", emit: fam
    path "${bed.SimpleName}.SamplesToInclude.txt", emit: sampleFile

    publishDir "${params.inputDir}/qc_logs", mode: 'copy', pattern: '*.txt'

    script:
    """
    mkdir -p ${params.inputDir}/final
    # 1. Do PCA and find outliers
    target_pca.R --target_bed ${bed} --outlier_threshold ${params.populationOutlierThreshold} --out ${bed.SimpleName}
    # 2. Remove outlier samples
    plink2 --bed ${bed} \
        --bim ${bim} \
        --fam ${fam} \
        --output-chr 26 \
        --keep ${bed.SimpleName}.SamplesToInclude.txt \
        --make-bed \
        --threads 4 \
        --out ${bed.SimpleName}.toImputation
    """
}


process finalSNPandGenotypeQC {
    time '6h'
    memory '8 GB'
    cpus 1

    input:
    path bed
    path bim
    path fam

    output:
    path "${bed.SimpleName}.bed"
    path "${bed.SimpleName}.bim"
    path "${bed.SimpleName}.fam"

    publishDir "${params.inputDir}/final_snp_qc/", mode: 'move'

    script:
    """
    mkdir -p ${params.inputDir}/final_snp_qc/
    plink2 --bed ${bed} \
        --bim ${bim} \
        --fam ${fam} \
        --maf ${params.maf} \
        --geno ${params.geno} \
        --mind ${params.mind} \
        --hwe ${params.hwe} \
        --autosome \
        --make-bed \
        --out chrAll_QC \
        --output-chr 26 \
        --not-chr 0 25-26 \
        --set-all-var-ids @:#[b38]\\\$r,\\\$a \
        --new-id-max-allele-len 10 truncate \
        --threads 4
    """
}

process finalPCA {
    time '6h'
    memory '16 GB'
    cpus 1

    input:
    path bed
    path bim
    path fam

    output:
    path '*.png'
    path '*.pdf'
    path '*.txt'

    publishDir "${params.inputDir}/qc_logs", mode: 'move', pattern: '*.txt'
    publishDir "${params.inputDir}/figures", mode: 'move', pattern: '*.{pdf,png}'

    script:
    """
    mkdir -p ${params.inputDir}/${bed.SimpleName}/final/
    final_pca.R --target_bed ${bed} 
    """
}


process shuffleSampleOrder {
    time '6h'
    memory '1 GB'
    cpus 1

    input:
    path bed
    path bim
    path fam
    path sampleFile

    output:
    path 'shuffled.bed', emit: bed
    path 'shuffled.bim', emit: bim
    path 'shuffled.fam', emit: fam

    script:
    """
    # 1. Create shuffled sample list
    create_shuffled_sample_list.R --sample_file ${sampleFile}
    # 2. Shuffle samples 
    plink2 --bed ${bed} \
        --bim ${bim} \
        --fam ${fam} \
        --indiv-sort f ShuffledSampleOrder.txt \
        --make-bed \
        --out shuffled \
        --threads 4
    """
}


workflow {
    concatCHRFiles(params.inputDir)
    fixGTAnnot(concatCHRFiles.output.sortedVcf)
    splitMultiAllelicVariants(fixGTAnnot.output.gtFixedVcf)
    // filterMultiAllelicVariants(splitMultiAllelicVariants.output) // optional step
    // createGraphs(splitMultiAllelicVariants.output.splitVcf) // optional step
    filterVariants(splitMultiAllelicVariants.output.splitVcf)
    convertToPlinkFormat(filterVariants.output.filteredVcf)
    calculateMissingness(convertToPlinkFormat.output, splitMultiAllelicVariants.output.splitVcf)
    createHetFile(convertToPlinkFormat.output)
    findMissingSamples(calculateMissingness.output.missing)
    filterMissingSamples(findMissingSamples.output.filteredSamples, convertToPlinkFormat.output)
    findHetSamples(createHetFile.output)
    filterHetSamples(findHetSamples.output.failedHetSamples, filterMissingSamples.output)
    filterLowAltFreq(filterHetSamples.output)
    filterRelated(filterLowAltFreq.output)
    // createMetricsFile(filterHetSamples.output) // mind the folder structure
    getVariantCount(splitMultiAllelicVariants.output.splitVcf, filterVariants.output.filteredVcf)
    getSampleCount(concatCHRFiles.output.sortedVcf, filterRelated.output.fam)
    convertBackToVCF(filterRelated.output.bed, filterRelated.output.bim, filterRelated.output.fam)
    runBeagle(convertBackToVCF.output.finalVcf)
    convertToPlinkFormatAlt(runBeagle.output.beagleVcf)
    popProject(filterRelated.output.bed, filterRelated.output.bim, filterRelated.output.fam)
    popProjectAfterBeagle(convertToPlinkFormatAlt.output.bed, convertToPlinkFormatAlt.output.bim, convertToPlinkFormatAlt.output.fam)
    targetPCA(filterRelated.output.bed, filterRelated.output.bim, filterRelated.output.fam)
    finalPCA(targetPCA.output.bed, targetPCA.output.bim, targetPCA.output.fam)  
}