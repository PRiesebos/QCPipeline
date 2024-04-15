"""
File:         main.nf
Created:      2024/04/11
Last Changed: 2024/04/15
Author:       Peter Riesebos & Orfeas Gkourlias
"""

nextflow.enable.dsl=2


process concatCHRFiles {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    // Parameterize the chromosome range
    int chrom_start = 1
    int chrom_end = 22

    input:
    path chr_files from Channel.fromPath(params.input_dir).filter { file -> file.name.endsWith('.vcf.gz') }

    output:
    path "${params.out_dir}/merged_output.vcf.gz" into merged_output

    script:
    '''
    bcftools concat $(printf "${chr_files} ") -Oz -o ${params.out_dir}/merged_output.vcf.gz
    '''
}

process splitMultiAllelicVariants {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path filteredVcfFile

    output:
    path "${filteredVcfFile.SimpleName}.vcf.gz"

    script:
    '''
    bcftools norm -m -any -o norm_${filteredVcfFile.SimpleName}.vcf.gz -Oz ${filteredVcfFile}
    '''
}


process filterVariants {
    publishDir "${params.out_dir}", mode: 'move'
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path vcfFile

    output:
    path "*-filtered.vcf.gz", emit: filteredVcfFile
    path "*.log"

    shell:
    '''
    # 1. Define command arguments
    commandArguments="--input !{vcfFile} \
    --output no_multi_allelic \
    --call_rate 0.5 \
    --filtered_depth 5 \
    --genotype_quality 10 \
    --minor_allele_frequency 0.01 \
    --no_snv_vqsr_check \
    --no_indel_vqsr_check \
    --remove_non_pass_snv \
    --remove_non_pass_indel \
    --replace_poor_quality_genotypes

    # 2. Run the custom VCF filter script
    python3 custom_vcf_filter.py ${commandArguments} \
    | tee custom_vcf_filter.log
    '''
}

process convertToPlinkFormat {
    publishDir "${params.out_dir}", mode: 'move'
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path vcfFile

    output:
    path "*.bed"
    path "*.bim"
    path "*.bam"

    script:
    '''
    plink2 --vcf ${vcfFile} --make-bed
    '''
}

process calculateMissingness {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path plinkData

    output:
    path "*.imiss"

    script:
    '''
    plink2 --bfile ${plinkData} --missing
    '''
}

// Create a txt file with samples where missingness => 50%
// Filter these samples. Output -> filtered plink files
process filterMissingSamples {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    file smiss from inputSmissFile
    file data from inputDataFile

    output:
    file 'filtered_samples.txt' into filteredSamplesFile
    file 'data_keep.bed' into bedFile
    file 'data_keep.bim' into bimFile
    file 'data_keep.fam' into famFile

    script:
    """
    python filter_samples.py ${smiss} filtered_samples.txt --threshold 0.5
    plink2 --bfile ${data} --keep filtered_samples.txt --make-bed --out data_keep
    """
}

process createHetFile {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    path smiss

    output:
    path plinkData

    script:
    '''
    plink --bfile your_data --het cols=hom,het,nobs,f --out heterozygosity_output
    '''
}

process filterHetFile {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 1

    input:
    file het_file
    path plinkData
    file failed_samples_txt from upstream

    output:
    file "*.png"
    file "*Failed.txt"
    file "*FailedSamplesOnly.txt"
    path plinkData

    script:
    '''
    python3 heterozygosityCheck.py ${het_file} ${output}
    plink2 --bfile ${plinkData}/input_data \
       --remove ${failed_samples_txt} \
       --make-bed \
       --out ${plinkData}/filtered_output
    '''
}

process filterRelated {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 4

    input:
    path vcf

    output:
    path 'RelatednessCheck.bed', emit: bedFile
    path 'RelatednessCheck.bim', emit: bimFile
    path 'RelatednessCheck.fam', emit: famFile
    path 'related.kin0'
    path 'RelatednessPassedSamples.txt'

    script:
    """
    # 1. Do relatedness check
    ~/plink2 --vcf ${plinkData} \
        --make-king-table \
        --king-table-filter ${params.kingTableFilter} \
        --out related \
        --threads 4 \

    # 2. Create sample list of non-related samples
    Rscript find_related_samples.R --kin_file related.kin0 --target_bed ${bedFile}

    # 3. Remove samples that are not on the list created above
    ~/plink2 --bed ${bedFile} \
        --bim ${bimFile} \
        --fam ${famFile} \
        --keep RelatednessPassedSamples.txt \
        --make-bed \
        --out RelatednessCheck
    """
}

process popProject {
    errorStrategy 'retry'
    maxRetries 1

    time '4h'
    memory '8 GB'
    cpus 4

    input:
    path vcf
    
    output:
    path '*.png'
    path '*.pdf'
    path '1000G_PC_projections.txt'
    path 'PopAssignResults.txt'
    
    script:
    """
    Rscript project_samples_to_superpop.R --ref_bed ${params.refPath}.bed \
        --target_bed ${bedFile} --ref_pop ${params.refPop}
    """
}

workflow {
    concatCHRFiles()
    splitMultiAllelicVariants()
    filterVariants() // Save output?
    convertToPlinkFormat()
    filterMissingSamples() // default threshold of 50%
    createHetFile()
    calculateMissingness() 
    popProject()
    // also show hwe / freq / etc ?
    // remove variants if sample has missigness >0.25
    // remove high heterozygosity samples (+/- 3 SD from the mean)
    // identity by state (IBS)
    // Bigsnpr (map sample genotypes PCs against 1000G to assign likely ancestry), wat moet hier nog mee gebeuren? Alleen EUR / alleen super-pop?

    // Add plots between steps!
}