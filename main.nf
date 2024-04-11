"""
File:         main.nf
Created:      2024/04/11
Last Changed: 2024/04/11
Author:       Peter Riesebos
"""

nextflow.enable.dsl=2


process concatCHRFiles {
    publishDir "${params.out_dir}", mode: 'move'
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
    publishDir "${params.out_dir}", mode: 'move'
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
    publishDir "${params.out_dir}", mode: 'move'
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

workflow {
    concatCHRFiles()
    splitMultiAllelicVariants()
    filterVariants()
    convertToPlinkFormat()
    calculateMissingness() // also show hwe / freq / etc ?
    // remove variants if sample has missigness >0.25
    // remove high heterozygosity samples (+/- 3 SD from the mean)
    // identity by state (IBS)
    // Bigsnpr (map sample genotypes PCs against 1000G to assign likely ancestry), wat moet hier nog mee gebeuren? Alleen EUR / alleen super-pop?
}