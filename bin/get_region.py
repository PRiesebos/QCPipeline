#!/usr/bin/env python3

"""
    usage:
        ./vcf_rompare.py -i [INPUT VCF.GZ FILE] -G [GTF REF FILE] -o [OUTPUT.TSV.GZ FILE]
"""

# METADATA VARIABLES
__author__ = "Orfeas Gkourlias"
__status__ = "WIP"
__version__ = "0.1"

# IMPORTS
import sys
import argparse
import csv
import gzip
import re
from GTFAnnotation import GTFAnnotation
from features.chromosome import Chromosome

# FUNCTIONS
def arg_parse():

    """
    Function to parse arguments.

    :file: Input VCF file.
    :reference: VCF file to be compared against.
    :input_header: Text file containing tab seperated headers.
    :comp_header: Text file containing tab seperated headers for compare file.
    :chr: Test file named as input chromosome.
    :return: Argparse namespace object.
    """

    # Parser init
    parser = argparse.ArgumentParser(prog="VcfCompare")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "-all",
        action="store_true",
        help="Run the program in Server mode; see extra options needed below",
    )
    mode.add_argument(
        "-comp",
        action="store_true",
        help="Run the program in Client mode; see extra options needed below",
    )
    parser.add_argument("-i", "--tsv", type=str)
    parser.add_argument("-g", "--gtf", type=str)
    parser.add_argument("-o", "--output", type=str)
    return parser.parse_args()

class get_region:
    def __init__(self, mode, tsv, gtf, output):
        self.mode = mode
        self.tsv = tsv
        self.output = output
        self.annotation = GTFAnnotation(gtf)
    
    def get_genes(self):
        self.genes_by_chr = self.annotation.getGenesByChromosome()
    
    def get_regions(self):
        with gzip.open(self.output, "wt") as fo:
            if self.mode == "comp":
                fo.write(
                "var_id\t"
                + "is_indel\t"
                + "sample_size\t"
                + "nrhoma_i\t"
                + "nrhets_i\t"
                + "nrhomb_i\t"
                + "nrhoma_r\t"
                + "nrhets_r\t" 
                + "nrhomb_r\t"
                + "matching\t"
                + "match_ratio\t"
                + "maf_i\t"
                + "maf_r\t"
                + "hwe_i\t"
                + "hwe_r\t"
                + "p_corr\t"
                + "s_corr\t"
                + "missing_i\t"
                + "missing_r\t"
                + "gt_total\t"
                + "maf_type_i\t"
                + "maf_type_r\t"
                + "pass_r\t"
                + "gq_i\t"
                + "gq_r\t"
                + "gene_types\t"
                + "transcripts\t"
                + "exons\t"
                + "var_status\n"
                )
            
            elif self.mode == "all":
                fo.write(
                "var_id\t"
                + "is_indel\t"
                + "sample_size\t"
                + "nr_aa\t"
                + "nr_ab\t"
                + "nr_bb\t"
                + "maf\t"
                + "hwe\t"
                + "missing\t"
                + "gt_total\t"
                + "maf_type\t"
                + "gq\t"
                + "ref_freq\t"
                + "alt_freq\t"
                + "af_info\t"
                + "maf_info\t"
                + "r2_info\t"    
                + "gene_types\t"
                + "transcripts\t"
                + "exons\t"
                + "var_status\n"        
                )

            self.get_genes()
            with gzip.open(self.tsv, "rt") as f:
                spamreader = csv.reader(f, delimiter='\t')
                next(spamreader, None)
                for row in spamreader:
                    variant = row
                    pos = int(variant[0].split(":")[1])
                    chr_str = str(variant[0].split(":")[0].replace("chr", ""))
                    chrom = Chromosome.parse(chr_str)
                    genes = self.genes_by_chr[chrom]
                    overlapGenes = []
                    gene_types = []
                    overlapTranscripts = []
                    overlapExons = []
                    region_status = "nan"
                    for gene in genes:
                        if pos >= gene.start and pos <= gene.stop:
                            overlapGenes.append(gene)
                        if pos > gene.stop:
                            continue

                    if len(overlapGenes) > 0:
                        for gene in overlapGenes:
                            gene_types.append(gene.type)
                            transcripts = gene.transcripts
                            for transcript in transcripts:
                                if pos >= transcript.start and pos <= transcript.stop:
                                    overlapTranscripts.append(transcript.name)
                                    if region_status != "exonic":
                                        region_status = "intronic"
                                    exons = transcript.exons
                                    for exon in exons:
                                        if pos >= exon.start and pos <= exon.stop:
                                            overlapExons.append(exon.name)
                                            region_status = "exonic"
                    
                    overlapTranscripts = ":".join(overlapTranscripts)
                    overlapExons = ":".join(overlapExons)
                    gene_types = ":".join(gene_types)
                    variant.append(gene_types)
                    variant.append(overlapTranscripts)
                    variant.append(overlapExons)
                    variant.append(region_status)
                    variant_str = "\t".join(variant) + "\n"
                    fo.write(variant_str) 
    
# MAIN
def main(args):
    """Main function"""
    # Get args.
    args = arg_parse()
    if args.comp:
        mode = "comp"

    elif args.all:
        mode = "all"

    getter = get_region(mode, args.tsv, args.gtf, args.output)
    getter.get_regions()
    # FINISH
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
