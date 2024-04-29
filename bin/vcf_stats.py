#!/usr/bin/env python3

"""
    usage:
        ./vcf_rompare.py -f [INPUT VCF FILE] -r [VCF COMPARISON FILE]
        -ih [INPUT HEADERS TEXT FILE] -ch [COMPARISON HEADER TEXT FILE]
        -chr [CHROMOSOME TEXT FILE] -n [CHUNKSIZE (Variant amount per run)]
"""

# METADATA VARIABLES
__author__ = "Orfeas Gkourlias"
__status__ = "WIP"
__version__ = "0.1"

# IMPORTS
import sys
import argparse
import gzip

# CLASSES
class VcfStat:
    """Class to handle VCF file and retrieve variant & dosgae statistics."""
    def __init__(self, input_handle, output_handle):
        self.input_handle = input_handle
        self.output_handle = output_handle

        self.tsv_head_str = (
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
            + "alt_freq\n"        
        )
        
    def walk(self):
        """Walk through the VCF file and retrieve the variants."""
        with gzip.open(self.input_handle, "rt") as vcf_f, gzip.open(self.output_handle, "wt") as tsv_f:
            tsv_f.write(self.tsv_head_str)
            ctr = 0
            for line in vcf_f:
                if line.startswith("#") != True:
                    tsv_f.write(self.record_extract(line))
                    if ctr % 10000 == 0:
                        print(f"{ctr} Variants Written.", end = "\r")
                    ctr += 1
            print(f"A Total of {ctr} Variants Were Written. Done!")

    def record_extract(self, variant):
        """Extracts the relevant intersection information from two tabix records.

        Args:
            i_record (list): list containing record entries for input file.
            headdx (dict): header indexes

        Returns:
            string: string to be written to output file.
        """
        # Get a column value for each row if header/column index is present in both files.
        variant = variant.split("\t")
        # sample_size and indel variables.
        sample_size = is_indel = 0
        nr_aa = nr_ab = nr_bb = a_freq = b_freq = 0
        missing = gt_total = 0
        # Genotype lists for use in pearson corr.
        gt_list = []
        gq = 0
        # Construct identifier string for the record/row.
        var_id = f"{variant[0]}:{variant[1]}:{variant[3]}_{variant[4]}"
        # Selection is same length because of column index selection.
        # Zip to iterate over both record values and track index.
        for i, entry in enumerate(variant):
            match i:  # Match the index value.
                case 0 | 1 | 2 | 5 | 7:  # Values on these indexes are not relevant.
                    continue
                
                case 3 | 4:  # Check whether ref or alt contains an indel.
                    if len(entry) > 1:
                        is_indel = 1
                
                case 6:
                    gq = entry

                case 8:  # Initialse the format column indicators/headers.
                    val_dict = {}

                    # Store format indication/header into a dictionary with empty values.
                    for ind in entry.split(":"):
                        val_dict[ind] = "./."

                case other:  # Sample rows.
                    gt_total += 1
                    # Assign values in format order.
                    for key, val in zip(
                        val_dict, entry.split(":")
                    ):
                        val_dict[key] = val

                    # Save genotypes in seperate var.
                    gt = val_dict["GT"]
                    ad = val_dict["AD"]

                    if "." in gt:
                        missing += 1

                    # Split and save genotype numbers for sum.
                    # Example: "1/0" becomes list: [1,0]
                    vals = gt.replace("|", "/").split("/")

                    # Check whether genotypes are the same between files.
                    # Also check whether they are not missing genotypes.
                    if vals[0] != ".":
                        sample_size += 1
                        gt = int(vals[0]) + int(vals[1])
                        gt_list.append(gt)

                        match (gt):
                            case 0:
                                nr_aa += 1
                            case 1:
                                nr_ab += 1
                            case 2:
                                nr_bb += 1
                        

        maf, maf_type = self.calc_maf(nr_aa, nr_ab, nr_bb)
        a_freq, b_freq = self.get_freqs(maf, maf_type)
        # Get HWE
        hwe = self.calc_hwe(nr_aa, nr_ab, nr_bb)

    # Return a constructed row string for immediate output writing.
        return (
            f"{var_id}\t{is_indel}\t{sample_size}\t"
            f"{nr_aa}\t{nr_ab}\t{nr_bb}\t"
            f"{maf}\t"
            f"{hwe}\t"
            f"{missing}\t{gt_total}\t"
            f"{maf_type}\t"
            f"{gq}\t"
            f"{a_freq}\t"
            f"{b_freq}\n"
        )

    def calc_maf(self, nr_aa, nr_ab, nr_bb):
        sample_size = nr_aa + nr_ab + nr_bb
        if sample_size == 0:
            return 0, "nan"
        nrB = (nr_bb * 2) + nr_ab
        nrAlleles = sample_size * 2
        maf = nrB / nrAlleles
        if maf > 0.5:
            return 1 - maf, "REF"  # reference allele is the minor allele
        return maf, "ALT"  # alternate allele is the minor allele
    
    def calc_hwe(self, obs_hom1, obs_hets, obs_hom2):
        obs_homc = obs_hom1
        if obs_hom1 < obs_hom2:
            obs_homc = obs_hom2
        obs_homr = obs_hom2
        if obs_hom1 < obs_hom2:
            obs_homr = obs_hom1

        rare_copies = 2 * obs_homr + obs_hets
        l_genotypes = obs_hets + obs_homc + obs_homr

        if l_genotypes == 0:
            return -1

        het_probs = [0] * (rare_copies + 1)
        mid = int(rare_copies) * int(2 * l_genotypes - rare_copies) / int(2 * l_genotypes)
        mid = int(mid)
        if mid % 2 != rare_copies % 2:
            mid += 1
        mid = int(mid)

        curr_hets = mid
        curr_homr = (rare_copies - mid) / 2
        curr_homc = l_genotypes - curr_hets - curr_homr
        het_probs[int(mid)] = 1.0
        sum = het_probs[int(mid)]

        curr_hets = int(mid)
        while curr_hets > 1:
            het_probs[curr_hets - 2] = (
                het_probs[curr_hets]
                * curr_hets
                * (curr_hets - 1.0)
                / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
            )
            sum += het_probs[curr_hets - 2]
            curr_homr += 1
            curr_homc += 1
            curr_hets -= 2

        curr_hets = int(mid)
        curr_homr = (rare_copies - mid) / 2
        curr_homc = l_genotypes - curr_hets - curr_homr
        while curr_hets <= (rare_copies - 2):
            het_probs[curr_hets + 2] = (
                het_probs[curr_hets]
                * 4.0
                * curr_homr
                * curr_homc
                / ((curr_hets + 2.0) * (curr_hets + 1.0))
            )
            sum += het_probs[curr_hets + 2]
            curr_homr -= 1
            curr_homc -= 1
            curr_hets += 2

        i = 0
        while i <= rare_copies:
            het_probs[i] /= sum
            i += 1

        p_hwe = 0.0
        i = 0
        while i <= rare_copies:
            if het_probs[i] <= het_probs[obs_hets]:
                p_hwe += het_probs[i]
            i += 1

        if p_hwe > 1:
            p_hwe = 1
        return p_hwe
    
    def get_freqs(self, maf, maf_type):
        # Write headers to output file.

        if maf_type == "ALT":
            b_freq = maf
            a_freq = 1 - maf

        else:
            a_freq = maf
            b_freq = 1 - maf
            
        
        return a_freq, b_freq
                    

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
    parser.add_argument("-i", "--input", type=str)
    parser.add_argument("-o", "--output", type=str)
    return parser.parse_args()


# MAIN
def main(args):
    """Main function"""
    # Get args.
    args = arg_parse()
    # Perform comparison and write output.
    vcf_stat = VcfStat(args.input, args.output)
    vcf_stat.walk()
    # FINISH
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
