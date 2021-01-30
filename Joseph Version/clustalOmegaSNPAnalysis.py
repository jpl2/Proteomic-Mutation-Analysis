import argparse
import time
import csv
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO


def main():
    # PARSE THE USER INPUT
    parser = argparse.ArgumentParser(description="Run Clustal Omega to perform MSA on sequences and identify the SNPs t"
                                     "hat occur with more than 5% mutation rate.")
    parser.add_argument("infile", metavar='I', type=str, help="The name of the fasta file to perform MSA on.")
    parser.add_argument("tree_out", metavar='G', type=str, help="The name of the output guidetree file.")
    parser.add_argument("align_out", metavar='A', type=str, help="The name of the output alignment file.")
    parser.add_argument("outfile", metavar='SO', type=str, help="The name of the output SNP analysis file (text).")
    parser.add_argument("snp_table", metavar='ST', type=str, help="The name of the output SNP analysis file (csv).")
    parser.add_argument("-PA", "--print_align", action="store_true", help="Print the MSA output to terminal.")
    parser.add_argument("-PR", "--print_result", action="store_true", help="Print the SNP analysis output to terminal.")
    args = parser.parse_args()

    # ESTABLISH AND RUN CLUSTAL OMEGA COMMAND
    print("RUNNING CLUSTAL OMEGA...")
    clustal_command = ClustalOmegaCommandline(infile=args.infile, outfile=args.align_out, guidetree_out=args.tree_out,
                                              outfmt='clustal', verbose=True)
    clustal_command()

    # READ IN THE ALIGNMENT
    print("READING ALIGNMENTS...")
    alignment = AlignIO.read(args.align_out, "clustal")

    # PRINT ALIGNMENT RESULTS TO TERMINAL
    if args.print_align:
        print(alignment)

    # IDENTIFY SNPS WITH MORE THAN 5% MUTATION RATE
    snps = []
    thresh = 0.05 * (len(alignment) - 1)
    for i in range(0, len(alignment[0].seq)):
        sn_aln = alignment[:, i]
        ref = sn_aln[0]
        bases = {'A': 0, 'G': 0, 'C': 0, 'T': 0, '-': 0, 'N': 0, 'Y': 0}
        for j in range(1, len(sn_aln)):
            bases[sn_aln[j]] += 1
        for base in bases:
            if base != ref and bases[base] > thresh:
                location = i + 1
                rate = (bases[base] / (len(alignment) - 1)) * 100
                snp = (location, ref, base, rate)
                snps.append(snp)

    # WRITE SNP ANALYSIS TO CSV
    print("WRITING RESULTS TO CSV...")
    with open(args.snp_table, "w") as file:
        table = csv.writer(file)
        table.writerow(['Location', 'Reference Base', 'Mutation Base', 'Percent of Mutation'])
        table.writerows(snps)

    # WRITE SNP ANALYSIS TO OUTPUT TEXT FILE
    print("WRITING RESULTS TO OUTFILE...")
    with open(args.outfile, "w") as out:
        out.write("SINGLE NUCLEOTIDE POLYMOPHISMS FOR " + args.infile)
        out.write("\nREFERENCE SEQUENCE: " + alignment[0].id)
        out.write("\nTOTAL SNPS IDENTIFIED: " + str(len(snps)))
        for i in range(0, len(snps)):
            out.write("\n\nSNP " + str(i))
            out.write("\n\tLOCATION: " + str(snps[i][0]))
            out.write("\n\tMUTATION: " + snps[i][1] + " -> " + snps[i][2])
            out.write("\n\tRATE OF MUTATION: " + str(snps[i][3]) + "%")

    # PRINT SNP ANALYSIS RESULTS TO TERMINAL
    if args.print_result:
        with open(args.outfile, "r") as out:
            print(out.read())


if __name__ == '__main__':
    # START THE PROGRAM EXECUTION TIMER
    start = time.time()

    # EXECUTE PROGRAM
    main()

    # REPORT THE FINAL EXECUTION TIME
    print("EXECUTION TOOK %s SECONDS" % (time.time() - start))
