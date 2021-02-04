import argparse
import time
import csv
from Bio import AlignIO
from Bio.Data import CodonTable


def main():
    # PARSE THE USER INPUT
    parser = argparse.ArgumentParser(description="Identify the SNPs in the alignment file that occur with more than 5% "
                                     "mutation rate.")
    parser.add_argument("align", metavar='A', type=str, help="The name of the input alignment file.")
    parser.add_argument("genome", metavar='G', type=str, help="The name of the input genome CSV file.")
    parser.add_argument("outfile", metavar='O', type=str, help="The name of the output nonsynonymous SNP analysis file "
                        "(text).")
    parser.add_argument("snp_table", metavar='ST', type=str, help="The name of the output nonsynonymous SNP analysis fi"
                        "le (csv).")
    parser.add_argument("-P", "--print", action="store_true", help="Print the output to terminal.")
    args = parser.parse_args()

    # READ IN THE ALIGNMENT
    print("READING ALIGNMENT...")
    alignment = AlignIO.read(args.align, "clustal")
    ref_seq = alignment[0].seq

    # PRINT ALIGNMENTS TO TERMINAL
    if args.print:
        print(alignment)
        print()

    # READ IN THE GENOME INFORMATION
    print("READING GENOME DATA...")
    genomes = []
    with open(args.genome, "r") as file:
        genome_reader = csv.reader(file)
        fields = next(genome_reader)
        for genome in genome_reader:
            genomes.append(genome)

    # PRINT GENOME INFORMATION TO TERMINAL
    if args.print:
        print(genomes)
        print()

    # IDENTIFY SNPS WITH MORE THAN 5% MUTATION RATE
    snps = []
    thresh = 0.05 * (len(alignment) - 1)
    for i in range(0, len(ref_seq)):
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

    # DISCARD SNPS TO GAP, N, AND Y
    valid_snps = []
    for snp in snps:
        if snp[2] != '-' and snp[2] != 'N' and snp[2] != 'Y':
            valid_snps.append(snp)
    snps.clear()

    # ANALYZE REMAINING SNPS
    proteome_snps = []
    for snp in valid_snps:
        protein = []
        gene_len = 0

        # CHECK IF SNP IN GENOME; KEEP TIGHTEST BOUND GENOME
        for genome in genomes:
            if not genome[3]:
                local_len = float('inf')
            else:
                local_len = float(genome[3])
            if float(genome[1]) <= snp[0] <= float(genome[2]) and (gene_len == 0 or local_len < gene_len):
                protein = genome
                gene_len = local_len

        # IF NO CONTAINING GENOME, SNP IS NONENCODING
        if gene_len == 0:
            p_snp = ['N/A', snp[0], 0, snp[1], snp[2], 'N/A', 'N/A', 'NONENCODING', snp[3]]
            proteome_snps.append(p_snp)

        # OTHERWISE, DETERMINE IF SNP IS SYNONYMOUS
        else:
            gene_loc = snp[0] - round(float(protein[1])) + 1

            # EXPAND REF CHARACTER FOR SNP INTO FULL CODON
            codon_pos = gene_loc % 3
            ref_codon = snp[1]
            if codon_pos == 1:
                ref_codon = ref_codon + ref_seq[snp[0]] + ref_seq[snp[0] + 1]
            elif codon_pos == 2:
                ref_codon = ref_seq[snp[0] - 2] + ref_codon + ref_seq[snp[0]]
            else:
                ref_codon = ref_seq[snp[0] - 3] + ref_seq[snp[0] - 2] + ref_codon

            # BUILD MUTATED CODON
            mut_codon = ref_codon
            mut_codon = mut_codon[:((codon_pos - 1) % 3)] + snp[2] + mut_codon[(((codon_pos - 1) % 3) + 1):]

            # TRANSLATE CODONS WITH CODON TABLE IF NOT STOP CODONS
            if ref_codon not in CodonTable.standard_dna_table.stop_codons:
                ref_amino = CodonTable.standard_dna_table.forward_table[ref_codon]
            else:
                ref_amino = "STOP"
            if mut_codon not in CodonTable.standard_dna_table.stop_codons:
                mut_amino = CodonTable.standard_dna_table.forward_table[mut_codon]
            else:
                mut_amino = "STOP"

            # BUILD ANALYZED SNP
            p_snp = [protein[0], snp[0], gene_loc, snp[1], snp[2], ref_amino, mut_amino]
            if ref_amino == mut_amino:
                p_snp.append("SYNONYMOUS")
            else:
                p_snp.append("NONSYNONYMOUS")
            p_snp.append(snp[3])
            proteome_snps.append(p_snp)

    # FREE MEMORY OF PREVIOUS LIST
    valid_snps.clear()

    # WRITE SNP ANALYSIS TO CSV
    print("WRITING RESULTS TO CSV...")
    with open(args.snp_table, "w") as file:
        table = csv.writer(file)
        table.writerow(['Protein', 'Exact Location', 'Protein Location', 'Reference Base', 'Mutation Base',
                        'Reference Amino Acid', 'Mutation Amino Acid', 'Type', 'Percent of Mutation'])
        table.writerows(proteome_snps)

    # WRITE SNP ANALYSIS TO OUTPUT TEXT FILE
    print("WRITING RESULTS TO OUTFILE...")
    with open(args.outfile, "w") as out:
        out.write("SINGLE NUCLEOTIDE POLYMOPHISMS (> " + str((thresh / (len(alignment) - 1)) * 100) + "%) FOR " +
                  args.align)
        out.write("\nREFERENCE SEQUENCE: " + alignment[0].id)
        out.write("\nTOTAL VALID SNPS IDENTIFIED: " + str(len(proteome_snps)))
        for i in range(0, len(proteome_snps)):
            out.write("\n\nSNP " + str(i+1))
            out.write("\n\tREFERENCE SEQUENCE LOCATION: " + str(proteome_snps[i][1]))
            out.write("\n\tCONTAINING GENOME: " + str(proteome_snps[i][0]))
            out.write("\n\tLOCATION IN GENOME: " + str(proteome_snps[i][2]))
            out.write("\n\tNUCLEOTIDE MUTATION: " + proteome_snps[i][3] + " -> " + proteome_snps[i][4])
            out.write("\n\tCORRESPONDING AMINO ACID MUTATION: " + proteome_snps[i][5] + " -> " + proteome_snps[i][6])
            out.write("\n\tTYPE OF MUTATION: " + proteome_snps[i][7])
            out.write("\n\tRATE OF MUTATION: " + str(proteome_snps[i][8]) + "%")

    # PRINT SNP ANALYSIS RESULTS TO TERMINAL
    if args.print:
        with open(args.outfile, "r") as out:
            print(out.read())


if __name__ == '__main__':
    # START THE PROGRAM EXECUTION TIMER
    start = time.time()

    # EXECUTE PROGRAM
    main()

    # REPORT THE FINAL EXECUTION TIME
    print("EXECUTION TOOK %s SECONDS" % (time.time() - start))
