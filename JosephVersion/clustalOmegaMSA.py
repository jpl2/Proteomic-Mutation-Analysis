import argparse
import time
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO


def main():
    # PARSE THE USER INPUT
    parser = argparse.ArgumentParser(description="Run Clustal Omega to perform MSA on sequences.")
    parser.add_argument("infile", metavar='I', type=str, help="The name of the fasta file to perform MSA on.")
    parser.add_argument("tree_out", metavar='G', type=str, help="The name of the output guidetree file.")
    parser.add_argument("align_out", metavar='A', type=str, help="The name of the output alignment file.")
    parser.add_argument("-P", "--print", action="store_true", help="Print the output to terminal.")
    args = parser.parse_args()

    # ESTABLISH AND RUN CLUSTAL OMEGA COMMAND
    print("RUNNING CLUSTAL OMEGA...")
    clustal_command = ClustalOmegaCommandline(infile=args.infile, outfile=args.align_out, guidetree_out=args.tree_out,
                                              outfmt='clustal', verbose=True)
    clustal_command()

    # PRINT RESULTS TO TERMINAL
    if args.print:
        print(AlignIO.read(args.align_out, "clustal"))


if __name__ == '__main__':
    # START THE PROGRAM EXECUTION TIMER
    start = time.time()

    # EXECUTE PROGRAM
    main()

    # REPORT THE FINAL EXECUTION TIME
    print("EXECUTION TOOK %s SECONDS" % (time.time() - start))
