# Generate cut site positions in genome from given restriction enzymes (Multiple enzymes and 'N's are supported!)
# Author: Xiang Zhou
# Affiliation: Arima Genomics Inc.

# Usage: generate_site_positions_Arima.py [-h] -i INPUT -e ENZYMES [ENZYMES ...] -o OUTPUT
# Arguments:
#   -h, --help                show this help message and exit
#   -i INPUT                  Input FASTA filename
#   -e ENZYMES [ENZYMES ...]  Enzyme sequences
#   -o OUTPUT                 Output filename

import sys
import re
import argparse
import os.path
from Bio import SeqIO

def main():
	parser = argparse.ArgumentParser(description = "Generate cut site positions in genome from given restriction enzymes (Multiple enzymes and 'N's are supported!)")

	parser.add_argument('-i', metavar='INPUT', type=str, dest='input', help='Input FASTA filename', required=True)
	parser.add_argument('-e', metavar='ENZYMES', type=str.upper, nargs='+', dest='enzymes', help='Enzyme sequences', required=True)
	parser.add_argument('-o', metavar='OUTPUT', type=str, dest='output', help='Output filename', required=True)

	args = parser.parse_args()

	with open(args.input, 'rU') as FASTA, open(args.output, 'w') as OUTFILE:
		for record in SeqIO.parse(FASTA, "fasta"):
			length = len(record.seq)
			positions = find_re_sites(str(record.seq), args.enzymes)

			if positions:
				OUTFILE.write("{} {} {}\n".format(record.id, " ".join(map(str, positions)), length))
			else:
				OUTFILE.write("{} {}\n".format(record.id, length))

def find_re_sites(sequence, enzymes):
	#positions = []
	regex = "|".join(enzymes)
	regex = regex.replace("N", "[ACGT]")
	print("Finding cut sites matching pattern: {} ...".format(regex))
	iter = re.finditer(regex, sequence, re.IGNORECASE)
	#print([1+m.start() for m in iter])

	return [1+m.start() for m in iter]

if __name__ == '__main__':
	main()
