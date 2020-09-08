
import sys

def handleArgs():
	
	if len(sys.argv) != 4:
		print("ERROR: Expected 3 arguments: input gff, input fasta, output bed", file=sys.stderr)
		sys.exit(1)
	
	input_gff_filename = sys.argv[1]
	input_fasta_filename = sys.argv[2]
	output_bed_filename = sys.argv[3]

	return input_gff_filename, input_fasta_filename, output_bed_filename

def getHeaderLengths(fasta_fn):
	
	lengths = {}

	with open(fasta_fn, 'r') as fasta_fd:
		line = fasta_fd.readline()
		while line != '':
			header = line.rstrip('\n')[1:]
			seqlen = 0
			line = fasta_fd.readline()
			while line != '' and line[0] != '>':
				seqlen += len(line.rstrip('\n'))
				line = fasta_fd.readline()

			lengths[header] = seqlen

	return lengths

if __name__ == "__main__":

	igff_fn, ifa_fn, obed_fn = handleArgs()

	header_lengths = getHeaderLengths(ifa_fn)

	with open(obed_fn, 'w') as ofd:
		with open(igff_fn, 'r') as ifd:
			for line in ifd:
				fields = line.rstrip('\n').split('\t')
				if len(fields) >= 5:
					header = fields[0]
					feature = fields[2]
					start = int(fields[3]) # 1-based, but we'll want it to be 0-based and inclusive
					stop = int(fields[4]) # 1-based, but we'll want it to be 0-based and exclusive

					if feature == "mRNA":

						# set start/stop to ~1000 bp on either side
						start -= 1000
						stop += 1000

						# modify start/stop to not exceed the limits of the sequence
						start = start if start >= 0 else 0

						if not header in header_lengths:
							print("ERROR: Houston, we have a problem.", file=sys.stderr)
							sys.exit(1)
						
						stop = stop if stop <= header_lengths[header] else header_lengths[header]

						ofd.write(f"{header}\t{start}\t{stop}\n")



