import sys
import simplesam
from simplesam import Reader, Writer

BAMFILE=sys.argv[1]
SAMFILE=sys.argv[2]

in_file = open(BAMFILE, 'r')
in_sam = Reader(in_file)
x = next(in_sam)
x.tags

with Reader(open(BAMFILE)) as in_bam:
	with Writer(open(SAMFILE, 'w'), in_bam.header) as out_sam:
		for read in in_bam:
			read["UB"] = read.qname.split(":")[2] # add the umi tag
			read["CB"] = read.qname.split(":")[1] # add the barcode tag
			out_sam.write(read)