#!/usr/bin/env python
'''
Convert a bed (https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file to the eland format (http://support.illumina.com/sequencing/documentation.ilmn) as used by the illumina CASAVA pipeline 1.8+ 
Author: Julius Muller 2010
'''

from argparse import ArgumentParser
from warnings import warn

parser=ArgumentParser(description='Convert a bed file to the eland format as used by the illumina CASAVA pipeline')
parser.add_argument('BEDFILE',type=argparse.FileType('r'),help='Input BED file (use - from standard input)')
parser.add_argument('--sequence_dummy','-s',action="store_true",help='A dummy sequence and quality string will be generated [default: blank to reduce file size]')
parser.add_argument('--old','-o',action="store_true",help="Set quality string to the maximal value as specified in CASAVA 1.8 and earlier [default: CASAVA 1.8+]")
parser.add_argument('--result_format','-r',action="store_true",help="Output the eland_result.txt format [default: eland_export.txt format]")

args=parser.parse_args();
dseq=""
qseq=""
qval="h" if args.old else "I"
for lnum,txt in enumerate(args.BEDFILE):    
  txt=txt.strip('\n').split('\t')
	if len(txt)<6:
	    print txt
	    warn(("Columns < 6, Line {} skipped").format(lnum+1))
	    continue
	mchr,pos,pos2,sname,score,strand=txt[0:6]
	seqlen=int(pos2)-int(pos)
	if args.sequence_dummy:dseq="".join(["A" for x in xrange(seqlen)])
	
	if strand=='+':strand='F' 
	elif strand=='-':strand='R'
	else:exit('Direction error! Must be + or - but is: %s'%strand)
	mchr=mchr+".fa"
	if args.result_format:eland_ff="{0}\t{1}\tU0\t1\t0\t0\t{2}\t{3}\t{4}\t..".format(sname,dseq,mchr,pos,strand)
	else:
		if args.sequence_dummy:qseq="".join([qval for x in xrange(seqlen)])
		eland_ff="{0}\t\t\t\t\t\t0\t1\t{1}\t{2}\t{3}\t\t{4}\t{5}\t{6}\t126\t\t\t\t\t\t".format(sname,dseq,qseq,mchr,pos,strand,seqlen)
	print(eland_ff)  



