#!/usr/bin/python

from Bio import SeqIO
import argparse #replaced the depreciated optparse
import re


def opts():
    usage = "python jseq_tools.py [options] <input file> <file type>\n"
    parser = argparse.ArgumentParser(usage=usage)
    group = parser.add_mutually_exclusive_group()
    parser.add_argument("--fmts", action="store_true", default=False,
                       help="prints a list of compatible formats for counting")
    group.add_argument("--count", "-c", dest="c_length", type=int, default=False,
                       help="counts all sequences over <int> length")
    group.add_argument("--filter", "-f", dest="f_length", type=int ,default=False,
                       help="saves file of all sequences over <int> length")
    group.add_argument("--extract", "-e", dest="ident", default=False,
                       help="saves a file of sequences extracted from a list of ids")
    parser.add_argument("input_file", type=str,
                        help="The input sequence file")
    parser.add_argument("file_type", type=str,
                        help="The input file format")
    args = parser.parse_args()

    if args.fmts:
        print usage
        print "\nFormats:"
        print "\
        fasta-\tThe generic sequence file format where each record starts with an identifer\n\
        \tline starting with a '>' character, followed by lines of sequence.\n\
        fastq-\tA 'FASTA like' format used by Sanger which also stores PHRED sequence quality\n\
        \tvalues (with an ASCII offset of 33).\n\
        genbank-The GenBank or GenPept flat file format.\n\
        swiss-\tPlain text Swiss-Prot aka UniProt format.\n\
        tab-\tSimple two column tab separated sequence files, where each line holds a record's\n\
        \tidentifier and sequence. For example, this is used as by Aligent's eArray software\n\
        \twhen saving microarray probes in a minimal tab delimited text file.\n\
        qual-\tA 'FASTA like' format holding PHRED quality values from sequencing DNA, but no\n\
        \tactual sequences (usually provided in separate FASTA files).\n\
        uniprot-xml-\tThe UniProt XML format (replacement for the SwissProt plain text format\n\
        \twhich we call 'swiss')\n\
        \n\tFor a full listing check Biopython SeqIO doumentation.\n"
        exit(0)

    """
    try:
      >>> import jseq_tools_argparse
  in_file, fmt = args
    except ValueError:
        print "\nIncorrect number of arguments!\n"
        print usage
        exit(1)
    """
    return args


def counting(records, length):
    """
    Takes records and optionally a length.
    Returns a count of all seqs in file and when a length is provided also a count of all\
    seqs greater than or equal to length. 
    """

    long_records = []
    count = 0
    count_1 = 0
    for record in records:
        count += 1 #counts every sequence in file
        if length and len(record.seq) >= length:
            count_1 += 1 #counts every sequence longer than "length"        
    if length:
        print "%i sequences longer than %i contained in this fasta!" % (count_1, length)
        print "%i sequences contained in this fasta!" % count
    else:
        print "%i sequences countained in this fasta!" % count


def filter_contigs(records, length, in_file, fmt):
    """
    Takes sequence records and an integer length.
    Writes a file of all sequences greater or equal to lengh as a new file.
    """

    long_records = []
    count = 0
    count_1 = 0
    for record in records:
        count += 1 #counts every sequence in file
        if len(record.seq) >= length:
            long_records.append(record)
            count_1 += 1 #counts every sequence longer than "length"
    outname = in_file.split(".")
    outname = "".join(outname[:-1])+"_large."+outname[-1]
    SeqIO.write(long_records, outname, fmt)
    print "%i sequences longer than %ibps extracted from a total of %i sequences!" % (
        count_1, length, count)


def id_parse(ids_handle):
    """
    Takes a file of ids and returns the ids as a list.
    Currently assumes each id is on a seperate line.
    """
    ##todo## be more flexible with input. Eg comma seperated list in file or cmd line.

    ids_handle = open(ids_handle, "r")
    ids = []
    for ID in ids_handle.readlines():
        ID = ID.strip()
        ids.append(ID)
    
    return ids


def extract_by_id(records, ids, in_file, fmt):
    """
    Takes a record file and a list of ids and extracts the records corresponding to those ids.
    That subset of records are then written to file using the file name as a base.
    """

    ##todo## 
    #generalise the regular expression so that it matches any unique string within the seq ID.
    #currently the regex is selective for contig_1 over contig_10 or contig_111 but requires
    #a change from digit type to non-digit ("\d" --> "\D") and would therefore fail to distinguish 
    #1_contig_illumina from 11_contig_illumina.
    #...maybe giving the user an option to enter their own regex would be the best solution...

    #extracts the records based on a regex using the id.
    #assumes the id is not followed by a digit character.
    interesting_records = []
    for rec in records:
        for ID in ids:
            regex = re.compile(ID + "\D" + "|" + ID + "$")
            match = re.search(regex, rec.id)
            if match:
                interesting_records.append(rec)
    #for rec in interesting_records:
    #    print rec.description

    outname = in_file.split(".")
    outname = "".join(outname[:-1])+"_extracted_ids."+outname[-1]
    SeqIO.write(interesting_records, outname, fmt)


def main():
    """
    jseq_tools is a collection of simple bioinformatics record manipulation tools.
    This function decides which functions to run based on user command line input.
    """

    args = opts() #parse args from command line

    handle = open(args.input_file, "rb") #opens file
    records = SeqIO.parse(handle, args.file_type) #parses file

    if args.f_length:
        filter_contigs(records, args.f_length, args.input_file, args.file_type)
    elif args.ident:
        ids = id_parse(args.ident)
        extract_by_id(records, ids, args.input_file, args.file_type)
    elif args.c_length:
        counting(records, args.c_length)
    else:
        print "No command given."


print __name__
if __name__ == '__main__':
    main()
