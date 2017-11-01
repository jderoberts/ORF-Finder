import argparse
import re

#define expected commandline arguments, inc. defaults and help/usage instructions
help_text = "Locate ORFs in uploaded .fasta file(s) and return predicted protein sequences.\n\n"
note = "Author: James Roberts <james.roberts-6@postgrad.manchester.ac.uk>"
parser = argparse.ArgumentParser(description=help_text, epilog=note, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument(
    '-n',
    '--nonstandard',
    dest='nonstandard',
    action='store_true',
    help='allow nonstandard bases (N) in ORFs',
)
parser.add_argument(
    '-v',
    '--overlap',
    dest='allow_overlap',
    action='store_true',
    help='allow overlapping ORFs',
)
parser.add_argument(
    '-s',
    '--summary',
    dest='summary',
    action='store_true',
    help='display summary information upon completion',
)
parser.add_argument(
    '-o',
    '--output',
    dest='output',
    action='store',
    type=str,
    default='',
    help='name of output file (default StdOut)',
    metavar='x'
)
parser.add_argument(
    '-f',
    '--frames',
    dest='frames',
    action='store',
    type=str,
    default='123456',
    help='frames to check for ORFs (default 123456, 4-6 == rc1-3)',
    metavar='x'
)
parser.add_argument(
    '-m',
    '--minlength',
    dest='minlength',
    action='store',
    type=int,
    default=150,
    help='minimum length for detected ORFs (default 150nt)',
    metavar='x'
)
parser.add_argument(
    '-c',
    '--maxlength',
    dest='maxlength',
    action='store',
    type=int,
    default=2500,
    help='maximum length for detected ORFs (default 2500nt)',
    metavar='x'
)
parser.add_argument(
    'filenames',
    metavar='input',
    type=str,
    nargs='+',
    help='input .fasta filename(s)')

arguments = parser.parse_args()

class DNASequence (object) :
    """Object to contain data for each sequence found."""
    def __init__(self, header, sequence) :
        #remove any trailing whitespace on header row
        self.header = header.rstrip().replace(" ","_")
        self.sequence = sequence.upper()
        self.rc_sequence = self.rc_seq()
    def __repr__(self) :
        #unique string representation for object
        return "DNASequence '{}'".format(self.header)
    def rc_seq(self) :
        #reverse sequence and translate
        rev_seq = self.sequence[::-1]
        rc_seq = rev_seq.translate(str.maketrans('ATGC', 'TACG'))
        return rc_seq

class ORF (object) :
    """Object to contain data for each ORF found."""
    known_base_pattern = re.compile(r'^[ATCG]*$')
    #stop codons as empty string, not represented in final peptide sequence
    codon_dict = {"AAA" : "K", "AAC" : "N", "AAG" : "K", "AAT" : "N",
	"ACA" : "T", "ACC" : "T", "ACG" : "T", "ACT" : "T",
	"AGA" : "R", "AGC" : "S", "AGG" : "R", "AGT" : "S",
	"ATA" : "I", "ATC" : "I", "ATG" : "M", "ATT" : "I",
	"CAA" : "Q", "CAC" : "H", "CAG" : "Q", "CAT" : "H",
	"CCA" : "P", "CCC" : "P", "CCG" : "P", "CCT" : "P",
	"CGA" : "R", "CGC" : "R", "CGG" : "R", "CGT" : "R",
	"CTA" : "L", "CTC" : "L", "CTG" : "L", "CTT" : "L",
	"GAA" : "E", "GAC" : "D", "GAG" : "E", "GAT" : "D",
	"GCA" : "A", "GCC" : "A", "GCG" : "A", "GCT" : "A",
	"GGA" : "G", "GGC" : "G", "GGG" : "G", "GGT" : "G",
	"GTA" : "V", "GTC" : "V", "GTG" : "V", "GTT" : "V",
	"TAA" :  "", "TAC" : "Y", "TAG" :  "", "TAT" : "Y",
	"TCA" : "S", "TCC" : "S", "TCG" : "S", "TCT" : "S",
	"TGA" :  "", "TGC" : "C", "TGG" : "W", "TGT" : "C",
	"TTA" : "L", "TTC" : "F", "TTG" : "L", "TTT" : "F"}
    def __init__(self, source, sequence, position, strand) :
        self.source = source
        self.sequence = sequence.upper()
        self.position = position
        self.strand = strand
        self.frame = self.get_frame()
        self.length = len(self.sequence)
        self.unknown = self.check_bases()
        self.protein = self.get_peptides()
    def __repr__(self) :
        #unique string representation for object
        return "ORF {}_{}{}".format(self.source,self.strand,self.position)
    def check_bases(self) :
        #return true if contains non-ACGT bases
        check = self.known_base_pattern.match(self.sequence)
        if check is not None :
            return False
        else :
            return True
    def get_frame(self) :
        #uses start position and strand to calculate frame
        strand_number = 0
        if self.strand == 'rc' : strand_number = 3
        return str((self.position % 3) + 1 + strand_number)
    def get_peptides(self) :
        #use dict.get() to specify default value 'X' if key not found
        #builds temporary list of codons, uses this to build a list of amino acids then joins to give peptide sequence
        protein = ''.join([self.codon_dict.get(x,'X') for x in [self.sequence[i:i+3] for i in range(0,len(self.sequence), 3)]] )
        return protein
    def fasta_peptides(self) :
        #print peptide sequence in fasta format
        #split peptide sequence to 60 aas per line
        split_protein = ''
        for i in range(0,len(self.protein),60):
            if split_protein != '' : split_protein += '\n'
            split_protein += self.protein[i:i+60]
        return('{} {} {} {}\n{}'.format(self.source, self.frame, self.length, self.position, split_protein,))

def orf_scanner (strand, code) :
    """Scan through sequence in each frame codon by codon to find ORFs."""
    counter = 1
    orf = False
    temp_sequence = ''
    orf_start = 0
    #check against each frame
    for i in [0,1,2]:
        while i < len(strand):
            while orf:
                #if stop codon found, stop recording ORF sequence and save
                if strand[i:i+3] in ('TAA','TGA','TAG'):
                    orf = False
                    temp_sequence += strand[i:i+3]
                    scanner_orf_list.append(ORF(sequence.header+'_'+str(counter),temp_sequence,orf_start,code))
                    counter += 1
                    temp_sequence = ''
                #if runs into end of sequence, break inner loop and discard temp_sequence
                elif i > len(strand):
                    orf = False
                    temp_sequence = ''
                    break
                else:
                    temp_sequence += strand[i:i+3]
                    i += 3
            #if start codon found, start recording ORF sequence
            if strand[i:i+3] == 'ATG' :
                orf = True
                temp_sequence += strand[i:i+3]
                orf_start = i+1
            i += 3

def remove_overlap (orf_list) :
    """Compares all ORFs and removes shortest where sequences overlap."""
    invalid_orf_list = []
    for outer_orf in orf_list:
        for inner_orf in orf_list:
            #do not compare with self
            if outer_orf.__repr__ == inner_orf.__repr__:
                continue
            #do not compare position across strands (200-300 on forward strand does not overlap with 200-300 on complementary strand)
            if outer_orf.strand != inner_orf.strand :
                continue
            if arguments.allow_overlap and outer_orf.frame != inner_orf.frame :
                continue
            outer_start = outer_orf.position
            inner_start = inner_orf.position
            outer_end = outer_orf.position + outer_orf.length
            inner_end = inner_orf.position + inner_orf.length
            outer_length = outer_orf.length
            inner_length = inner_orf.length
            if inner_start >= outer_start and inner_start < outer_end :
                #assign the shorter ORF to variable and add to list of invalid ORFs
                shorter = outer_orf if outer_length < inner_length else inner_orf
                invalid_orf_list.append(shorter)
    return invalid_orf_list

sequence_list = []
scanner_orf_list = []
input_sequence = ''

#loop over each input file given, with IOError exception raised if file not found
for infile in arguments.filenames :
    try:
        with open(infile,'r') as myfile :
            lines = myfile.readlines()
            for line in lines :
                if '>' in line :
                    input_sequence += "\n"+line
                else :
                    input_sequence += line.rstrip()
            #search within file contents for sequences and save as DNASequence objects
            pattern = re.compile(r'(>.*)\n([A-Za-z]*)')
            matches = pattern.finditer(input_sequence)
            for m in matches:
                sequence_list.append(DNASequence(m.group(1),m.group(2)))
            #now build up list of orfs for all frames
            for sequence in sequence_list:
                orf_scanner(sequence.sequence, 'seq')
                orf_scanner(sequence.rc_sequence, 'rc')
    except IOError as error:
        print("Error - specified file '{}' not found.".format(infile))

#initial orf filtering based on command line arguments - length, nonstandard bases, allowed frames
scanner_valid_orfs = [x for x in scanner_orf_list if x.length > arguments.minlength and x.length <= arguments.maxlength and x.unknown == arguments.nonstandard and x.frame in arguments.frames]

#remove overlapping ORFs
scanner_invalid_orfs = remove_overlap(scanner_valid_orfs)
scanner_valid_count = 0
scanner_no_overlap = []
for orf in scanner_valid_orfs :
    if orf not in scanner_invalid_orfs:
        scanner_no_overlap.append(orf)
        scanner_valid_count += 1

#output results as file or to stdout
if arguments.output != '' :
    with open(arguments.output,'w') as outfile :
        #write to file - use object's formatted output funct
        for orf in scanner_valid_orfs :
            outfile.write(orf.fasta_peptides()+'\n')
    print("Output saved to file: {}".format(arguments.output))
else :
    print("Valid protein ORFs found:")
    #write to stdout - use output funct
    for orf in scanner_valid_orfs :
        print(orf.fasta_peptides())
        
#only runs if summary selected on command line
if arguments.summary == True :
    print("#Summary:")
    print("#Files read in: {}".format(str(len(arguments.filenames))))
    print("#Sequences found: {}".format(str(len(sequence_list))))
    print("#Total ORFs found: {}".format(str(len(scanner_orf_list))))
    print("#Valid ORFs found: {}".format(str(scanner_valid_count)))
    #dict comprehension for dict of frame # and ORFs in that frame, internal list compr for list of orfs in each frame
    scanner_framecounts = {str(i) : len([x for x in scanner_no_overlap if x.frame == str(i)]) for i in range(1,7)}
    print('#ORF counts by frame:')
    for i in range(1,7):
        print("#{}: {}" .format(str(i), scanner_framecounts[str(i)]))