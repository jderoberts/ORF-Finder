# ORF-Finder

This is a Python programming project completed as part of my Bioinformatics MSc at the University of Manchester. The command line tool searches an input genetic sequence for Open Reading Frames (potential protein-coding sequences).

## Motivation

DNA sequencing outputs a list of As, Ts, Cs and Gs that make up the genetic code.  To predict what proteins are encoded by a given sequence we need to find the Open Reading Frames.  In eukaryotic organisms (i.e. not bacteria) these match the pattern:  
ATG([ATGC]{3})+(TAA|TGA|TAG)  
A universal start sequence ATG, then one or more codons (triplets of bases), followed by a termination signal.  Each codon is matched to a specific amino acid, the subunits that make up proteins - by mapping the codons in the located ORF to their corresponding amino acids we can calculate the encoded protein sequence.

## Assignment
A single Python script to detect ORFs and predict protein sequences as the first part of a pipeline of scripts to predict proteins for virtual mass spectrometry.

## Input: 
 - DNA sequence(s) in .fasta format.  Tested with C. elegans genome.
## Output: 
 - All protein sequences potentially encoded by the input DNA, written to standard output or .fasta file.

## Documentation:
orf_finder.py reads in sequence files in standard .fasta format, splits files into individual  
sequences and locates all open reading frames (ORFs) in each sequence.  Multiple .fasta files  
can be specified, and each file can contain multiple sequences.  Each sequence is scanned for  
ORFs matching the pattern (start codon, x codons, stop codon).  Default behaviour is to only  
return ORFs between 150 and 2500nt in length, from all reading frames, containing no nonstandard  
bases and with no overlapping ORFs.  This behaviour can be overridden using the appropriate  
optional flags, e.g. -f 123 if only ORFs on the forward strand are of interest.  Valid ORFs  
are assigned a unique identifier and converted to a peptide sequence using the standard genetic  
code; if nonstandard bases are permitted the resulting unknown codons will be marked 'X'.  
All output in fasta format is returned to standard output by default; if a filename is specified  
the output will instead be written to file (e.g. -o output.fasta).  
Summary information of number of sequences, total number of ORFs, ORFs per reading frame etc.  
can be written to StdOut by using the -s flag.

## Learning outcomes:
 - Command line Python using argparse
 - Practice with regular expressions
 - Object orientation - practical use of classes

