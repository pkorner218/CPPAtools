# CPPAtools

****
# A computation toolset for Chimera Protein Prediction and Analysis 

please note that the lib folders should be decompressed (it is .tar compressed) in place before starting pipeline

s4pred can be downloaded from 
`https://github.com/psipred/s4pred`

InterProScan can be obtained from 
`https://interproscan-docs.readthedocs.io/en/v5/HowToDownload.html`


# 1 get_validated_fasta.py
Since many translation fasta files show problematic inconsistencies to their transcript counterparts, this script ensures that both are the same.
Further problems were found with regards to annotated coding sequences not starting with ATG or ending with stop codons. The script removes transcripts with wrongly annotated CDS.
An optional filter possibility is `-C ` which keeps only the longest transcript per protein (here taken to be canonical)

usage: `python3 get_validated_fasta.py -t [transcript referencefile] -p [protein referencefile] -o [outfilename] -c [canonical] -aao [amino acid as outlevel]`

Further analysis can also be performed with a non validated file if that is preferred. 
In that case please format you headers to a format of `>ID|Gene|CDS:1-X`

# 2 main_aberrant_transls_pep.py
The main script of the pipeline. it creates the aberrant proteins based on a validated reference fasta file. Alternative give a single sequence as input. aberrant events can be chosen and consequenctly the main input is where the event should take place. first input is a level of input `mca` which can be mRNA, codon or aminoacid. Then the more detailed information of the aberrant event position and the nature of the event (how many repeats etc) can be given. for a detailed explanation please consider --help functions

usage: main_aberrant_transls_pep.py [-h] {frameshift,repeat,skip_deletion,truncation,alternative_start,insertion,reversion,mutation,vcf_file}

# 2.1 Usage for frameshifts

The frameshift direction option *-fsd* can be either *m1*, *p1* or *both*. 

`python3 main_aberrant_transls_pep.py frameshift -fsd [both|m1|p1] -i [input fastafile] -pst [codon] -ps [TGG] -mca [codon] -o [outfilename] -trp [trypsin]` 

The output fasta file further named #aberrant fasta# is the input for steps 3 and 6 and also for InterProScan.

# 3 define_chemproperties.py
uses the aberrant fasta file and calcuates the differential chemical properties between aberrant and wt sequence. This can then in the following be visualized. 

usage: `main_get_chemproperties.py [-h] -i INPUT -o OUTFILENAME [-p PROPERTIES]`

# 4 separate_fasta.py
uses the aberrant fasta file and splits it into a set of numbers. Allows to make subsequent analysis faster and reduce memory problems on certain analysis such as InterProScan.
usage: `separate_fasta.py [-h] [-n NUMBER_ENTRIES] -i INPUTFASTAFILE -ipp INTERPROSCANPATH`

# 5  count_domains.py 
counts domains provided by InterProscan. It references them back to a given aberrant input fastafile and compares domain numbers of aberrant and WT proteins.  
usage: `count_domains.py [-h] -idf INPUTDOMAINFILE -iff INPUTFASTAFILE -o OUTFILENAME`

# 6 form_secondary_structures.py
calculates secondary structure differences for a given aberrant fasta file using s4pred. these can then be visualized and be further analyzed.
#usage: form_secondary_structures.py [-h] -i INPUT -o OUTFILENAME [-of OUTPUTFORMAT]
