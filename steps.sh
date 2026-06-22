#!/bin/sh

python3 create_paths.py 
## basedir, mains, lib, s4pred are given
## creates the subfolders in a new basedir folder (current work directory) # takes no inputs

# mains  -lib
#        - s4pred

#        - validated_fasta/
#        - chem/
#        - dom/
#        - sec_struc/
#        - aberr_fasta/   -metadata (optional created by subscript later)

######################################

python3 ./mains/validate_fasta.py
## takes a transcript and protein reference file and validates the sequences and trasforms headers in a unified format
## outfile written into  basedir/validated_fasta/ folder

#usage: main_get_validated_fasta.py [-h] -t TRANSCRIPTREF -p PROTEINREF -o OUTFILENAME [-c] [-aao]
#
#options:
#  -h, --help            show this help message and exit
#  -t TRANSCRIPTREF, --transcriptREF TRANSCRIPTREF
#                        transcript reference file (fasta format)
#  -p PROTEINREF, --proteinREF PROTEINREF
#                        protein reference file (fasta format)
#  -o OUTFILENAME, --outfilename OUTFILENAME
#                        name of the outfile to create (ending .fa or .fasta)
#  -c, --CanonicalFilter
#                        Only use longest transript per gene
#  -aao, --aminoacid_outlevel
#                        Level of your output file sequence default is mRNA. Option changes
#                        fasta output format to aminoacid/peptide sequence. Please be aware
#                        that this level can not be used for frameshifting. Other aberrant
#                        translation options are possible though.

######################################

python3 ./main_aberrant_transls_pep.py
## input from basedir/validated_fasta/
## outfile written into  basedir/aberr_fasta/ folder
## metadataoutfile written into  basedir/aberr_fasta/metadata folder

#usage: main_aberrant_transls_pep.py [-h] {frameshift,repeat,skip_deletion,truncation,alternative_start,insertion,reversion,mutation,vcf_file}

#positional arguments:
#  {frameshift,repeat,skip_deletion,truncation,alternative_start,insertion,reversion,mutation,vcf_file}
#    frameshift          mRNA or Amino Acid sequence should be frameshifted
#    repeat              mRNA or Amino Acid sequence should contain repeats
#    skip_deletion       mRNA or Amino Acid sequence should contain a skip or a deletion
#    truncation          mRNA or Amino Acid sequence should be truncated
#    alternative_start   mRNA or Amino Acid sequence should have an alternative start site
#    insertion           mRNA or Amino Acid sequence should have an insertion
#    reversion           mRNA or Amino Acid sequence should contain a reverted piece
#    mutation            mRNA mutation should be translated
#    vcf_file            vcf file with mutations should be translated

#usage: main_aberrant_transls_pep.py frameshift [-h] -fsd FRAMESHIFT_DIRECTION [-rt] -i INPUT -pst
#                                               POSITION_STRING_TYPE [-ff] [-ca] -ps POSITION_STRING -mca
#                                               MRNA_CODON_AMINOACID [-ref REFERENCE] [-nb NRBEFORE] [-na NRAFTER]
#                                               [-sl SEQUENCELENGTH] [-trm TRIMM] [-trp] -o OUTFILENAME [-aao]
#                                               [-chr] [-ct CONTROLTYPE] [-cn CONTROLNUMBER] [-meta]

#options:
#  -h, --help            show this help message and exit
#  -fsd FRAMESHIFT_DIRECTION, --frameshift_direction FRAMESHIFT_DIRECTION 
#                        In which direction should the frameshift be considered? [p1],[m1],[both])
#  -rt, --remove_truncation
#                        Flag whether frameshifts resulting in immediate stop codons should be kept (default) or be
#                        removed (if -rt).
#  -i INPUT, --input INPUT
#                        Name of input. Either a fasta file [file.fasta], tab separated pos file (ENSTID/GeneID pos
#                        or chr:1-2) [file.tsv], a mutation file [file.vcf] or a single input sequence (mRNA or
#                        amino acid level) [ACGTGGAGACGAT]
#  -pst POSITION_STRING_TYPE, --position_string_type POSITION_STRING_TYPE
#                        The type of position where aberrant translation starts. options are ([codon], [aminoacid],
#                        [sequence], [startATGpos] (from start codon), [STOPpos] (from stop codon), [FILEpos] (from
#                        given tsv file or vcf file), [regex] regular expression in quotation marks ''
#  -ff, --force_frame    In case you gave a sequence as -pst argument, this will only consider occurrences were the
#                        sequence was found to be in the normal reading frame
#  -ca, --codon_aware    If you chose the -spt aminoacid option but still want the output sequences to be aware and
#                        include the different codons
#  -ps POSITION_STRING, --position_string POSITION_STRING
#                        The actual position where aberrance starts? (e.g [TTA] for codon, aminoacid [L], sequence
#                        [CTGGTGATTGATGCG], a position [9], based on input file [FILE])
#  -mca MRNA_CODON_AMINOACID, --mRNA_codon_aminoacid MRNA_CODON_AMINOACID
#                        is the event to be on mRNA [mRNA], codon [codon] or amino acid [aminoacid] level. This
#                        determines also the distances and lengths of other options to be counted either in
#                        nucleotides, or codons/amino acids! If you provide any vcf or chr based input mRNA level is
#                        required.
#  -ref REFERENCE, --Reference REFERENCE
#                        REquired if running with a positionbased file. REFERENCE on mRNA or peptide (based on -mca
#                        choice) level to map the positions and mutations against. Given as file [reference.fa] or
#                        format [hg19,hg38].
#  -nb NRBEFORE, --NRbefore NRBEFORE
#                        Minimum distance before the aberrant translation site [100]
#  -na NRAFTER, --NRafter NRAFTER
#                        Minimum distance between aberrant translation site and canonical stop codon [50]
#  -sl SEQUENCELENGTH, --sequencelength SEQUENCELENGTH
#                        Minimum length of the created aberrant peptide after the aberrant translation site [100]
#  -trm TRIMM, --trimm TRIMM
#                        If true all sequences are trimmed to length of NR before and NR after, this can be used for
#                        oligo creation. The trimming is centered around the site of aberrant translation. NOT
#                        around the later aberrant sequence. please take this into account. The returned trimmed
#                        sequence is given on the level of the position string type.
#  -trp, --trypsin       If this flag is used output amino acid sequences will be trypsinized into peptides at K/R.
#                        Peptide length is 5-55 amino acid
#  -o OUTFILENAME, --outfilename OUTFILENAME
#                        name of the outfile
#  -aao, --aminoacid_outlevel
#                        Level of your output file sequence default is aminoacid. Option changes fasta output format
#                        to mRNA sequence. Please be aware that this level can not be used in further pipeline
#  -chr                  If this flag is used the positions in the input tsv are recognized as chromosome based
#  -ct CONTROLTYPE, --controltype CONTROLTYPE
#                        Extra control type to be included, default random positions are always included. options
#                        are [controlscramble], [controlreverse], [controlrandom], [controlcharged], [all]
#  -cn CONTROLNUMBER, --controlnumber CONTROLNUMBER
#                        Number of controls per type to be included
#  -meta, --metadata     Flag to turn on metadata. Will write a metadata file of the aberrant fasta file

######################################

python3 define_chemproperties.py
## input from basedir/aberr_fasta/
## outfile written into  basedir/chem/ folder
## uses the aberrant fasta file and calcuates the differential chemical properties between aberrant and wt sequence 

#usage: main_get_chemproperties.py [-h] -i INPUT -o OUTFILENAME [-p PROPERTIES]
#
#options:
#  -h, --help            show this help message and exit
#  -i INPUT, --input INPUT
#                        inputfile (fasta format)
#  -o OUTFILENAME, --outfilename OUTFILENAME
#                        name of the outfile to create (ending .fa or .fasta
#  -p PROPERTIES, --properties PROPERTIES
#                        names of the properties to use (comma , separated without spaces
#                        inbetween e.g boman,charge,all)

######################################

python3 separate_fasta.py

## input from basedir/aberr_fasta/
## filtered outfile fasta written into  basedir/aberr_fasta/
## separates filtered_fasta into separate fastas in basedir/dom/ folder with ending batch1.fa 
## batch fasta files are then input for interproscan.sh

#usage: separate_fasta.py [-h] [-n NUMBER_ENTRIES] -i INPUTFASTAFILE -ipp INTERPROSCANPATH
#
#options:
#  -h, --help            show this help message and exit
#  -n NUMBER_ENTRIES, --number_entries NUMBER_ENTRIES
#                        number of entries to split fasta file into. Number of lines / 2
#                        (>header\sequence)
#  -i INPUTFASTAFILE, --inputfastafile INPUTFASTAFILE
#                        inputfile aberrant fasta


######################################

for fasta in dom/*fa; do echo $fasta; bash interproscan.sh; done



######################################


python3 form_secondary_structures.py

## input from basedir/aberr_fasta/
## outfile written into  basedir/sec_struc/ folder
## uses the aberrant fasta file and calcuates the differential secondary structures between aberrant and wt sequence 
## needs s4pred installed in same directory basedir/s4pred.

#usage: form_secondary_structures.py [-h] -i INPUT -o OUTFILENAME [-of OUTPUTFORMAT]

#options:
#  -h, --help            show this help message and exit
#  -i INPUT, --input INPUT
#                        inputfile (fasta format)
#  -o OUTFILENAME, --outfilename OUTFILENAME
#                        name of the outfile to create
#  -of OUTPUTFORMAT, --outputformat OUTPUTFORMAT
#                        which format should the output be in? default = 'diff', alternatively include only
#                        'original' confidence values

######################################

python3 property_statistics.py

# takes input from previous scripts (domains, secondary structure, pdb align, chemical properties )
# output is global and per oof statistics

#usage: property_statistics.py [-h] -i INPUT -o OUTFILENAME -t
#                              {chemproperties,secondary_structure,secondary_structure_diff,pdb_tmalign,domains}

#options:
#  -h, --help            show this help message and exit
#  -i INPUT, --input INPUT
#                        inputfile
#  -o OUTFILENAME, --outfilename OUTFILENAME
#                        name of the outfile to create (automatic ending will be attached)
#  -t {chemproperties,secondary_structure,secondary_structure_diff,pdb_tmalign,domains}, --inputtype {chemproperties,secondary_structure,secondary_structure_diff,pdb_tmalign,domains}
#                        indicate the datatype of your inputfile [chemproperties],[secondary_structure],
#                        [secondary_structure_diff], [pdb_tmalign]


