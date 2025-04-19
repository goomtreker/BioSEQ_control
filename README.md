# BioSeq_control

## О BioSeq_control
***Bioseq_control*** is small python program designed for parsing and processing sequences.
The softwere has two python modules and main script.


- **BioSEQ.py**


### Installing

```bash
git clone git@github.com:goomtreker/BioSEQ_control.git
```
Adding modules:
Run python interpretator that you prefer
```python
import sys
sys.path.append("path/to/BioSEQ_control")
```
### Description
Features \

Classes for handling biological sequences: 

    BiologicalSequence — Base abstract class.

    NucleicAcidSequence — Class for working with nucleic acids (DNA/RNA).

    DNASequence — Class for DNA sequences, supports transcription to RNA.

    RNASequence — Class for RNA sequences.

    AminoAcidSequence — Class for amino acid sequences, supports charge calculation.

FASTQ file filtering:

The FilterFastQC function allows filtering records in a FASTQ file based on:

    Sequence length

    GC content

    Average Phred quality score



```

```filter_fastq``` input accept a dictionary {Sequence name: (read, read quality symbols)}. Example:
```
EXAMPLE_FASTQ = {
    '@SRX079801': ('ACAGCAACATAAACATGA....AA', 'FGGGFGGGFGGGFG...GGD'),
    '@SRX079802': ('ATTA...TG', 'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D'),
    ...
    '@SRX079812': ('AGT...AATGACCCG', '<98;...BBC')
    }      
```
As well as thresholds in length, gc composition, and read phread quality control.
example:
```python
filter_fastq(
        seqs: EXAMPLE_FASTQ,
        gc_bounds=44,
        length_bounds=(10000, 100000), quality_threshold=34):
```
The module now contains only one scipt

**Note**:
Questions, comments and suggestions can be asked/suggested in the issues tab.P.S. 
We often receive question:
 Do we know anything about the Bio.SeqIO package?
 Our response:
 We've never heard of it. :)
