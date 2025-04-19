from abc import ABC
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import numpy as np


class BiologicalSequence(ABC):
    def __init__(self, sequence: str):
        self.seq = sequence.upper()

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, index):
        return self.seq[index]

    def __str__(self):
        return self.seq

    def __repr__(self):
        return f"{self.__class__.__name__}({self.seq})"

    _alphabet = set()

    def check_sequence(self):
        return set(list(self.seq)).issubset(self._alphabet)


class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, sequence):
        super().check_sequence
        super().__init__(sequence)
        if not super().check_sequence():
            raise ValueError(f"Invalid sequence: {sequence}")

    _complement_rule = {}

    _alphabet = set("ACGTU")

    def complement(self):
        return self.__class__("".join(self._complement_rule[nucl] for nucl in self.seq))

    def reverse(self):
        return self.seq[::-1]

    def reverse_complement(self):
        return self.__class__(self.complement().seq[::-1])


class DNASequence(NucleicAcidSequence):
    def transcribe(self):
        return RNASequence(self.seq.replace("T", "U").replace("t", "u"))
    _complement_rule = {"A": "T", "G": "C", "C": "G", "T": "A"}
    _alphabet = set(_complement_rule.keys())


class RNASequence(NucleicAcidSequence):
    _complement_rule = {"A": "U", "G": "C",
                        "C": "G", "U": "A"}
    _alphabet = {"A", "U", "G", "C"}

class AminoAcidSequence(BiologicalSequence):
    _alphabet = set(list("ARNDCEQGHILKMFPSTWYV"))

    _charges = {
        "R": +1, "K": +1, "H": +0.5,
        "D": -1, "E": -1
    }

    def charge(self):
        return sum(self._charges.get(aa, 0) for aa in self.seq)


# fastq filrator
def make_bounds(limit: float | int | tuple) -> tuple:
    """Define limits, if it
    specified with one case"""
    if isinstance(limit, (int, float)):
        return (0, limit)
    return limit


# Function for filter fastq file format
def FilterFastQC(
        fastq_input: str,
        fastq_ouput: str,
        gc_bounds=(0, 100),
        len_bounds=(0, 2**32),
        quality_threshold=0) -> None:
    gc_bounds, length_bounds = make_bounds(gc_bounds), make_bounds(len_bounds)
    with open(fastq_ouput, 'w') as out_handle:
        for record in SeqIO.parse(fastq_input, 'fastq'):
            seq_len = len(record.seq)
            avg_quality = np.mean(record.letter_annotations['phred_quality'])
            gc_content = gc_fraction(record.seq)
            if gc_bounds[0] <= gc_content*100 <= gc_bounds[1] and \
               length_bounds[0] <= seq_len <= length_bounds[1] and \
               avg_quality >= quality_threshold:
                SeqIO.write(record, out_handle, 'fastq')
