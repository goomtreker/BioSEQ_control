from abc import ABC


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

    _complement_rule = {"A": "T", "G": "C",
                        "C": "G", "T": "A"}

    _alphabet = _complement_rule.keys()

    def complement(self):
        return DNASequence("".join(self._complement_rule[nucl] for nucl in self.seq))

    def reverse(self):
        return self.seq[::-1]

    def reverse_complement(self):
        return DNASequence(self.complement()[::-1])


class DNASequence(NucleicAcidSequence):
    def transcribe(self):
        return RNASequence(self.seq.replace("T", "U").replace("t", "u"))


class RNASequence(NucleicAcidSequence):
    _alphabet = {"A", "U", "G", "C"}


class AminoAcidSequence(BiologicalSequence):
    _alphabet = set(list("ARNDCEQGHILKMFPSTWYV"))

    _charges = {
        "R": +1, "K": +1, "H": +0.5,
        "D": -1, "E": -1
    }

    def charge(self):
        return sum(self._charges.get(aa, 0) for aa in self.seq)
