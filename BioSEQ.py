from abc import ABC, abstractmethod

class BiologicalSequence(ABC):
    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __repr__(self):
        pass

    @abstractmethod
    def __getitem__(self, index):
        pass

    def check_sequence(self):
        pass

class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, sequence):
        super().__init__(sequence)

    def complement(self):
