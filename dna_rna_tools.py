#!/usr/bin/env python
# Словарь для комплементарности
COMPL_DNA = {"A": "T", "G": "C", "C": "G", "T": "A",
            "a": "t", "g": "c", "c": "g", "t": "a"}
COMPL_RNA = {"A": "U", "G": "C", "C": "G", "U": "A",
            "a": "u", "g": "c", "c": "g", "u": "a"}
stop_codons = ["UAA", "UAG", "UGA"]
start_codon = 'AUG'

def complement(seq: str) -> str:  # Возвращает комплементаруню цепь
    seq_result = ''
    compl_rule = COMPL_RNA if check_acid_type(seq) == "RNA" else COMPL_DNA
    for nucl in seq:
        seq_result += compl_rule[nucl]
    return seq_result


def transcribe(seq: str) -> str:  # Транскрипция ДНК, возрващает РНК
    if check_acid_type(seq) == 'DNA':
        return seq.replace('T', "U").replace('t', 'u')
    raise ValueError("cant transcribe 'U'racil")


def reverse(seq: str):  # Разворачивает нуклеотидную цепь
    return seq[::-1]


def reverse_complement(seq: str):  # Разворачивает коплементарную нуклеотидную цепль
    return reverse(complement(seq))


# Поиск возможной РС, если РС есть возвращает координаты,
# если нет возвращает ноль
def find_possible_ORF(seq: str, start=start_codon, stop=stop_codons) -> list|int:
    if check_acid_type(seq) == 'RNA':
        seq = seq.upper()
        stop_list = [seq.find(codon)
                     for codon in stop_codons]
        stop_list = [i for i in stop_list if i != -1]
        start = seq.find(start_codon)
        if len(stop_list) == 0:
            return 0
        start = seq.find(start_codon)
        stop = min(stop_list)
        if (start != -1 and stop != -1):
            if ((stop - start)**2)**0.5 % 3 == 0:
                return start, stop    # Начало, Конец РС
            return 0
        return 0
    return find_possible_ORF(transcribe(seq))


def GC_status(seq: str) -> float:  # ГЦ состав
    seq = seq.upper() # type: ignore
    result = (seq.count('G') + seq.count('C'))/len(seq)
    return result*100


def check_acid_type(seq: str) -> str:  # Чек на валидный тип НК
    seq_up = seq.upper()
    if set(seq_up).issubset(('A', 'U', 'G', 'C')):
        return "RNA"
    if set(seq_up).issubset(('A', 'T', 'G', 'C')):
        return "DNA"
    raise ValueError("Incorrect sequence")


