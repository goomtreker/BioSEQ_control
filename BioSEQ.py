# Импорт модулей
from modules_HW4 import dna_rna_tools as drt
from modules_HW4 import fastq_filter_m as fq


# Главная функция для модуля fastq_filter_m
def filter_fastq(
        seqs: dict[str: tuple[str, str]],
        gc_bounds=(0, 100),
        length_bounds=(0, 2**32), quality_threshold=0):
    '''
    The main function of fastq_filter_m, filter the sequecnes,
    according to specified parameters: 
    gc_bounds, length_bounds, quality threshold.
    For input it takes dictionary : {sequence name: (sequence, quality)}
    return same dictionary after filtration
    '''
    filtered_seqs = {}
    for name, (seq, qual_str) in seqs.items():
        length, quality, gc = len(seq), fq.phread_score(qual_str), drt.GC_status(seq)
        gc_bounds, length_bounds = fq.make_bounds(gc_bounds), fq.make_bounds(length_bounds)
        filtration = (
            gc_bounds[0] <= gc <= gc_bounds[1],
            length in range(length_bounds[0], length_bounds[1] + 1),
            quality_threshold < quality)
        if all(filtration):
            filtered_seqs[name] = (seq, qual_str)
    return filtered_seqs




# Главная функция для модуля dna_rna_tools,
def run_dna_rna_tools(*args: list[str]) -> list | str:
    """
    For inputs take list of sequences, last argument must be
    one of the dna_rna_tools_function
    """
    *seqs, action = args
    action = eval("drt." + action)
    collector = []
    for seq in seqs:
        drt.check_acid_type(seq)
        collector.append(action(seq))
    if len(collector) == 1:
        return collector[0]
    return collector
