# Импорт модулей
import dna_rna_tools as drt
import fastq_filter_m as fq


# Главная функция для модуля fastq_filter_m
def filter_fastq(
        seqs: dict[str | tuple[str | str]],
        gc_bounds=(0, 100),
        length_bounds=(0, 2**32), quality_threshold=0):
    '''
    Главная функция для модуля fastq_filter,
    Фильтрует словарь из {Имя сиквенса:(сиквенс;строка с символами phread)}
    На вход получает словарь и пороги для фильтрации
    выходе получает словарь,
    '''
    filtered_seqs = {}
    for name, (seq, field_4) in seqs.items():
        length, quality, gc = len(seq), fq.phread_score(field_4), drt.GC_status(seq)
        gc_bounds, length_bounds = fq.check_limits(gc_bounds, length_bounds)
        filtration = (
            gc_bounds[0] <= gc <= gc_bounds[1],
            length in range(length_bounds[0], length_bounds[1] + 1),
            quality_threshold < quality)
        if filtration == (True, True, True):
            filtered_seqs[name] = (seq, field_4)
        continue
    return filtered_seqs


# Главная функция для модуля dna_rna_tools,
def run_dna_rna_tools(*args: list[str]) -> list|str:
    """
    На вход получает агрументы листом, вначале должны быть
    сиквенсы, последней командой
    прописывается манипуляция с сиквенсами
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