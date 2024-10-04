def mean(array: list[int|float]) -> float:
    return sum(array)/len(array)


def phread_score(field_4: str) -> float:  #Функция для расчета 
    symbol_collector = []
    for symbol in field_4:
        symbol_collector.append(ord(symbol) - 33)
    return mean(symbol_collector)


def check_limits(gc_bounds: tuple|int, length_bounds: tuple|int) -> tuple:
    if type(gc_bounds) is int:
        gc_bounds = (0, gc_bounds)
    elif type(length_bounds) is int:
        length_bounds = (0, length_bounds)
    return gc_bounds, length_bounds


def filter_fastq(seqs, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0):
    filtered_seqs = {}
    for name, (seq, field_4) in seqs:
        length, quality, gc = len(seq), phread_score(quality), GC_status(seq)
        gc_bounds, length_bounds = check_limits(gc_bounds, length_bounds)
        filtration =  (gc in range(gc_bounds[0], gc_bounds[1] + 1),
                          length in range(length_bounds[0], length_bounds[1] + 1),
                          quality_threshold < quality)
        if filtration == (True, True, True):
            filtered_seqs[name] = (seq, field_4)
        continue