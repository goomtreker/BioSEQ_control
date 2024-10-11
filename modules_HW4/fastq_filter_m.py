# Импорт модулей
from os import mkdir
from os.path import isdir
from os.path import join as path_join
from BioSEQ import filter_fastq


def phread_score(qual: str) -> float:  # Функция для расчета phread score
    return sum(ord(q) - 33 for q in qual) / len(qual)


def make_bounds(limit):
    """Определение лимитов,
    в том случае если подали
    одним числом"""
    if isinstance(limit, (int, float)):
        return (0, limit)
    return limit


def transform_to_dict(fastq_file):
        """
        Функци для превразщения 4 строк к словорю
        """
        seq_read = fastq_file.readline
        seq_dict = {
                seq_read().strip():
                (seq_read().strip(),
                seq_read().strip(),
                seq_read().strip())[::2]
                    }
        return seq_dict


# Функция для фильтрации fastq файла
def Record_filt_fastq(
        fastq_input,
        fastq_ouput,
        gc_bounds=(0, 100),
        len_bounds=(0, 2**32),
        quality_threshold=0):
    if isdir('filtered'):
        with open(fastq_input, 'r') as input, open(path_join('filtered', fastq_ouput), 'w') as output:
                while True:
                    dict = transform_to_dict(input)
                    if '' not in dict.keys():
                        result_dict = filter_fastq(dict, gc_bounds, len_bounds, quality_threshold)
                        for key, values in result_dict.items():
                            output.write(key + '\n' + values[0] + '\n')
                            output.write(key.replace('@', '+') + '\n' + values[1] + '\n')
                    else:
                        return
    mkdir('filtered')
    return Record_filt_fastq(fastq_input, fastq_ouput, gc_bounds, len_bounds, quality_threshold)
