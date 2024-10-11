# import modules
from os import mkdir
from os.path import isdir
from os.path import join as path_join
from BioSEQ import filter_fastq


def phread_score(qual: str) -> float:  # estimate the read quality
    return sum(ord(q) - 33 for q in qual) / len(qual)


def make_bounds(limit: float|int|tuple) -> float:
    """Define limits, if it 
    specified with one case"""
    if isinstance(limit, (int, float)):
        return (0, limit)
    return limit


def transform_to_dict(fastq_file: __file__) -> dict:
        """
        Transform 4 readline to dict
        for filter_fastq
        """
        seq_read = fastq_file.readline
        seq_dict = {
                seq_read().strip():
                (seq_read().strip(),
                seq_read().strip(),
                seq_read().strip())[::2]
                    }
        return seq_dict


# Function for filter fastq file format
def Record_filt_fastq(
        fastq_input,
        fastq_ouput,
        gc_bounds=(0, 100),
        len_bounds=(0, 2**32),
        quality_threshold=0) -> __file__:
    """
    Function for record filt fastq_file, for input takes file fastq format
    for output return file with same format, with filtered seqs
    """
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
