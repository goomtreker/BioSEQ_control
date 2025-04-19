def convert_multiline_fasta_to_oneline(path_to_fasta: str) -> None:
    with open(path_to_fasta, 'r') as input_fastq, open(path_to_fasta.split('/')[-1] + '_oneline', 'w') as out:
        flag = False
        for lines in input_fastq:
            if '>' in lines and flag is False:
                out.write(lines)
                flag = True
            elif '>' in lines and flag:
                out.write('\n' + lines)
            else:
                out.write(lines.strip())


def parse_blast_output(input_fasta: str, output_fasta: str) -> None:
    with open(input_fasta, 'r') as input_data, open(output_fasta, 'w') as output:
        read_next_line = input_data.readline
        list_of_description = []
        for line in input_data:
            if 'Description' in line:
                target_line = read_next_line()
                target_line = target_line.split('  ')[0]
                list_of_description.append(target_line.split('[')[0])
        for descr in sorted(list_of_description):
            output.write(descr + '\n')