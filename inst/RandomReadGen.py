import argparse
import random


def random_read_generator(chunk_size=None, length_of_reads=None):
    DNA_nucleotides = ["G", "C", "A", "T"]
    if chunk_size is not None and length_of_reads is not None:
        random_reads_gen = (random.choices(DNA_nucleotides, weights=[1, 1, 1, 1], k=length_of_reads) for i in
                            range(chunk_size))
        return random_reads_gen


def write_random_reads_in_file(number_of_reads, length_of_reads):
    if number_of_reads is not None and length_of_reads is not None:
        chunk_size = 10000
        f = open(f'random_reads_n_{number_of_reads}_l_{length_of_reads}.fastq', 'w')
        f.close()
        number_of_circles = number_of_reads / chunk_size
        num_of_circles_at_beginning = number_of_reads / chunk_size
        while number_of_circles > 0:
            print(f"write progress: {chunk_size * (num_of_circles_at_beginning - number_of_circles)}|{number_of_reads}")
            random_reads = random_read_generator(chunk_size, length_of_reads=length_of_reads)
            with open(f'random_reads_n_{number_of_reads}_l_{length_of_reads}.fastq', 'a') as file:
                for read_as_lst in random_reads:
                    read = ''.join(read_as_lst)
                    file.write('@' + '\n')
                    file.write(read + '\n')
                    file.write('+' + '\n')
                    file.write(''.join(['X'] * length_of_reads) + '\n')
            number_of_circles -= 1
        print(f"write progress: {number_of_reads}|{number_of_reads}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=int, metavar='number of generated reads')
    parser.add_argument('-l', type=int, metavar='length of reads')
    args = parser.parse_args()
    write_random_reads_in_file(number_of_reads=args.n, length_of_reads=args.l)
