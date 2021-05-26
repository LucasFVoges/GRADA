import argparse


def get_length_of_reads(file, n_reads_base):
    count = 0
    sum_of_read_length = 0
    num_of_rows = n_reads_base * 4
    with open(file, "r") as file:
        while count < num_of_rows:
            line = file.readline()
            if count in range(1, num_of_rows, 4):
                sum_of_read_length += len(line.strip())
            count += 1
    avg_length = sum_of_read_length / (num_of_rows // 4)
    return avg_length


def write_avg_length_in_file(file, n_reads_base):
    avg_length = get_length_of_reads(file=file, n_reads_base=n_reads_base)
    with open("avg_length.txt", "w") as file:
        file.write(str(avg_length))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, metavar='fastq file')
    parser.add_argument('-nb', type=int, metavar='number of reads for the basis of calculation')
    args = parser.parse_args()
    write_avg_length_in_file(file=args.f, n_reads_base=args.nb)
