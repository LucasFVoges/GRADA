import argparse


def load_sequences(file_name):
    f = open(file_name, 'r')
    lines = f.readlines()
    f.close()
    return lines


def make_rc(seq):
    seq_rc_lst = []

    seq_reversed = seq[::-1]

    for nuc in seq_reversed:
        if nuc == "C":
            seq_rc_lst.append("G")
        if nuc == "T":
            seq_rc_lst.append("A")
        if nuc == "G":
            seq_rc_lst.append("C")
        if nuc == "A":
            seq_rc_lst.append("T")
        if nuc == "N":
            seq_rc_lst.append("N")

    seq_rc = ''.join(seq_rc_lst)

    return seq_rc


def make_tuples_out_of_name_and_seq_from_file(lines):
    list_of_tuples = []
    for index, line in enumerate(lines):
        if line.startswith(">"):
            name = line
            seq = lines[index + 1]
            list_of_tuples.append((name, seq))
    return list_of_tuples


def make_tuples_of_name_and_seq_in_rc_version(list_of_tuples):
    list_of_tuples_rc = []
    for name, seq in list_of_tuples:
        name_rc = '_'.join([name.strip(), 'rc'])
        seq_rc = ''.join(make_rc(seq=seq))
        list_of_tuples_rc.append((name.strip(), seq.strip()))
        list_of_tuples_rc.append((name_rc, seq_rc))
    return list_of_tuples_rc


def write_new_file_with_rc_content(list_of_tuples_rc):
    f = open('seq_rc.txt', 'w')
    for name, seq in list_of_tuples_rc:
        f.write(name + '\n')
        f.write(seq + '\n')
    f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', type=str, metavar='path to file')
    args = parser.parse_args()
    lines = load_sequences(file_name=args.f)
    list_of_tupls = make_tuples_out_of_name_and_seq_from_file(lines=lines)
    list_of_tpls_rc = make_tuples_of_name_and_seq_in_rc_version(list_of_tuples=list_of_tupls)
    write_new_file_with_rc_content(list_of_tuples_rc=list_of_tpls_rc)
