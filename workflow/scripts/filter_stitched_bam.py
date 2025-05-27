import pysam
import argparse

def trim_mol(mol):
    while True:
        if mol.cigartuples[0][0] == 3:
            mol.cigartuples = mol.cigartuples[1:]
            continue
        if mol.cigartuples[-1][0] == 3:
            mol.cigartuples = mol.cigartuples[:-1]
            continue
        if len(mol.cigartuples) <= 2:
            return mol
        if mol.is_reverse:
            if mol.cigartuples[0][0] == 0 and mol.cigartuples[1][0] == 3 and mol.cigartuples[0][1] <= 25:
                q = mol.query_qualities
                mol.query_sequence = mol.query_sequence[mol.cigartuples[0][1]:]
                mol.query_qualities = q[mol.cigartuples[0][1]:]
                mol.reference_start += (mol.cigartuples[0][1]+mol.cigartuples[1][1])
                mol.cigartuples = mol.cigartuples[2:]
            elif mol.cigartuples[-1][0] == 0 and mol.cigartuples[-2][0] == 3 and mol.cigartuples[-1][1] <= 25:
                q = mol.query_qualities
                mol.query_sequence = mol.query_sequence[:-mol.cigartuples[-1][1]]
                mol.query_qualities = q[:-mol.cigartuples[-1][1]]
                mol.cigartuples = mol.cigartuples[:-2]
            elif mol.cigartuples[0][0] == 0 and mol.cigartuples[1][0] == 2 and mol.cigartuples[0][1] <= 25:
                q = mol.query_qualities
                mol.query_sequence = mol.query_sequence[mol.cigartuples[0][1]:]
                mol.query_qualities = q[mol.cigartuples[0][1]:]
                mol.reference_start += (mol.cigartuples[0][1]+mol.cigartuples[1][1])
                mol.cigartuples = mol.cigartuples[2:]
            elif mol.cigartuples[-1][0] == 0 and mol.cigartuples[-2][0] == 2 and mol.cigartuples[-1][1] <= 25:
                q = mol.query_qualities
                mol.query_sequence = mol.query_sequence[:-mol.cigartuples[-1][1]]
                mol.query_qualities = q[:-mol.cigartuples[-1][1]]
                mol.cigartuples = mol.cigartuples[:-2]
            else:
                return mol
        else:
            if mol.cigartuples[-1][0] == 0 and mol.cigartuples[-2][0] == 3 and mol.cigartuples[-1][1] <= 25:
                q = mol.query_qualities
                mol.query_sequence = mol.query_sequence[:-mol.cigartuples[-1][1]]
                mol.query_qualities = q[:-mol.cigartuples[-1][1]]
                mol.cigartuples = mol.cigartuples[:-2]
            elif mol.cigartuples[0][0] == 0 and mol.cigartuples[1][0] == 3 and mol.cigartuples[0][1] <= 25:
                q = mol.query_qualities
                mol.query_sequence = mol.query_sequence[mol.cigartuples[0][1]:]
                mol.query_qualities = q[mol.cigartuples[0][1]:]
                mol.reference_start += (mol.cigartuples[0][1]+mol.cigartuples[1][1])
                mol.cigartuples = mol.cigartuples[2:]
            elif mol.cigartuples[-1][0] == 0 and mol.cigartuples[-2][0] == 2 and mol.cigartuples[-1][1] <= 25:
                q = mol.query_qualities
                mol.query_sequence = mol.query_sequence[:-mol.cigartuples[-1][1]]
                mol.query_qualities = q[:-mol.cigartuples[-1][1]]
                mol.cigartuples = mol.cigartuples[:-2]
            elif mol.cigartuples[0][0] == 0 and mol.cigartuples[1][0] == 2 and mol.cigartuples[0][1] <= 25:
                q = mol.query_qualities
                mol.query_sequence = mol.query_sequence[mol.cigartuples[0][1]:]
                mol.query_qualities = q[mol.cigartuples[0][1]:]
                mol.reference_start += (mol.cigartuples[0][1]+mol.cigartuples[1][1])
                mol.cigartuples = mol.cigartuples[2:]
            else:
                return mol
def fix_cigar(mol):
    new_tuples = []
    prev_m = False
    new_len = 0
    for (op, length) in mol.cigartuples:
        if op == 0:
            if prev_m:
                new_len += length
            else:
                new_len = length
            prev_m = True
        else:
            if prev_m:
                new_tuples.append((0, new_len))
            new_tuples.append((op, length))
            prev_m = False
    if prev_m:
        new_tuples.append((0, new_len))
    mol.cigartuples = new_tuples
    return mol
def hamming_distance(sub_s1, sub_s2):
    return sum(c1 != c2 for c1, c2 in zip(sub_s1, sub_s2))
complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
def reverse_complement(dna_sequence):

    reversed_sequence = dna_sequence[::-1]
    reverse_complement_sequence = ''.join(complement[base] for base in reversed_sequence)
    
    return reverse_complement_sequence

    

def main():
    parser = argparse.ArgumentParser(description='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input stitched bam file')
    parser.add_argument('--molecules-out', type=str, help='Molecules out bam file')
    args = parser.parse_args()

    bamfile = args.input
    molecules_out_bamfile = args.molecules_out


    bam_in = pysam.AlignmentFile(bamfile, 'rb')

    bam_out = pysam.AlignmentFile(molecules_out_bamfile, 'wb', template=bam_in)

    for mol in bam_in.fetch(until_eof=True):
        if mol.query_name.split(':')[-1][0] == '_' or mol.query_name.split(':')[-1][:5] == 'Unass':
            continue
        mol = trim_mol(mol)
        mol = fix_cigar(mol)

        bam_out.write(mol)
    
    bam_in.close()
    bam_out.close()


if __name__ == '__main__':
    main()