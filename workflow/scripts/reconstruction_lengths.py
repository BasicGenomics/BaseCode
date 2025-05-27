import pysam
import argparse
import pandas as pd
from joblib import delayed,Parallel
from multiprocessing import Process, Manager
import warnings

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Reconstruction length statistics', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bam file')
    parser.add_argument('-o','--outfile', help='Output file')
    args = parser.parse_args()


    bamfile = args.input
    outfile = args.outfile

    res = []
    bam = pysam.AlignmentFile(bamfile, 'rb')
    for mol in bam.fetch(until_eof=True):
        qname = mol.query_name
        sample_name = '_'.join(qname.split(':')[0].split('_')[:-1])
        res.append([mol.query_name, sample_name, mol.get_tag('BC'),mol.get_tag('TC'), mol.get_tag('IC'), mol.get_tag('FC'), mol.get_tag('XT'), mol.get_tag('CC'), mol.get_tag('SC'), mol.query_length, mol.cigarstring.count('D'), mol.get_tag('F1'), mol.get_tag('T1')])
    df = pd.DataFrame(res, columns = ['MOL_NAME','SM','BC', 'TC', 'IC', 'FC', 'XT','CC','SC','QL', 'GAPS', 'F1', 'T1'])
    df['ANY_TC'] = df['TC'] > 0
    df['ANY_IC'] = df['IC'] > 0
    df['ANY_FC'] = df['FC'] > 0
    df['ANY_GAPS'] = df['GAPS'] > 0
    df.to_csv(outfile)