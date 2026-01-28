from collections import deque
import multiprocessing as mp
import pandas as pd
import pysam
import numpy as np
import os
from operator import itemgetter
from pathlib import Path
_getitem0 = itemgetter(0)

def parse_gtf_gene_name(gtffile, contig, ban_set=set()):
    
    gene_list = []

    ext = Path(gtffile).suffix 
    counter = 0
    
    with open(gtffile, 'r') as f:
        for line in f:
            l = line.split('\t')
            if len(l) < 8:
                continue
            if l[2] == 'gene':
                if l[2] in ban_set:
                    continue
                if contig is not None:
                    if l[0] == contig:
                        if ext == '.gff3':
                            try:
                                gene_list.append({'gene_id': l[8].split(';')[3].split('=')[-1], 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4]), 'strand': l[6]})
                            except:
                                gene_list.append({'gene_id': l[8].split(';')[0].split('=')[-1], 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4]), 'strand': l[6]})
                        elif ext == '.gtf':
                            try:
                                gene_list.append({'gene_id': l[8].split(' ')[5].replace('"', '').strip(';'), 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4]), 'strand': l[6]})
                            except:
                                gene_list.append({'gene_id': l[8].split(' ')[1].replace('"', '').strip(';'), 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4]), 'strand': l[6]})
                    else:
                        continue
                else:
                    if ext == '.gff3':
                        try:
                            gene_list.append({'gene_id': l[8].split(';')[3].split('=')[-1], 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4]), 'strand': l[6]})
                        except:
                            gene_list.append({'gene_id': l[8].split(';')[0].split('=')[-1], 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4]), 'strand': l[6]})
                    elif ext == '.gtf':
                        try:
                            gene_list.append({'gene_id': l[8].split(' ')[5].replace('"', '').strip(';'), 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4]), 'strand': l[6]})
                        except:
                            gene_list.append({'gene_id': l[8].split(' ')[1].replace('"', '').strip(';'), 'seqid':l[0], 'start':int(l[3]), 'end':int(l[4]), 'strand': l[6]})
            counter +=1
            if counter >10:
                break 
                
    gene_dict = {g['gene_id']: g for g in gene_list}
    
    return gene_dict

def return_UB(moldict, BC, GE, UX):
    try:
        UB = moldict[GE][BC][UX]
    except KeyError:
        UB = UX
    return(UB)
def correct_tags(inpath, threads, chr):
    global mols
    #nreads = 0
    if chr == '*':
        chrlabel = 'unmapped'
    else:
        chrlabel = chr
    outpath = inpath+".tmp."+chrlabel+".bam"
    inp = pysam.AlignmentFile(inpath, 'rb', threads = threads)
    out = pysam.AlignmentFile(outpath, 'wb', template = inp, threads = threads)
    for read in inp.fetch(chr):
        if read.get_tag('XX') == 'threep_BCUMI_read':
            umi = read.get_tag('UB')
            umi_new = umi
            cell = get_cellBC(read)
            if read.has_tag('GE'):
                gene = read.get_tag('GE')
                umi_new = return_UB(moldict = mols, BC = cell, GE = gene, UX = umi)
            read.set_tag(tag = 'UX', value = umi, value_type = 'Z')
            read.set_tag(tag = 'UB', value = umi_new, value_type = 'Z')
        out.write(read)
    inp.close()
    out.close()
def hamming_distance(x, y):
    """Calculate the hamming distance (number of bits different) between the
    two integers given.
    >>> [hamming_distance(x, 15) for x in [0, 8, 10, 12, 14, 15]]
    [4, 3, 2, 2, 1, 0]
    """
    return bin(x ^ y).count('1')


class BKTree(object):
    """BK-tree data structure that allows fast querying of matches that are
    "close" given a function to calculate a distance metric (e.g., Hamming
    distance or Levenshtein distance).
    Each node in the tree (including the root node) is a two-tuple of
    (item, children_dict), where children_dict is a dict whose keys are
    non-negative distances of the child to the current item and whose values
    are nodes.
    """
    def __init__(self, distance_func, items=[]):
        """Initialize a BKTree instance with given distance function
        (which takes two items as parameters and returns a non-negative
        distance integer). "items" is an optional list of items to add
        on initialization.
        >>> tree = BKTree(hamming_distance)
        >>> list(tree)
        []
        >>> tree.distance_func is hamming_distance
        True
        >>> tree = BKTree(hamming_distance, [])
        >>> list(tree)
        []
        >>> tree = BKTree(hamming_distance, [0, 4, 5])
        >>> sorted(tree)
        [0, 4, 5]
        """
        self.distance_func = distance_func
        self.tree = None

        _add = self.add
        for item in items:
            _add(item)

    def add(self, item):
        """Add given item to this tree.
        >>> tree = BKTree(hamming_distance)
        >>> list(tree)
        []
        >>> tree.add(4)
        >>> sorted(tree)
        [4]
        >>> tree.add(15)
        >>> sorted(tree)
        [4, 15]
        """
        node = self.tree
        if node is None:
            self.tree = (item, {})
            return

        # Slight speed optimization -- avoid lookups inside the loop
        _distance_func = self.distance_func

        while True:
            parent, children = node
            distance = _distance_func(item, parent)
            node = children.get(distance)
            if node is None:
                children[distance] = (item, {})
                break

    def find(self, item, n):
        """Find items in this tree whose distance is less than or equal to n
        from given item, and return list of (distance, item) tuples ordered by
        distance.
        >>> tree = BKTree(hamming_distance)
        >>> tree.find(13, 1)
        []
        >>> tree.add(0)
        >>> tree.find(1, 1)
        [(1, 0)]
        >>> for item in [0, 4, 5, 14, 15]:
        ...     tree.add(item)
        >>> sorted(tree)
        [0, 0, 4, 5, 14, 15]
        >>> sorted(tree.find(13, 1))
        [(1, 5), (1, 15)]
        >>> sorted(tree.find(13, 2))
        [(1, 5), (1, 15), (2, 4), (2, 14)]
        >>> sorted(tree.find(0, 1000)) == [(hamming_distance(x, 0), x) for x in tree]
        True
        """
        if self.tree is None:
            return []

        candidates = deque([self.tree])
        found = []

        # Slight speed optimization -- avoid lookups inside the loop
        _candidates_popleft = candidates.popleft
        _candidates_extend = candidates.extend
        _found_append = found.append
        _distance_func = self.distance_func

        while candidates:
            candidate, children = _candidates_popleft()
            distance = _distance_func(candidate, item)
            if distance <= n:
                _found_append((distance, candidate))

            if children:
                lower = distance - n
                upper = distance + n
                _candidates_extend(c for d, c in children.items() if lower <= d <= upper)

        found.sort(key=_getitem0)
        return found

    def __iter__(self):
        """Return iterator over all items in this tree; items are yielded in
        arbitrary order.
        >>> tree = BKTree(hamming_distance)
        >>> list(tree)
        []
        >>> tree = BKTree(hamming_distance, [1, 2, 3, 4, 5])
        >>> sorted(tree)
        [1, 2, 3, 4, 5]
        """
        if self.tree is None:
            return

        candidates = deque([self.tree])

        # Slight speed optimization -- avoid lookups inside the loop
        _candidates_popleft = candidates.popleft
        _candidates_extend = candidates.extend

        while candidates:
            candidate, children = _candidates_popleft()
            yield candidate
            _candidates_extend(children.values())

    def __repr__(self):
        """Return a string representation of this BK-tree with a little bit of info.
        >>> BKTree(hamming_distance)
        <BKTree using hamming_distance with no top-level nodes>
        >>> BKTree(hamming_distance, [0, 4, 8, 14, 15])
        <BKTree using hamming_distance with 3 top-level nodes>
        """
        return '<{} using {} with {} top-level nodes>'.format(
            self.__class__.__name__,
            self.distance_func.__name__,
            len(self.tree[1]) if self.tree is not None else 'no',
        )

def makeDNAbits(instring):
    lookupdict = {'A': '00001', 'C': '00010', 'G': '00100', 'T': '01000', 'N': '10000'}
    outbits = '1'
    for c in instring:
        outbits += lookupdict[c]
    return int(outbits,2)

def makeDNAstring(inbits):
    lookupdict = {'00001': 'A', '00010': 'C', '00100': 'G', '01000': 'T', '10000': 'N'}
    outstring = ''
    inbitsbin = '{0:08b}'.format(inbits)
    inbitsbinlist = [inbitsbin[i:i+5] for i in range(1,len(inbitsbin),5)]
    for c in inbitsbinlist:
        outstring += lookupdict[c]
    return(outstring)

def make_unique(inlist):
    joinchar = "."
    out = []
    dict = {}
    for val in inlist:
        if val not in dict.keys():
            dict[val] = 0
        else:
            dict[val] += 1
            val = joinchar.join([val,str(dict[val])])
        out.append(val)
    return(out)

def fix_genenames(genenames, geneids):
    if '' in genenames:
        failindices = [i for i, x in enumerate(genenames) if x == '']
        for idx in failindices:
            genenames[idx] = geneids[idx]
    if len(genenames) != len(np.unique(genenames)):
        genenames = make_unique(genenames)
    return genenames

def directional_adjacency_corrections(umistrings):
    umihamdist = 1 #this could become an option
    umibits = [makeDNAbits(u) for u in np.unique(umistrings)]
    tree = BKTree(hamming_distance, umibits)
    close_hits = [tree.find(umi, umihamdist*2) for umi in umibits]
    bklist = []
    for b in close_hits:
        b_0 = b[0][1]
        b_rest = b[1:]
        bklist.extend([(makeDNAstring(b_0), int(b_1[0]/2), makeDNAstring(b_1[1])) for b_1 in b_rest])
    bkdf = pd.DataFrame(bklist)
    if len(bkdf) == 0:
        return None
    bkdf.columns=['umi1','hamdist','umi2']
    bkdf['pairing'] = [" ".join(sorted(bkdf.loc[i,['umi1','umi2']])) for i in range(len(bkdf))]
    #bkdf.drop(['umi1','umi2'], axis = 1, inplace = True)
    bkdf.drop_duplicates('pairing', inplace = True)
    #bkdf[['umi1', 'umi2']] = bkdf['pairing'].str.split(' ')
    bkdf['umi1_obs'] = [umistrings.count(x) for x in bkdf.loc[:,"umi1"]]
    bkdf['umi2_obs'] = [umistrings.count(x) for x in bkdf.loc[:,"umi2"]]
    #directional adjacency, make sure abundance is higher for more observed UMI
    bkdf = bkdf[(bkdf["umi1_obs"] <= bkdf["umi2_obs"]/2) | (bkdf["umi2_obs"] <= bkdf["umi1_obs"]/2)]
    if len(bkdf) == 0:
        return None
    bkdf['falseUMI'] = np.where((bkdf['umi1_obs']>bkdf['umi2_obs']), bkdf['umi2'], bkdf['umi1'])
    bkdf['falseUMI_obs'] = np.where((bkdf['umi1_obs']>bkdf['umi2_obs']), bkdf['umi2_obs'], bkdf['umi1_obs'])
    bkdf['trueUMI'] = np.where((bkdf['umi1_obs']<bkdf['umi2_obs']), bkdf['umi2'], bkdf['umi1'])
    bkdf['trueUMI_obs'] = np.where((bkdf['umi1_obs']<bkdf['umi2_obs']), bkdf['umi2_obs'], bkdf['umi1_obs'])
    bkdf.drop(['umi1_obs','umi1','umi2_obs','umi2','pairing'], axis = 1, inplace = True)
    bkdf = bkdf[bkdf['falseUMI'] != bkdf['trueUMI']]
    if bkdf['falseUMI'].duplicated().any(): #is the collapsing ambigous?
        bkdf = bkdf.sort_values(by='trueUMI_obs', ascending=False)  # order by most abundant true observation
        bkdf = bkdf.drop_duplicates('falseUMI')  #select the first row

    bkdf = bkdf.sort_values(by='trueUMI_obs', ascending=False)
    non_true_UMIs = bkdf[bkdf['trueUMI'].isin(bkdf['falseUMI'])].loc[:,'trueUMI']
    real_true_UMIs = bkdf[~bkdf['trueUMI'].isin(bkdf['falseUMI'])].loc[:,'trueUMI']
    if len(non_true_UMIs) > 0:
        #skipping here the search for alternate pairings as implemented in zUMIs or UMIcountR for now -> can be added later
        #?umi_out[, diff := n.true-n.false]?
        bkdf = bkdf[bkdf["trueUMI"].isin(real_true_UMIs)]
    bkdf.drop(['hamdist','falseUMI_obs','trueUMI_obs'], axis = 1, inplace = True)
    return bkdf
def get_cellBC(read):
    cell = read.get_tag("SM")
    cell += read.get_tag("BC")
    return cell
def get_umi_mappings(inbam, chr, start, end, gene):
    #phase 1: load UMI strings from bam
    dat = {gene: {}}
    bam = pysam.AlignmentFile(inbam, 'rb')
    if chr not in bam.references:
        return None
    for read in bam.fetch(chr,start,end):
        if read.is_read2:
            continue
        if read.get_tag("XX") != "threep_BCUMI_read":
            continue
        cell = get_cellBC(read)
        umi = read.get_tag("UB")
        try:
            genetag = read.get_tag("GE")
        except:
            try:
                genetag = read.get_tag("GI")
            except:
                continue #read not tagged with a gene -> skip counting it
        if genetag != gene:
            continue #read tag inconsistent with gene to count -> skip
        if cell not in dat[gene]:
            dat[gene][cell] = []
        dat[gene][cell].append(umi)
    bam.close()
    #phase2: loop over cells, construct umi string mapping table
    out = {gene: {}}
    for cell in dat[gene].keys():
        if len(dat[gene][cell]) == 1:
            continue
        umi_mapping = directional_adjacency_corrections(dat[gene][cell])
        if umi_mapping is None:
            continue
        if cell not in out[gene].keys():
            out[gene][cell] = umi_mapping
        else:
            out[gene][cell] = pd.concat([out[gene][cell],umi_mapping])
    return out

def mol_pandas_to_dict(all_umi_mappings):
    #combine list to one dictionary
    all_umi_mappings_dict = {}
    for d in all_umi_mappings:
        if d != None:
            gene = list(d.keys())[0]
            if len(d[gene]) > 0:
                all_umi_mappings_dict[gene] = {}
                for cell in d[gene].keys():
                    all_umi_mappings_dict[gene][cell] = {}
                    d[gene][cell].reset_index(inplace = True)
                    for index, row in d[gene][cell].iterrows():
                        all_umi_mappings_dict[gene][cell][row['falseUMI']] = row['trueUMI']
    return(all_umi_mappings_dict)


def collect_bam_chunks(inpath, chrs, outpath):
    allpaths = [inpath+".tmp."+c+".bam" for c in chrs[:-1]]
    allpaths.append(inpath+".tmp."+"unmapped"+".bam")
    cat_args = ['-o', outpath]+allpaths
    pysam.cat(*cat_args)
    x = [os.remove(f) for f in allpaths]



def bam_update_umis(inbam, umi_mappings, threads):
    outbam = inbam.replace(".aligned_trimmed_genetagged_sorted.bam",".aligned_trimmed_genetagged_sorted_umicorrected.bam")
    global mols
    mols = umi_mappings
    chrs = pysam.idxstats(inbam).split('\n')
    chrs = [c.split('\t')[0] for c in chrs[:-1]]

    if threads > 8:
        pysam_workers = 2
        n_jobs = int(threads/2)
    else:
        pysam_workers = 1
        n_jobs = threads

    pool = mp.Pool(n_jobs)
    results = [pool.apply_async(correct_tags, (inbam,pysam_workers,chr, )) for chr in chrs]
    x = [r.get() for r in results]

    collect_bam_chunks(inpath = inbam, chrs = chrs, outpath = outbam)

    return outbam



def error_correct_umis(inbam, gtfpath, threads):
    gtf = parse_gtf_gene_name(gtfpath, None)
    starts = []
    ends = []
    seqnames = []
    genenames = []
    for gene_name, d in gtf.items():
        starts.append(d['start'])
        ends.append(d['end'])
        seqnames.append(d['seqid'])
        genenames.append(gene_name)

    pool = mp.Pool(threads)
    results = [pool.apply_async(get_umi_mappings, (inbam, seqnames[x], starts[x], ends[x], genenames[x], )) for x in range(len(genenames))]
    all_umi_mappings = [r.get() for r in results]
    pool.close()

    #combine list to one dictionary
    all_umi_mappings_dict = mol_pandas_to_dict(all_umi_mappings)

    print("Write corrected UMIs to bam file..")
    outbam = bam_update_umis(inbam, all_umi_mappings_dict, threads)
    return gtf, outbam
import sys
if __name__ == '__main__':
    error_correct_umis(sys.argv[1], sys.argv[2], int(sys.argv[3]))