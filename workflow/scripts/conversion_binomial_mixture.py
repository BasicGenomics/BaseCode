import pysam
from collections import Counter
import pandas as pd
import sys
import os
import numpy as np
import argparse
import time
from collections import Counter
from joblib import delayed,Parallel
from multiprocessing import Process, Manager
from pyfaidx import Fasta
import itertools
from tools import parse_gtf, filterGeneDict, determine_gene_tag

import pandas as pd
import numpy as np
import os

import numpy as np
import torch
from torch.distributions.uniform import Uniform
from torch.distributions.binomial import Binomial

if torch.cuda.is_available():
    import torch.cuda as t
    device = torch.device('cuda:0')
else:
    import torch as t
    device = torch.device('cpu')


# if torch.cuda.is_available():
#     print("cuda is available")
#     import torch.cuda as t
# else:
#     print("cuda is unavailable")
#     import torch as t

class BinomialMixture():
    '''
    Binomial Mixture Model.
    The BinomialMixture object implements the expectation-maximization (EM) algorithm for
    fitting mixture-of-binomial-distributions models.
    The BinomialMixture.fit() method fits the given data, a number of (n,N) pairs, by a
    mixture of n Binomial Distributions, where n needs to be specified by the parameter
    `n_components` by the user.
    The BinomialMixture.predict() method, taking in (n, N) pairs, predicts the probabilities
    for each (n,N) pair to belong to any of the fitted Binomial-Distribution components.
    '''
    def __init__(self, n_components, tolerance=1e-5, max_step=1000, verbose=True):
        '''
        Initialize the object.
        Input
        -----
        n_components: int
                      number of Binomial Distributions in the Mixture.
        tolerance: Float
                   the object fits the data by EM iteration. tolerance is used to define one
                   of the conditions for the iteration to stop. This condition says that
                   the iteration continues until the change of the parameters within the
                   current iteration is greater than tolerance.
        max_step: int
                  the maximum number of iteration steps.
        verbose: Boolean
                 whether to print out the information of the fitting process.
        '''
        self.K = n_components # int, number of Binomial distributions in the Mixture
        self.tolerance = tolerance
        self.max_step = max_step
        self.verbose = verbose

        # initialize the pi_list
        pi_list = Uniform(low=1e-6, high=1e0-1e-6).sample([self.K-1]).to(device)
        pi_K = t.FloatTensor([1e0]) - pi_list.sum()
        self.pi_list = torch.cat([pi_list, pi_K], dim=0)

        # initialize the theta_list
        self.theta_list = Uniform(low=1e-6, high=1e0-1e-6).sample([self.K])

    def calc_logL(self, N_ls, n_ls):
        '''
        Calculate the log likelihood.
        Input
        -----
        N_ls:  a [S] shape tensor = [N1, N2, ..., NS]
        n_ls:  a [S] shape tensor = [n1, n2, ..., nS]
        Output
        ------
        log_likelihood: the conditional probability of the parameters
                        (pi and theta) given the observed data (N, n)
        '''
        S = len(N_ls)
        pi_list = self.pi_list
        theta_list = self.theta_list
        #K = self.K

        # log_binom_mat has shape (S,K), element_{i,l} = log_Binomial(ni|Ni, theta_l)
        # log with natural base.
        log_binom_mat = Binomial(N_ls.reshape(S,1), theta_list.reshape(1,self.K)).log_prob(n_ls.reshape(S,1))

        # mean_log_binom, the mean value of all elements in log_binom_mat.
        c = torch.mean(log_binom_mat)

        # binom_mat has shape (S,K), element_{i,l} = Binomial(ni|Ni, theta_l)
        binom_mat = torch.exp(log_binom_mat - c)

        # log_likelihood = sum_{i=1}^{S} log(prob_i), this is a real number
        log_likelihood = S*c + torch.sum(torch.log(torch.matmul(binom_mat, pi_list)))

        return log_likelihood

    def calc_Posterior(self, N_ls, n_ls):
        '''
        Calculate the posterior probabilities.
        Input
        -----
        N_ls:  a [S] shape tensor = [N1, N2, ..., NS]
        n_ls:  a [S] shape tensor = [n1, n2, ..., nS]
        Output
        ------
        Posterior: a tensor with shape (K,S)
                   its element_{m,i} = P(zi=m|ni,Theta_old) which is
                   the posterior probability of the i-th sample belonging
                   to the m-th Binomial distribution.
        '''
        pi_list = self.pi_list
        theta_list = self.theta_list
        # K = self.K
        S = len(N_ls)

        # shape = (K,K) with theta_ratio_{m,l} = theta_l/theta_m, m-th row, l-th column
        theta_ratio = torch.div(theta_list.reshape(1,self.K), theta_list.reshape(self.K,1))

        # shape = (K,K), element_{ml} = (1-theta_l)/(1-theta_m)
        unity_minus_theta_ratio = torch.div((1e0 - theta_list).reshape(1,self.K), (1e0 - theta_list).reshape(self.K,1))

        # shape = (K,K), element_{m,l} = (theta_l/theta_m) * [(1-theta_l)/(1-theta_m)]
        mixed_ratio = torch.mul(theta_ratio, unity_minus_theta_ratio)

        # shape = (K,K,S) with element_{m,l,i} = [(theta_l/theta_m)*(1-theta_l)/(1-theta_m)]^ni
        # its element won't be either 0 or infty no matther whether theta_l > or < theta_m
        mixed_ratio_pow = torch.pow(theta_ratio.reshape(self.K, self.K, 1), n_ls)
        mixed_ratio_pow = torch.clamp(mixed_ratio_pow, min=0.0, max=1e15)

        # shape = (K,K,S) with element_{m,l,i} = [ (1-theta_l)/(1-theta_m) ]^(Ni-2ni)
        # its element may be infty if theta_l<<theta_m, or 0 if theta_l >> theta_m
        unity_minus_theta_ratio_pow = torch.pow(unity_minus_theta_ratio.reshape(self.K,self.K,1), N_ls-2.0*n_ls)
        unity_minus_theta_ratio_pow = torch.clamp(unity_minus_theta_ratio_pow, min=0.0, max=1e15)

        # In below, we multiply the element of mixed_ratio_pow and the element of unity_minus_theta_ratio_pow,
        # and there won't be nan caused by 0*infty or infty*0 because the element in mixed_ratio_pow won't be 0 or infty.
        # Thus we make sure there won't be nan in Posterior.

        # element-wise multiply, pow_tensor has shape(K,K,S), element_{m,l,i} = (theta_l/theta_m)^ni * [(1-theta_l)/(1-theta_m)]^(Ni-ni).
        # Note that torch.mul(a, b) would broadcast if a and b are different in shape & they are
        # broadcastable. If a and b are the same in shape, torch.mul(a,b) would operate element-wise multiplication.

        pow_tensor = torch.mul(mixed_ratio_pow, unity_minus_theta_ratio_pow)

        # pi_ratio has shape (K,K) with element_{m,l} = pi_l/pi_m
        pi_ratio = torch.div(pi_list.reshape(1,self.K), pi_list.reshape(self.K,1))

        # posterior probability tensor, Pzim = P(zi=m|ni,Theta_old)
        # shape (K,S), element_{m,i} = P(zi=m|ni,Theta_old)
        Posterior = torch.pow(torch.matmul(pi_ratio.reshape(self.K,1,self.K), pow_tensor), -1e0).reshape(self.K,S)

        return Posterior

    def calc_params(self, N_ls, n_ls, Posterior):
        '''
        Calculate the parameters, pi_list and theta_list, using EM iteration.
        Input
        -----
         N_ls: tensor of shape [S]
               input data
         n_ls: tensor of shape [S]
               input data
         Posterior: tensor of shape (K,S)
                    contains the posterior probabilites for the ith sample belonging
                    to the m-th Binomial component. Here i=1,2..,K and m=1,2,..,K
        '''
        # update pi_list
        # torch.sum(tensor, n) sum over the n-th dimension of the tensor
        # e.g. if tensor'shape is (K,S) and n=1, the resulting tensor has shape (K,)
        # the m-th element is the sum_{i=1}^{S} tensor_{m,i}
        pi_list = torch.sum(Posterior,1)/len(N_ls)

        # update theta_list
        theta_list = torch.div(torch.matmul(Posterior, n_ls), torch.matmul(Posterior, N_ls))

        return pi_list, theta_list

    def initialize_pi_theta(self, N_ls, n_ls):
        '''
        Initialize pi_list and theta_list based on the input data
        N_ls and n_ls.
        This is a better initialization than completely random
        initilaization, since it takes into account information from N_ls and
        n_ls, and thus a better initialization should help the EM iterations
        in the .fit() method to converge faster.
        Input
        -----
        N_ls: 1D tensor
        n_ls: 1D tensor
        Output
        ------
        None
        '''

        # Conver to numpy.array
        N_array = np.array(N_ls.cpu())
        n_array = np.array(n_ls.cpu())

        # list of n/N sorted in ascending order
        ratio_ls = np.sort(n_array/N_array)
        ratio_ls += 1e-5
        # reproducibility
        np.random.seed(seed=123)

        # pick K random integer indice, [index_1, ..., index_K]
        random_indice = np.sort(np.random.choice(len(ratio_ls),self.K))
        # theta are the ratio at the random indice, { ratio_ls[index_k] | k=1,2,...,K }
        if self.K == 2:
            theta_array = np.array([0.01, 0.05])
        elif self.K == 1:
            theta_array = np.array([0.01])
        #theta_array = ratio_ls[random_indice]

        # the proportion of the midpoint of each pair of consecutive indice
        # (index_k + index_{k+1})/2/len(ratio_ls), for k=1,2,...,K-1
        acc_portion_ls = (random_indice[1:]+random_indice[:-1])/2.0/len(ratio_ls)
        acc_portion_ls = np.append(acc_portion_ls, 1.0)
        # initialize pi_list using the portions of indice
        pi_array = np.insert(acc_portion_ls[1:] - acc_portion_ls[:-1], obj=0, values=acc_portion_ls[0])
        

        # convert numpy arrays to torch tensors.
        self.theta_list = t.FloatTensor(theta_array)
        self.pi_list = t.FloatTensor(pi_array)

    def fit(self, N_ls, n_ls):
        '''
        Fit the input data, N_ls and n_ls, by a mixture of K Binomial Distributions.
        It uses the Expectation-Maximization algorithm to iteratively fit the data by
        a Binomial-Mixture described by parameters, pi_list and theta_list. The distribution's
        probability mass function is
        BM(n|N, pi_list, theta_list) = sum_{m=1}^{K} pi_m * Binomial(n|N, theta_m),
        where pi_m is the m-th element in pi_list, theta_m is the m-th element in
        theta_list, and Binomial(n|N,theta) = N!/(n!*(N-n)!) * theta^n *(1-theta)^(N-n).
        Input
        -----
         N_ls: tensor of shape [S]
               input data
         n_ls: tensor of shape [S]
               input data
        Output
        ------
        None
        Note: this method updates the self.pi_list and self.theta_list using the EM
        Algorithm iterations until they converge within the specified self.tolerance
        or until the iteration steps exceed the self.max_step.
        '''

        # Initialize pi and theta
        self.initialize_pi_theta(N_ls, n_ls)

        # calculate the initial log_likelihood
        log_likelihood = self.calc_logL(N_ls, n_ls)

        # initialize the change of log_likelihood named `delta_log_likelihood` and the iteration step called `iter_step`
        delta_log_likelihood = torch.norm(log_likelihood)
        iter_step = 0

        # old theta and pi
        pi_list_old = self.pi_list
        theta_list_old = self.theta_list
        params_old = torch.stack([pi_list_old, theta_list_old], dim=1)
        delta_params = torch.norm(params_old.reshape(2*self.K,))

        # The iteration stops when either of the following two conditions is broken first
        cond_step = t.BoolTensor([iter_step < self.max_step])
        # the condition below is that "the params are still changing much"
        cond_params = (delta_params >self.tolerance)

        import time
        start_time = time.time()

        while cond_step & cond_params:

            # posterior probability tensor, Pzim = P(zi=m|ni,Theta_old)
            # shape (K,S), element_{m,i} = P(zi=m|ni,Theta_old)
            Posterior = self.calc_Posterior(N_ls, n_ls)

            # calculate the new pi_list and theta_list
            self.pi_list, self.theta_list =  self.calc_params(N_ls, n_ls, Posterior)

            # update params and delta_params
            params = torch.stack([self.pi_list, self.theta_list], dim=1)
            delta_params = torch.norm((params-params_old).reshape(2*self.K,))

            # calculate the new log_likelihood
            log_likelihood_new = self.calc_logL(N_ls, n_ls)

            # calculate the change of the log-likelihood
            delta_log_likelihood = torch.norm(log_likelihood_new - log_likelihood)

            # update params
            params_old = params

            # update log_likelihood
            log_likelihood = log_likelihood_new

            # increase iter_step by 1
            iter_step += 1

            # update the conditions for the while loop
            # cond_params = (delta_params > epsilon)
            cond_params = (delta_params > self.tolerance)
            cond_step = t.BoolTensor([iter_step < self.max_step])

            if self.verbose & (iter_step % 10 == 0):
                print(f"Iteration {iter_step}:")
                print(f"delta_log_likelihood = {delta_log_likelihood}")
                print(f"log_likelihood ={log_likelihood}")
                print(f"{params}")

        if self.verbose:
            print(f"used {time.time()-start_time}")

    def predict(self, N_ls, n_ls):
        '''
        Calculate the posterior probabilities.
        Input
        -----
        N_ls:  1D tensor of any arbitrary length = len.
               It could be new data that BinomialMixture object has
               never seen before or the data that the BinomialMixture
               object has fitted on.
        n_ls:  1D tensor of the same length as N_ls.
        Output
        ------
        Posterior: a tensor with shape [K,len]
                   its element_{m,i} = P(zi=m|ni,Theta_old) which is
                   the posterior probability of the i-th sample belonging
                   to the m-th Binomial distribution.
        '''
        Posterior = self.calc_Posterior(N_ls, n_ls)
        return Posterior

    def calc_AIC_BIC(self, N_ls, n_ls):
        '''
        Calculate the Akaike Information Criterion (AIC) and
        the Bayesian Information Criterion (BIC)
        Input
        -----
        N_ls: 1D tensor of shape [S]
        n_ls: 1D tensor of shape [S]
        Output
        ------
        AIC: torch.FloatTensor of shape [1]
        BIC: torch.FloatTensor of shape [1]
        '''
        # the number of samples in N_ls
        S = len(N_ls)

        # calculate the initial log_likelihood
        log_likelihood = self.calc_logL(N_ls, n_ls)

        # calculate Akaike Information Criterion (AIC)
        AIC = -2.0/float(S)*log_likelihood + 2.0*(2.0*float(self.K)+1.0)/float(S)
        # Bayesian Information Criterion
        BIC = -2.0*log_likelihood + np.log(float(S))*(2.0*float(self.K)+1.0)

        return AIC, BIC

    def p_value(self, N_ls_new, n_ls_new, side='right'):
        '''
        Calculate the p-value of new (N,n) pairs.
        After the BinomialMixture object has fitted on N_ls and n_ls,
        we can use the fitted BinomialMixture distribution to calculate
        the p-value for new (N, n) pairs. The p-value indicates how surprisingly
        rare the new (N,n) pair is if this pair follows the Binomial-Mixture-
        Distribution that is fitted on N_ls and n_ls. The probability of n given
        N and the Binomial-Mixture-model is
        BM(n|N, pi_list, theta_list) = sum_{m=1}^{K} pi_m * Binomial(n|N, theta_m),
        where pi_m is the m-th element in pi_list, theta_m is the m-th element in
        theta_list, and Binomial(n|N,theta) = N!/(n!*(N-n)!) * theta^n *(1-theta)^(N-n).
        BM.cmf(n0) is the cumulative mass function, i.e. the probability for n <= n0.
        BM.sf(n0) = 1 - BM.cmf(n0) = the probability for n > n0.
        BM.pmf(n0) = the probability for n = n0.
        Input
        -----
        N_ls_new: 1D tensor
                  new data that the BinomialMixture object has never seen before.
        n_ls_new: 1D tensor
                  new data
        side: str
              'right', 'left', 'both'
              If side = 'right', calculate the right-side p-value, namely
              p-value = probability for n >= n0 = BM.sf(n0) + BM.pmf(n0)
              If side = 'left', calculate the left-side p-value, namely
              p-value = BM.cmf(n0)
              If side = 'both', calculate the double-sided p-value, namely
              p-value = 2 * min[BM.cmf(n0), BM.sf(n0)+BM.pmf(n0)]
        Output
        ------
        p_value_ls: 1D tensor of the same length as N_ls and n_ls.
                    contains the p-value for each (N,n) pair.
        '''

        from scipy.stats import binom

        # save to host.cpu and then convert to numpy
        # we need to do it because scipy needs it.
        N_ls_new = N_ls_new.cpu().numpy()
        n_ls_new = n_ls_new.cpu().numpy()
        theta_list = self.theta_list.cpu()
        pi_list = self.pi_list

        p_value_ls = []
        for i in range(len(N_ls_new)):
            N = N_ls_new[i]
            n = n_ls_new[i]

            if side == 'right':
                sf_list = binom.sf(n, N, theta_list)
                pmf_list = binom.pmf(n, N, theta_list)
                p_value = torch.matmul(pi_list, t.FloatTensor(sf_list + pmf_list))
                p_value_ls.append(p_value)

            elif side == 'left':
                cmf_list = binom.cdf(n, N, theta_list)
                p_value = torch.matmul(pi_list, t.FloatTensor(cmf_list))
                p_value_ls.append(p_value)

            elif side == 'both':

                sf_list = binom.sf(n, N, theta_list)
                pmf_list = binom.pmf(n, N, theta_list)
                p_value_right = torch.matmul(pi_list, t.FloatTensor(sf_list + pmf_list))

                cmf_list = binom.cdf(n, N, theta_list)
                p_value_left = torch.matmul(pi_list, t.FloatTensor(cmf_list))

                p_value = 2.0*torch.min(p_value_right, p_value_left)
                p_value_ls.append(p_value)

        return t.FloatTensor(p_value_ls)

def base_conversion_complement(base_conversions):
    complement_dict_lower = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    complement_dict_upper = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement = [''.join([complement_dict_lower[bc[0]], complement_dict_upper[bc[1]]]) for bc in base_conversions]
    return complement

def get_whole_sequence(fasta_file, g_dict):

    fasta_ref = Fasta(fasta_file)
    whole_sequence = fasta_ref[g_dict['seqid']][g_dict['start']:g_dict['end']].seq.upper()

    return whole_sequence





def get_insertions_locs(cigtuples):
    insertion_locs = []
    l = 0
    for c in cigtuples:
        if c[0] == 0:
            l += c[1]
        elif c[0] == 1:
            for i in range(c[1]):
                insertion_locs.append(l)
                l += 1
    return insertion_locs

def compare_to_reference(read_seq, read_qual, ref_seq):
    total_content = Counter({'A': 0, 'T': 0, 'C': 0, 'G': 0})
    specific_content = Counter()
    for i, (read_base, base_qual, ref_base) in enumerate(zip(read_seq, read_qual, ref_seq)):
        if base_qual < 15:
            continue
        if read_base == ref_base:
            total_content.update({ref_base: 1})
        else:
            k = ''.join((ref_base.lower(),read_base))
            specific_content.update({k: 1})
    return specific_content, total_content

def find_base_conversions(read, whole_sequence, start_offset):

    intervals = read.get_blocks()
    offset = 0
    read_seq = read.query_alignment_sequence
    read_qual = read.query_alignment_qualities

    cigtuples = read.cigartuples
    insertion_locs = get_insertions_locs(cigtuples)

    if len(insertion_locs) > 0:
        read_seq = "".join([char for idx, char in enumerate(read_seq) if idx not in insertion_locs])
        read_qual = [qual for idx, qual in enumerate(read_qual) if idx not in insertion_locs]
    total_content = Counter()
    specific_content = Counter()
    for interval in intervals:
        interval_start = interval[0]
        interval_end = interval[1]-1
        interval_length = interval_end-interval_start
        start_idx = interval_start-start_offset
        end_idx = interval_end-start_offset
        if start_idx < 0:
            interval_start -= start_idx
            offset -= start_idx
            interval_length += start_idx
            start_idx = 0
        if end_idx > len(whole_sequence):
            end_idx = len(whole_sequence)
            interval_length = end_idx - offset
        ref_seq = whole_sequence[start_idx:end_idx]
        read_seq_interval = read_seq[offset:offset + interval_length]
        read_qual_interval = read_qual[offset:offset + interval_length]
        offset += interval_length + 1
        specific_content_interval, total_content_interval = compare_to_reference(read_seq_interval, read_qual_interval, ref_seq)
        specific_content.update(specific_content_interval)
        total_content.update(total_content_interval)
    return specific_content, total_content



def get_consistency(read, g_dict, read_type):
    if read.is_reverse and (g_dict['strand'] == '-'):
        if read_type == 'FP_read':
            return True
        else:
            return False
    elif not read.is_reverse and (g_dict['strand'] == '+'):
        if read_type == 'FP_read':
            return True
        else:
            return False
    else:
        if read_type == 'FP_read':
            return False
        else:
            return True

def func(bamfile,fasta_file, g_dict, conversion_plus, conversion_minus, coverage_plus, coverage_minus):
    whole_sequence = get_whole_sequence(fasta_file, g_dict)
    bam = pysam.AlignmentFile(bamfile)
    res_dict = {}
    for read in bam.fetch(g_dict['seqid'], g_dict['start'], g_dict['end']):
        gene_name = determine_gene_tag(read)
        if gene_name != g_dict['gene_id'] or read.mapping_quality < 1:
            continue
        sample = read.get_tag('SM')
        read_type = read.get_tag('XX')
        if sample not in res_dict:
            res_dict[sample] = {}
        if read_type not in res_dict[sample]:
            res_dict[sample][read_type] = {}
            res_dict[sample][read_type][True] = {}
            res_dict[sample][read_type][True][True] = Counter()
            res_dict[sample][read_type][True][False] = Counter()
            res_dict[sample][read_type][False] = {}
            res_dict[sample][read_type][False][True] = Counter()
            res_dict[sample][read_type][False][False] = Counter()
        consistent_read = get_consistency(read, g_dict, read_type)
        specific_conversions, total_content = find_base_conversions(read, whole_sequence, g_dict['start'])
        if g_dict['strand'] == '+':
            if conversion_plus in specific_conversions:
                res_dict[sample][read_type][consistent_read][read.is_read1].update({(specific_conversions[conversion_plus], total_content[coverage_plus]): 1})
            else:
                if coverage_plus in total_content:
                    res_dict[sample][read_type][consistent_read][read.is_read1].update({(0, total_content[coverage_plus]): 1})
                else: 
                    res_dict[sample][read_type][consistent_read][read.is_read1].update({(0, 0): 1})
        else:
            if conversion_minus in specific_conversions:
                res_dict[sample][read_type][consistent_read][read.is_read1].update({(specific_conversions[conversion_minus], total_content[coverage_minus]): 1})
            else:
                if coverage_minus in total_content:
                    res_dict[sample][read_type][consistent_read][read.is_read1].update({(0, total_content[coverage_minus]): 1})
                else: 
                    res_dict[sample][read_type][consistent_read][read.is_read1].update({(0, 0): 1})
    return res_dict
def make_arrays(counter):
    sc_list = []
    tc_list = []
    for (sc, tc), count in counter.items():
        sc_list.extend([sc for i in range(count)])
        tc_list.extend([tc for i in range(count)])
    sc_array = np.array(sc_list)
    tc_array = np.array(tc_list)
    return sc_array, tc_array
def make_pc_and_pe_estimation(counter, sample, read_type, consistency, read_pair):
    sc_array, tc_array = make_arrays(counter)
    tc_array_filtered = tc_array[tc_array>0]
    sc_array_filtered = sc_array[tc_array>0]
    filtering = sc_array_filtered <= tc_array_filtered
    tc_array_filtered = tc_array_filtered[filtering]
    sc_array_filtered = sc_array_filtered[filtering]
    tc_tensor = torch.Tensor(tc_array_filtered)
    sc_tensor = torch.Tensor(sc_array_filtered)

    BM_1 = BinomialMixture(n_components=1)
    BM_2 = BinomialMixture(n_components=2)
    BM_1.fit(tc_tensor, sc_tensor)
    BM_2.fit(tc_tensor, sc_tensor)

    BM_2_pi = pd.Series([t.numpy() for t in BM_2.pi_list])

    BM_2_pi.index = ['BM_2_pi_0', 'BM_2_pi_1']

    BM_2_theta = pd.Series([t.numpy() for t in BM_2.theta_list])
    BM_2_theta.index = ['BM_2_theta_0', 'BM_2_theta_1']



    BM_1_pi = pd.Series(BM_1.pi_list)
    BM_1_pi.index = ['BM_1_pi']

    BM_1_theta = pd.Series(BM_1.theta_list)
    BM_1_theta.index = ['BM_1_theta']

    BM_1_selection = pd.Series([t.numpy() for t in BM_1.calc_AIC_BIC(tc_tensor, sc_tensor)])
    BM_1_selection.index = ['AIC_BM1', 'BIC_BM1']

    BM_2_selection = pd.Series([t.numpy() for t in BM_2.calc_AIC_BIC(tc_tensor, sc_tensor)])
    BM_2_selection.index = ['AIC_BM2', 'BIC_BM2']
    
    info_series = pd.Series([sample, read_type, consistency, read_pair])
    info_series.index = ['sample', 'read_type', 'consistency', 'read_pair']

    full_series = pd.concat([BM_2_pi, BM_2_theta, BM_2_selection, BM_1_pi, BM_1_theta, BM_1_selection, info_series])
    full_series.name = sample
    return full_series

def get_n_reads(counter, sample, read_type, consistency, read_pair):
    sc_array, tc_array = make_arrays(counter)
    tc_array_filtered = tc_array[tc_array>0]
    sc_array_filtered = sc_array[tc_array>0]
    filtering = sc_array_filtered <= tc_array_filtered
    tc_array_filtered = tc_array_filtered[filtering]
    sc_array_filtered = sc_array_filtered[filtering]
    info_series = pd.Series([sample, read_type, consistency, read_pair, len(tc_array_filtered)])
    info_series.index = ['sample', 'read_type', 'consistency', 'read_pair', 'n_reads']
    info_series.name = sample
    return info_series

def count_reads_with_conversions(counter, sample, read_type, consistency, read_pair, tc_cutoff, fraction_cutoff):
    sc_array, tc_array = make_arrays(counter)
    tc_array_filtered = tc_array[tc_array>tc_cutoff]
    sc_array_filtered = sc_array[tc_array>tc_cutoff]
    filtering = sc_array_filtered <= tc_array_filtered
    tc_array_filtered = tc_array_filtered[filtering]
    sc_array_filtered = sc_array_filtered[filtering]
    fraction_converted = sc_array_filtered/tc_array_filtered
    more_than_fraction_cutoff = fraction_converted > fraction_cutoff
    mean_conversion_rate_above = np.mean(fraction_converted[more_than_fraction_cutoff])
    mean_conversion_rate_below = np.mean(fraction_converted[~more_than_fraction_cutoff])
    info_series = pd.Series([sample, read_type, consistency, read_pair, len(tc_array_filtered), np.sum(more_than_fraction_cutoff), mean_conversion_rate_above, mean_conversion_rate_below])
    info_series.index = ['sample', 'read_type', 'consistency', 'read_pair', 'n_reads', 'n_converted', 'rate_above', 'rate_below']
    info_series.name = sample
    #plt.hist(fraction_converted[more_than_fraction_cutoff], label='above', bins=30)
    #plt.hist(fraction_converted[~more_than_fraction_cutoff], label='below', bins=30)
    #plt.legend()
    #plt.show()
    return info_series

def main():
    parser = argparse.ArgumentParser(description='Calculate conversion rates', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--input',metavar='input', type=str, help='Input .bam file')
    parser.add_argument('-f','--fasta', metavar='fasta', type=str, help='Reference fasta file')
    parser.add_argument('-g','--gtf', metavar='gtf', type=str, help='gtf file with gene information')
    parser.add_argument('-p','--prefix', help='Run prefix')
    parser.add_argument('-o','--out_dir', help='Output directory')
    parser.add_argument('-t', '--threads', metavar='threads', type=int, default=1, help='Number of threads')
    parser.add_argument('-bc', '--base-conversions', nargs='+', type=str, default=['aG'], help='Base conversions used to reconstruct')

    args = parser.parse_args()

    bam_infile = args.input
    fastafile = args.fasta
    gtffile = args.gtf
    prefix = args.prefix
    threads = args.threads
    out_dir = args.out_dir
    isExist = os.path.exists(out_dir)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(out_dir)

    
    
    conversion_plus = 'gA'
    conversion_minus = 'cT'
    coverage_plus = 'G'
    coverage_minus = 'C'

    gene_dict = filterGeneDict(parse_gtf(gtffile, None), bam_infile)

    res = Parallel(n_jobs=threads, verbose = 3, backend='multiprocessing')(delayed(func)(bam_infile,fastafile, g_dict, conversion_plus, conversion_minus, coverage_plus, coverage_minus) for g,g_dict in gene_dict.items())

    full_res_dict = {}
    for gene_res in res:
        for sample, sample_dict in gene_res.items():
            if sample not in full_res_dict:
                full_res_dict[sample] = {}
            for read_type, read_type_dict in sample_dict.items():
                if read_type not in full_res_dict[sample]:
                    full_res_dict[sample][read_type] = {}
                for consistency, consistency_dict in read_type_dict.items():
                    if consistency not in full_res_dict[sample][read_type]:
                        full_res_dict[sample][read_type][consistency] = {}
                    for read_pair, read_pair_counter in consistency_dict.items():
                        if read_pair not in full_res_dict[sample][read_type][consistency]:
                            full_res_dict[sample][read_type][consistency][read_pair] = Counter()
                        full_res_dict[sample][read_type][consistency][read_pair].update(read_pair_counter)
    estimation_list = []

    for sample, sample_dict in full_res_dict.items():
        for read_type, read_type_dict in sample_dict.items():
            for consistency, consistency_dict in read_type_dict.items():
                for read_pair, read_pair_counter in consistency_dict.items():
                    try:
                        estimation_list.append(make_pc_and_pe_estimation(read_pair_counter, sample, read_type, consistency, read_pair))
                    except:
                        continue
    df_estimation = pd.DataFrame(estimation_list)
    df_estimation['consistency'] = df_estimation.apply(lambda row: not row['consistency'] if not row['read_pair'] else row['consistency'], axis=1)

    df_estimation_consistent = df_estimation[df_estimation['consistency']]

    df_estimation_consistent['BM_2_theta_1'] = df_estimation_consistent['BM_2_theta_1'].astype(float)

    df_estimation_consistent['BM_2_theta_0'] = df_estimation_consistent['BM_2_theta_0'].astype(float)

    df_estimation_consistent['BM_2_pi_1'] = df_estimation_consistent['BM_2_pi_1'].astype(float)

    df_estimation_consistent['BM_2_pi_0'] = df_estimation_consistent['BM_2_pi_0'].astype(float)

    df_estimation_consistent.to_csv('{}/{}_conversion_binomial_mixture.csv'.format(out_dir, prefix))

if __name__ == "__main__":
    main()
