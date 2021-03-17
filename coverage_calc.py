import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.signal as signal
import argparse
import os

kmer_dict = {}
def get_kmers(k, string):
    """Performs kmer extraction from reads and counts the number of repeats of that particular kmer

    Args:
        k (int): The desired length of the kmer
        string (str): Takes the read data input as a string
    
    Returns:
        dict: key:kmer, value:count
    """
    temp_reads = string.replace("\n", " ")
    end = len(temp_reads) - k + 1
    for start in range(0, end):
        kmer = temp_reads[start:start + k]
        if kmer.isalpha():
            if kmer in kmer_dict:
                kmer_dict[kmer] += 1
            else:
                kmer_dict[kmer] = 1
    return kmer_dict

cov_threshold = []
def plot_coverage(dict_kmer):
    """Plots kmer coverage plot and calculates the minimum threshold required to eliminate infrequent kmers, returns 
    true kmers. Basically, helps in removing errors from genomic sequencing data. 

    Args:
        dict_kmer (dict): key: kmer, value: count
    
    Returns:
        Histogram-density plot
        Threshold line marked by -- 
    """
    count_values = list(dict_kmer.values())
    count_log = np.log(count_values)

    mask = (count_log >= 0) & (count_log <= 4)

    plt.subplots(figsize=(12,8))
    sns.distplot(count_log, bins='auto')

    minindices = signal.argrelextrema(count_log[mask], np.less)
    local_minima = count_log[mask][minindices]
    threshold = (local_minima > 0) & (local_minima < 1)
    cov_threshold.append(sorted(local_minima[threshold])[0])

    plt.axvline(x=cov_threshold,color='red',linestyle='--')
    plt.xlabel("Coverage")
    plt.ylabel("Count")
    plt.autoscale()
    plt.savefig(plot_dir, dpi=1000)

kmer_temp = {}
def true_kmers(kmer_dict):
    """Extracts the number of true kmers by eliminating the kmers below the calculated threshold value.

    Args:
        dict_kmer (dict): key: kmer, value: counting
    
    Returns:
        int: Total number of kmers after removal of infrequenctly occuring kmers 
    """
    for k,v in kmer_dict.items():
        t = np.log(v)
        if t > cov_threshold:
            kmer_temp[v] = v
    return "The number of kmers generated after slicing off the less frequent kmers is {}".format(len(kmer_temp))


def main():
    print('Starting coverage plot calculation, will take upto 118 to 130s')

    data = open(read_file)
    reads = data.read()

    kmers = get_kmers(kmer_length, reads)
    plot_coverage(kmers)
    
    print("The number of kmers generated using the sliding window is {}".format(len(kmer_dict)))
    
    true_kmer_num = true_kmers(kmers)
    print(true_kmer_num)
    
if __name__ == "__main__":
    kmer_list = [31,41,51,61,71,81,91,101]

    parser = argparse.ArgumentParser()
    parser._action_groups.pop()
    required = parser.add_argument_group('Required Argument')
    optional = parser.add_argument_group('Optional Argument')

    required.add_argument('-reads', help='Submit read file', required=True)

    required.add_argument('-kmer_len', help='Insert the desired length of kmer', 
                            type=int, required=True, choices=kmer_list)

    optional.add_argument('-plot_name', help='Desired name for the coverage plot', 
                            required=False, default='coverage_plot.png')
    
    args = parser.parse_args()
    
    read_file = args.reads
    kmer_length = args.kmer_len
    plot_name = args.plot_name

    home_dir = os.environ['HOME']
    output_dir = os.path.join(home_dir, 'coverage plot')
    plot_dir = os.path.join(output_dir, plot_name)

    os.makedirs(output_dir, exist_ok=True)

    main()