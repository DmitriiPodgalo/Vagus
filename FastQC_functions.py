import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def read_file(file_path):
    """
    Function that processes input fastq file.
    This function reads fastq file in a temporary list of 4 strings
    (identifier, read sequence, plus, and quality score sequence).
    Then it adds this list into  the main 'parsed file' list.
    """
    parsed_file = []
    temp_list = []

    with open(file_path) as inf:

        for line in inf:
            temp_list.append(line.rstrip())

            if len(temp_list) == 4:
                parsed_file.append(temp_list)
                temp_list = []

    return parsed_file


# 1. Per base sequence quality

def max_read_length(lst):
    """
    Function that calculates the greatest read length, using a list of reads.
    """
    l_max = 0  # We set maximal length equal to zero

    # And then for each sequence in a list we compare its length to current l_max
    for k in lst:
        if len(k) > l_max:
            # Update l_max if len(k) is greater:
            l_max = len(k)

    return l_max


def calculate_quality_per_base(parsed_file):
    """
    Function that calculates quality score per base.
    The positions(bp) will be the keys in this dictionary
    and the values will contain the lists of quality scores in each position.
    We process the 3rd element of 'reads' list (quality score sequences)

    In the next step (plot drawing) we will need to turn dictionary into dataframe.
    Therefore, we need all values (which are lists) to contain equal number of elements.
    This number will be equal to the max read length.
    """
    qualities_per_base = dict()  # the new empty dictionary is created

    qualities = []
    for read in parsed_file:
        qualities.append(read[3])  # the 3rd element of 'read' list  is quality score sequences
    n = max_read_length(qualities)

    for qual_seq in qualities:
        # For each quality score sequence we add 1 value for each key of 'qualities_per_base' dictionary
        for i in range(n):

            # If this specific read contains >= i symbols, we add specific score to the list:
            if i < len(qual_seq):
                if i + 1 not in qualities_per_base:
                    # For the symbol in quality score line, we calculate its number in ascii
                    # table and subtract 33
                    qualities_per_base[i + 1] = [ord(qual_seq[i]) - 33]
                else:
                    qualities_per_base[i + 1].append(ord(qual_seq[i]) - 33)

            # If this specific read contains less then i symbols, we add None to the list:
            else:
                if i + 1 not in qualities_per_base:
                    qualities_per_base[i + 1] = [None]
                else:
                    qualities_per_base[i + 1].append(None)

    return qualities_per_base


def calculate_mean_quality_per_base(qualities_per_base):
    """
    Function that calculates average quality per base.
    """
    dict_mean_qual = dict()  # The new empty dictionary is created

    # We use 'quality_per_base' dictionary. It is a result of previous function,
    # containing the positions(bp) as the keys and the lists of quality scores in each position as values
    for key in qualities_per_base:
        # For each key we calculate average of value-list (we do not include None, while calculating mean):
        dict_mean_qual[key-1] = np.mean([i for i in qualities_per_base[key] if i is not None])

    return dict_mean_qual


def plot_per_base_seq_quality(parsed_file, DEFAULT_OUTPUT_DIR='./Report_data/'):
    """
    Function for 'Per base sequence quality' plot drawing.
    """
    # We turn 'qualities_per_base' dictionary into pandas data frame
    qualities_per_base = calculate_quality_per_base(parsed_file)
    dict_mean_qual = calculate_mean_quality_per_base(qualities_per_base)
    base = pd.DataFrame.from_dict(qualities_per_base)

    # Set default theme and remove margins:
    plt.margins(0)

    # Create colored background:
    plt.axhspan(0, 20, facecolor='#D14139', alpha=0.1)
    plt.axhspan(20, 28, facecolor='yellow', alpha=0.1)
    plt.axhspan(28, 42, facecolor='green', alpha=0.1)

    # And, finally, draw boxplot
    sns.boxplot(data=base, color="yellow", whis=[10, 90], showfliers=False,
                medianprops=dict(linewidth=0.5, color='#D14139'))

    # We process dictionary in order to set x and y for the plot:
    lists = sorted(dict_mean_qual.items())  # sorted by key, return a list of tuples
    x2, y2 = zip(*lists)  # unpack a list of pairs into two tuples

    # And also draw the mean plot:
    plt.plot(x2, y2, color="#1D2DD8", linewidth=0.5)

    # Some improvements to make the plot easier to read:
    plt.yticks(np.arange(0, 42, step=2))
    plt.xticks(np.arange(0, len(qualities_per_base) + 1, step=10))
    plt.title('Quality scores across all bases')
    plt.xlabel('Position in read')
    plt.gcf().set_size_inches(8, 6)
    plt.savefig(DEFAULT_OUTPUT_DIR + 'Per_base_sequence_quality.png', dpi=100, bbox_inches='tight')
    plt.close()

    # Checker:
    qualities_list = []
    for quality in qualities_per_base.values():
        quality_without_none = [j for j in quality if j is not None]
        qualities_list.append(quality_without_none)
    medians = []
    quartiles = []
    for lst in qualities_list:
        lst = np.array(lst)
        median = np.median(lst)
        medians.append(median)
        quartile = np.quantile(lst, 0.25)
        quartiles.append(quartile)
    medians = np.array(medians)
    quartiles = np.array(quartiles)
    if np.any(quartiles < 5) or np.any(medians < 20):
        return 'failure'
    elif np.any(quartiles < 10) or np.any(medians < 25):
        return 'warning'
    return 'good'


# 2. Per sequence quality scores

def mean_quality(qual):
    """
    Function that calculate mean quality score per sequence.
    """
    quality_score = 0
    #  For each symbol in quality score line, we calculate its number in ascii
    #  table and subtract 33
    for i in qual:
        quality_score += ord(i) - 33
    # Returns round value of mean quality score:
    return round(quality_score / len(qual))


def per_sequence_quality(parsed_file):
    """
    We process the 3rd element of 'reads' list - the list of quality scores sequences.
    Then we put it to function that calculate mean quality score per sequence.
    """
    qualities = []
    for read in parsed_file:
        qualities.append(read[3])  # the 3rd element of 'read' list  is quality score sequences

    # The new empty dictionary for average quality scores (keys)
    # and the number of sequences with that average (values):
    qual_and_numbers = dict()

    for qual_seq in qualities:
        # For each quality score sequence in the list we calculate the average quality score:
        n = mean_quality(qual_seq)
        if n in qual_and_numbers:
            # If the quality score is in the dictionary, we add 1 to its value
            qual_and_numbers[n] += 1
        else:
            # If it is not, we set 1 its value
            qual_and_numbers[n] = 1

    return qual_and_numbers


def plot_per_seq_quality_scores(parsed_file, DEFAULT_OUTPUT_DIR='./Report_data/'):
    """
    Function for 'Per sequence quality scores' plot drawing and checker.
    Checker returns 'failure' if the most frequently observed mean quality is below 27
    Returns 'warning' if the most frequently observed mean quality is below 20
    """
    d = per_sequence_quality(parsed_file)

    # We process dictionary in order to set x and y for the plot:
    lists = sorted(d.items())  # sorted by key, return a list of tuples
    x, y = zip(*lists)  # unpack a list of pairs into two tuples

    # And also draw the plot:
    plt.plot(x, y, color="#D14139", label='Average quality per read')
    # Some improvements to make the plot easier to read:
    plt.xlabel('Mean Sequence Quality (Pherd Score)')
    plt.title('Quality score distributions over all sequences')
    plt.xticks(np.arange(0, 40, step=2))
    plt.yticks(np.arange(0, max(d.values())+10000, step=10000))
    plt.legend(loc='upper right')
    plt.grid(alpha=0.5)

    plt.gcf().set_size_inches(8, 6)
    plt.savefig(DEFAULT_OUTPUT_DIR + 'Per_sequence_quality_scores.png', dpi=100, bbox_inches='tight')
    plt.close()

    # Checker:
    max_freq = max(list(d.values()))
    scores = np.array([score for score, freq in d.items() if freq == max_freq])
    if np.any(scores < 20):
        return 'failure'
    elif np.any(scores < 27):
        return 'warning'
    return 'good'


# 3. Per base sequence content

def per_base_nucleotides_proportion(parsed_file):
    """
    Function that creates and returns a list of 4 dictionaries (1 for each nucleotide)
    that contain number of base pair in the sequence as a key
    and the proportion of this nucleotide as a value.
    """
    seqs = []
    for read in parsed_file:
        seqs.append(read[1])  # We process the 1st element in the 'read' list - seq
    # We create 4 empty dictionaries:
    a_count, t_count, g_count, c_count = dict(), dict(), dict(), dict()

    # And count the number of each nucleotide for each base pair:
    for seq in seqs:
        for i in range(len(seq)):

            if seq[i] == 'A':
                if i+1 in a_count:
                    a_count[i + 1] += 1
                else:
                    a_count[i + 1] = 1
                    t_count[i + 1] = 0
                    g_count[i + 1] = 0
                    c_count[i + 1] = 0

            elif seq[i] == 'T':
                if i+1 in t_count:
                    t_count[i + 1] += 1
                else:
                    t_count[i + 1] = 1
                    a_count[i + 1] = 0
                    g_count[i + 1] = 0
                    c_count[i + 1] = 0

            elif seq[i] == 'G':
                if i+1 in g_count:
                    g_count[i + 1] += 1
                else:
                    g_count[i + 1] = 1
                    a_count[i + 1] = 0
                    t_count[i + 1] = 0
                    c_count[i + 1] = 0

            elif seq[i] == 'C':
                if i+1 in c_count:
                    c_count[i + 1] += 1
                else:
                    c_count[i + 1] = 1
                    a_count[i + 1] = 0
                    t_count[i + 1] = 0
                    g_count[i + 1] = 0

    # We create 4 new empty dictionaries for proportion:
    a_proportion, t_proportion, g_proportion, c_proportion = dict(), dict(), dict(), dict()

    # And calculate nucleotide proportion for each base:
    for count in a_count:
        general_count = a_count[count] + t_count[count] + g_count[count] + c_count[count]
        a_proportion[count] = a_count[count] / general_count * 100
        t_proportion[count] = t_count[count] / general_count * 100
        g_proportion[count] = g_count[count] / general_count * 100
        c_proportion[count] = c_count[count] / general_count * 100

    return [a_proportion, t_proportion, g_proportion, c_proportion]


def plot_per_base_seq_content(parsed_file, DEFAULT_OUTPUT_DIR='./Report_data/'):
    """
    Function for 'Per base sequence content' drawing and checker.
    The checker issues a 'warning' if the difference between A and T, or G and C is greater than 10% in any position.
    Returns 'failure' if the difference between A and T, or G and C is greater than 20% in any position.
    """
    lst_proportions = per_base_nucleotides_proportion(parsed_file)
    a_proportion = lst_proportions[0]
    t_proportion = lst_proportions[1]
    g_proportion = lst_proportions[2]
    c_proportion = lst_proportions[3]

    # The y's are different for each nucleotide, it is the proportion of nucleotide in the read (dict values):
    ay = list(a_proportion.values())
    ty = list(t_proportion.values())
    gy = list(g_proportion.values())
    cy = list(c_proportion.values())

    # The x is the same for each nucleotides, it is the number of base in the read (dict keys):
    x = list(a_proportion.keys())

    # And draw the plot for each nucleotide:
    plt.plot(x, ay, color="#D14139", linewidth=0.5, label='%A')
    plt.plot(x, ty, color="#1D2DD8", linewidth=0.5, label='%T')
    plt.plot(x, gy, color="green", linewidth=0.5, label='%G')
    plt.plot(x, cy, color="black", linewidth=0.5, label='%C')

    # Some improvements to make the plot easier to read:
    plt.xticks(np.arange(0, len(x)+2, step=10))
    plt.yticks(np.arange(0, 110, step=10))
    plt.xlabel('Position in read(bp)')
    plt.title('Sequence content across all bases')
    plt.legend(loc='upper right')
    plt.grid(alpha=0.5)
    plt.gcf().set_size_inches(8, 6)
    plt.savefig(DEFAULT_OUTPUT_DIR + 'Per_base_sequence_content.png', dpi=100, bbox_inches='tight')
    plt.close()

    # Checker:
    a, t, g, c = np.array(ay), np.array(ty), np.array(gy), np.array(cy)
    difference = np.append(np.absolute(a - t), np.absolute(g - c))
    if np.any(difference > 20):
        return 'failure'
    elif np.any(difference > 10):
        return 'warning'
    return 'good'


def main():
    # Sample input
    input_fastq = './Test_data/amp_res_2_passed.fastq'
    reads = read_file(input_fastq)
    DEFAULT_OUTPUT_DIR = './Report_data/'

    # 1
    plot_per_base_seq_quality(reads, DEFAULT_OUTPUT_DIR)

    # 2
    plot_per_seq_quality_scores(reads, DEFAULT_OUTPUT_DIR)

    # 3
    plot_per_base_seq_content(reads, DEFAULT_OUTPUT_DIR)


if __name__ == '__main__':
    main()
