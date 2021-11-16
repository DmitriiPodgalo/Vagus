# Function that processes input fastq file
def process_file(input_fastq):
    with open(input_fastq, 'r') as f:
        reads = f.readlines()
        # This function reads fastq file in a list and make
        # 4 different new lists:
        identifier = [sr.strip() for sr in reads[0::4]]
        seq = [sr.strip() for sr in reads[1::4]]
        plus = [sr.strip() for sr in reads[2::4]]
        quality = [sr.strip() for sr in reads[3::4]]
        # Then it makes list of identifiers, read sequences, pluses, and quality score sequences
        reads = [identifier, seq, plus, quality]
    return reads

# 1. Per base sequence quality


# Function that calculates the greatest read length, using a list of reads
def max_read_length(lst):
    # We set maximal length equal to zero
    l_max = 0
    # And then for each sequence in a list we compare its length to current l_max
    for k in lst:
        if len(k) > l_max:
            # Update l_max if len(k) is greater:
            l_max = len(k)
    return l_max


# Function that calculates quality score per base
def calculate_quality_per_base(reads):
    # The new empty dictionary is created:
    qualities_per_base = dict()
    # The positions(bp) will be the keys in this dictionary
    # And the values will contain the lists of quality scores in each position
    # We process the 3rd element of 'reads' list (quality score sequences)
    qualities = reads[3]
    # In the next step (plot drawing) we will need to turn dictionary into dataframe
    # Therefore, we need all values (which are lists) to contain equal number of elements
    # This number will be equal to the max read length:
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


# Function that calculates average quality per base
def calculate_mean_quality_per_base(qualities_per_base):
    import numpy as np
    # The new empty dictionary is created:
    dict_mean_qual = dict()
    # We use 'quality_per_base' dictionary. It is a result of previous function,
    # containing the positions(bp) as the keys and the lists of quality scores in each position as values
    for key in qualities_per_base:
        # For each key we calculate average of value-list (we do not include None, while calculating mean):
        dict_mean_qual[key-1] = np.mean([i for i in qualities_per_base[key] if i is not None])
    return dict_mean_qual


# Function for plot drawing
def plot_per_base_seq_quality(qualities_per_base, dict_mean_qual):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    # We turn 'qualities_per_base' dictionary into pandas data frsme
    base = pd.DataFrame.from_dict(qualities_per_base)
    # Set suitable size:
    plt.figure(figsize=(15, 10))
    # Set default theme and remove margins:
    sns.set_theme()
    plt.margins(0)
    # Create colored background:
    plt.axhspan(0, 20, facecolor='red', alpha=0.1)
    plt.axhspan(20, 28, facecolor='yellow', alpha=0.1)
    plt.axhspan(28, 42, facecolor='green', alpha=0.1)
    # And, finally, draw boxplot
    sns.boxplot(data=base, color="yellow", whis=[10, 90], showfliers=False,
                medianprops=dict(linewidth=0.5, color='red'))
    # We process dictionary in order to set x and y for the plot:
    lists = sorted(dict_mean_qual.items())  # sorted by key, return a list of tuples
    x2, y2 = zip(*lists)  # unpack a list of pairs into two tuples
    # And also draw the mean plot:
    plt.plot(x2, y2, color="blue", linewidth=0.5)
    # Some improvements to make the plot easier to read:
    plt.yticks(np.arange(0, 42, step=2), size=8)
    plt.xticks(np.arange(0, len(qualities_per_base)+2, step=2), size=7)
    plt.title('Quality scores across all bases', size=8)
    plt.xlabel('Position in read', size=8)
    plt.show()
    # plt.savefig('Per_base_sequence_quality.png')


# 2. Per sequence quality scores


# Function that calculate mean quality score per sequence:
def mean_quality(qual):
    quality_score = 0
    #  For each symbol in quality score line, we calculate its number in ascii
    #  table and subtract 33
    for i in qual:
        quality_score += ord(i) - 33
    # Returns round value of mean quality score:
    return round(quality_score / len(qual))


def per_sequence_quality(reads):
    # We process the 3rd element of 'reads' list - the list of quality scores sequences
    quality = reads[3]
    # The new empty dictionary for average quality scores (keys)
    # and the number of sequences with that average (values):
    qual_and_numbers = dict()
    for qual_seq in quality:
        # For each quality score sequence in the list we calculate the average quality score:
        n = mean_quality(qual_seq)
        if n in qual_and_numbers:
            # If the quality score is in the dictionary, we add 1 to its value
            qual_and_numbers[n] += 1
        else:
            # If it is not, we set 1 its value
            qual_and_numbers[n] = 1
    return qual_and_numbers


# Function for plot drawing
def plot_per_seq_quality_scores(d):
    import matplotlib.pyplot as plt
    import numpy as np
    # We process dictionary in order to set x and y for the plot:
    lists = sorted(d.items())  # sorted by key, return a list of tuples
    x, y = zip(*lists)  # unpack a list of pairs into two tuples
    # Set suitable size:
    plt.figure(figsize=(15, 10))
    # And also draw the plot:
    plt.plot(x, y, color="red", label='Average quality per read')
    # Some improvements to make the plot easier to read:
    plt.xlabel('Mean Sequence Quality (Pherd Score)', size=8)
    plt.title('Quality score distributions over all sequences', size=8)
    plt.xticks(np.arange(0, 40, step=2), size=8)
    plt.yticks(np.arange(0, max(d.values())+10000, step=10000), size=8)
    plt.legend(loc='upper right')
    plt.grid(alpha=0.5)
    plt.show()
    # plt.savefig('Per_sequence_quality_scores.png')


# 3. Per base sequence content


# Function that creates and returns a list of 4 dictionaries (1 for each nucleotide)
# that contain number of base pair in the sequence as a key
# and the proportion of this nucleotide as a value
def per_base_nucleotides_proportion(reads):
    # We process the 1st element in the 'reads' list - seq
    seqs = reads[1]
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


# Function for plot drawing
def plot_per_base_seq_content(lst_proportions):
    import matplotlib.pyplot as plt
    import numpy as np
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
    # Set the size
    plt.figure(figsize=(18, 10))
    # And draw the plot for each nucleotide:
    plt.plot(x, ay, color="red", linewidth=0.5, label='%A')
    plt.plot(x, ty, color="blue", linewidth=0.5, label='%T')
    plt.plot(x, gy, color="green", linewidth=0.5, label='%G')
    plt.plot(x, cy, color="black", linewidth=0.5, label='%C')
    # Some improvements to make the plot easier to read:
    plt.xticks(np.arange(0, len(x)+2, step=2), size=7)
    plt.yticks(np.arange(0, 110, step=10), size=8)
    plt.xlabel('Position in read(bp)', size=8)
    plt.title('Sequence content across all bases', size=8)
    plt.legend(loc='upper right')
    plt.grid(alpha=0.5)
    plt.show()
    # plt.savefig('Per_base_sequence_content.png')


# Sample input
input_fastq = 'amp_res_1.fastq'
reads = process_file(input_fastq)
# 1
qualities_per_base = calculate_quality_per_base(reads)
dict_mean_qual = calculate_mean_quality_per_base(qualities_per_base)
plot_per_base_seq_quality(qualities_per_base, dict_mean_qual)
# 2
d = per_sequence_quality(reads)
plot_per_seq_quality_scores(d)
# 3
lst_proportions = per_base_nucleotides_proportion(reads)
plot_per_base_seq_content(lst_proportions)
