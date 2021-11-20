import matplotlib.pyplot as plt
from scipy import stats
import numpy as np


def read_file(file_path):
    parsed_file = []
    temp_list = []

    with open(file_path) as inf:

        for line in inf:
            temp_list.append(line.rstrip())

            if len(temp_list) == 4:
                parsed_file.append(temp_list)
                temp_list = []

    return parsed_file


def kde(x, steps):
    return np.array([np.sum((x >= steps[i]) & (x <= steps[i + 1])) for i in range(len(steps) - 1)])


def count_gc_content(line):
    length = len(line)
    G_count = line.upper().count('G')
    C_count = line.upper().count('C')
    percent = ((G_count + C_count) / length) * 100
    return percent


def draw_gc_content(parsed_file, DEFAULT_OUTPUT_DIR='./Report_data/'):
    gc_content = [count_gc_content(read[1]) for read in parsed_file]    
    gc_content = np.array(gc_content)
    mean = np.median(gc_content)
    sd = np.std(gc_content)
    xx = np.arange(0, 100, 1)
    y = kde(gc_content, xx)
    fig, ax = plt.subplots(figsize=(8,6))
    ax.plot(xx[:-1], y, color='red', label='GC count per read')

    theoretical_y = stats.norm.pdf(xx[:-1], loc = mean, scale = sd) * len(parsed_file)
    ax.plot(xx[:-1], theoretical_y, color='blue', label='Theoretical Distribution')
    plt.xticks(range(0, 100, 10))
    plt.title('GC distribution over all sequences')
    plt.xlabel('Mean GC content (%)')

    plt.legend()
    plt.grid(alpha = 0.5)
    plt.gcf().set_size_inches(8, 6)
    plt.savefig(DEFAULT_OUTPUT_DIR+'gc_content.png', dpi=100, bbox_inches='tight')
    plt.close()

    total_deviation = np.sum(np.abs(y - theoretical_y)) / len(parsed_file) * 100
    if total_deviation > 30:
        return 'Failure'

    if total_deviation > 15:
        return 'Warning'

    return 'Good'


def draw_N_content(parsed_file, DEFAULT_OUTPUT_DIR='./Report_data/'):
    length_dict = {}
    
    for read in parsed_file:
        length = len(read[1])
        if length in length_dict.keys():
            length_dict[length] += 1
        else:
            length_dict[length] = 1
            
    max_length = max(length_dict.keys())
    N_counter = [0 for i in range(max_length)]
    Read_counter = [0 for i in range(max_length)]
    for read in parsed_file:
        for i in range(len(read[1])):
            if read[1][i] == 'N':
                N_counter[i] += 1
                Read_counter[i] += 1
            else:
                Read_counter[i] += 1

    N_content = [N_counter[i] / Read_counter[i] for i in range(len(N_counter))]    
    N_content = np.array(N_content)
    
    fig, ax = plt.subplots(figsize=(8,6))
    ax.plot(range(len(N_content)), N_content * 100, color='red', label='%N')
    plt.yticks(range(0, 100, 10))
    plt.title('N content across all bases')
    plt.xlabel('Position in read (bp)')

    plt.legend()
    plt.grid(alpha = 0.5)
    plt.gcf().set_size_inches(8, 6)
    plt.savefig(DEFAULT_OUTPUT_DIR+'N_content.png', dpi=100, bbox_inches='tight')
    plt.close()

    max_content = np.max(N_content * 100) 
    if max_content > 20:
        return 'Failure'
    
    if max_content > 5:
        return 'Warning'

    return 'Good'


def draw_deduplicated(parsed_file, DEFAULT_OUTPUT_DIR='./Report_data/'):
    reads_count = {}
    for read in parsed_file:
        if read[1] in reads_count:
            reads_count[read[1]] += 1
        else:
            reads_count[read[1]] = 1
    
    counts = [a for a in reads_count.values()]
    total_number = len(parsed_file)
    distinct_number = len(reads_count)
    
    percent = round(100 * distinct_number / total_number, 2)
    
    lower = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 50, 100, 500, 1000, 5000, 10000]
    upper = lower[1:] + [1e10]
    
    y_sum = []
    y_len = []
    
    for i in range(len(lower)):
        low = lower[i]
        hight = upper[i]
        values = [a for a in counts if low <= a < hight]
        y_sum.append(sum(values))
        y_len.append(len(values))
        
    fig, ax = plt.subplots(figsize=(8,6))
    x_values = range(len(lower))
    x_ticks = map(str, lower)
    y_sum = [100 * x / total_number for x in y_sum]
    y_len = [100 * x / distinct_number for x in y_len]

    plt.plot(y_len, color='red', label='Deduplicated sequences')
    plt.plot(y_sum, color='blue', label='Total sequences')
    plt.xticks(x_values, x_ticks)
    plt.yticks(range(0, 100, 10))
    plt.title('Percent of seq remaining if deduplicated ' + str(percent) + '%')
    plt.xlabel('Sequence duplication level')
    plt.legend()
    plt.grid(alpha = 0.5)
    plt.gcf().set_size_inches(8, 6)
    plt.savefig(DEFAULT_OUTPUT_DIR+'deduplication.png', dpi=100, bbox_inches='tight')
    plt.close()


    counts = np.array(counts)
    non_unique_seq_frac = np.sum(counts[counts > 1]) / total_number * 100

    if non_unique_seq_frac > 50:
        return 'Failure'

    if non_unique_seq_frac > 20:
        return 'Warning'
    
    return 'Good'


def main():
    parsed_file = read_file('23769689.fastq')
    draw_gc_content(parsed_file)
    draw_N_content(parsed_file)
    draw_deduplicated(parsed_file)


if __name__ == '__main__':
    main()
