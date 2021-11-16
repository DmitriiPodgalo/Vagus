# 1 и 3 функция возвращают check и записывают рисунок
# 2 функция возвращает check и записывает csv файл
# сначала вызываем read_file() и передаем путь до fastq файла
# сохраняем в переменную и передаем ее для остальных функций

from collections import Counter
import matplotlib.pyplot as plt
import csv


def sequence_length_distribution(parsed_file):
    lines_reads = [i[1] for i in parsed_file]
    seq_dict = Counter(sorted(map(len, lines_reads)))
    counter = 1

    if len(seq_dict) == 1:
        key = list(seq_dict.keys())[0]
        value = list(seq_dict.values())[0]
        xs = [key - 1, key, key + 1]
        ys = [0, value, 0]
    else:
        xs = seq_dict.keys()
        ys = seq_dict.values()

    while True:
        if (max(xs) - min(xs) + 1) // counter < 10:
            xticks = (list(range(min(xs), max(xs) + 1, counter)))
            break
        counter += 1

    color = '#CB382E'
    label = 'Sequence length'
    xlabel = 'Sequence Length (bp)'
    title = 'Distribution of sequence lengths over all sequences'
    yticks = None

    draw_plot(xs, ys, color, label, xlabel, title, xticks, yticks)

    # сохраняем картинку и возвращаем check

    plt.savefig('sequence_length_distribution.png', bbox_inches='tight')
    plt.close()

    if seq_dict[0] != 0:
        return 'failure'
    elif len(seq_dict) != 1:
        return 'warning'
    else:
        return 'good'


def overrepresented_sequences(parsed_file):
    lines_reads = [i[1] for i in parsed_file]
    number_reads = len(lines_reads)
    counter = 0

    d = Counter([i[:50] if len(i) > 75 else i for i in lines_reads])
    d = list(map(list, sorted(d.items(), key=lambda x: x[1], reverse=True)))

    for i in range(0, len(d)):
        percents = d[i][1] / number_reads

        if percents < 0.001:
            d = d[:counter]
            break

        d[i].append(percents * 100)
        counter += 1

    if len(d) == 0:
        return 'good'

    # сохраняем таблицу и возвращаем check

    with open('or_seq.csv', 'w') as f:
        writer = csv.writer(f)
        for read in d:
            writer.writerow(read)

    if d[0][2] > 1:
        return 'failure'
    else:
        return 'warning'


def adapter_content(parsed_file):
    lines_reads = [i[1] for i in parsed_file]
    adapters = {'Illumina Universal Adapter': 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT',
                'Illumina Small RNA 3\' Adapter': 'TGGAATTCTCGGGTGCCAAGG',
                'Illumina Small RNA 5\' Adapter': 'GUUCAGAGUUCUACAGUCCGACGAUC',
                'Nextera Transposase Sequence 1': 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG',
                'Nextera Transposase Sequence 2': 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG',
                'SOLID Small RNA Adapter': 'CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT'}

    adap_check = {'Illumina Universal Adapter': [],
                  'Illumina Small RNA 3\' Adapter': [],
                  'Illumina Small RNA 5\' Adapter': [],
                  'Nextera Transposase Sequence 1': [],
                  'Nextera Transposase Sequence 2': [],
                  'SOLID Small RNA Adapter': []}

    colors = {'Illumina Universal Adapter': '#D14139',
              'Illumina Small RNA 3\' Adapter': '#1D2DD8',
              'Illumina Small RNA 5\' Adapter': '#73DF57',
              'Nextera Transposase Sequence 1': 'purple',
              'Nextera Transposase Sequence 2': '#6D6D6D',
              'SOLID Small RNA Adapter': '#EC69F8'}
    max_val = 0
    min_key = max(map(len, lines_reads))
    max_key = 0

    for read in lines_reads:
        for adapter in adapters:
            in1 = read.upper().rfind(adapters[adapter]) + 1
            in2 = read.upper().rfind(adapters[adapter][::-1].translate(str.maketrans('ATCG', 'TAGC'))) + 1

            if in1 != 0:
                adap_check[adapter] += range(in1, len(read) + 1)
            elif in2 != 0:
                adap_check[adapter] += range(in2, len(read) + 1)

    for adap in adap_check:
        c_dict = Counter(adap_check[adap])
        adap_check[adap] = {k: v for k, v in sorted(c_dict.items(), key=lambda x: x[0])}

    for keys, vals in adap_check.items():
        for key, val in vals.items():
            if val > max_val:
                max_val = val
            if key > max_key:
                max_key = key
            if key < min_key:
                min_key = key

    thre = max_val / len(lines_reads)

    for key, vals in adap_check.items():
        if len(vals) == 0:
            if max_key == 0:
                xs = range(1, min_key + 1)
                ys = [0] * min_key
            else:
                xs = range(1, max_key + 1)
                ys = [0] * max_key
        else:
            xs = list(range(1, min_key)) + list(vals.keys())
            ys = [0] * (min_key - 1) + [i / len(lines_reads) * 100 for i in list(vals.values())]

        color = colors[key]
        label = key
        xlabel = 'Position in read (bp)'
        title = '% Adapter'
        xticks = None
        yticks = range(0, 101, 10)

        draw_plot(xs, ys, color, label, xlabel, title, xticks, yticks)

    # сохраняем картинку и возвращаем check

    plt.savefig('adapter_content.png', bbox_inches='tight')
    plt.close()

    if thre > 0.1:
        return 'failure'
    elif thre > 0.05:
        return 'warning'
    else:
        return 'good'


def read_file(file_path):
    parsed_file = []
    temp_list = []

    with open(file_path) as inf:

        for line in inf:
            temp_list.append(line.rstrip())

            if len(temp_list) == 4:
                parsed_file.append(temp_list)
                temp_list = []

        if len(temp_list) != 0:
            raise Exception('Invalid number of file\'s lines')

    return parsed_file


def draw_plot(xs, ys, color, label, xlabel, title, xticks, yticks):

    plt.plot(xs, ys, color=color, label=label)

    plt.xticks(xticks)
    plt.yticks(yticks)

    plt.xlabel(xlabel, size=8)
    plt.title(title, size=8)

    plt.grid(alpha=0.5)
    plt.legend()
