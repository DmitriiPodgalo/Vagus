import FastQC_G


def sequence_length(parsed_file):
    '''
    Return min and max read length
    '''
    length_reads = [len(read[1]) for read in parsed_file]
    return min(length_reads), max(length_reads)


def count_all_GC(parsed_file):
    '''
    Count GC content over all sequenses.
    '''
    GC_count, line_lenght = int(), int()

    for read in parsed_file:
        GC_lenght = FastQC_G.count_gc(read[1])
        GC_count += GC_lenght[0]
        line_lenght += GC_lenght[1]

    return round(GC_count * 100 / line_lenght, 1)


def encoding_detector(parsed_file):
    '''
    We calculate the minimum and maximum ASCII quality in all reads. 
    Based on this, we choose the encoding.
    '''
    min_quality, max_qulity = 0, 0

    for read in parsed_file:
        qulity_line = [ord(quality_element) for quality_element in read[3]]
        min_qual_line, max_qual_line = min(qulity_line), max(qulity_line)

        if min_qual_line < min_quality or min_quality == 0:
            min_quality = min_qual_line
        if max_qual_line > max_qulity or max_qulity == 0:
            max_qulity = min_quality

    if min_quality < 59 and min_quality >= 33 and max_qulity <= 74:
        return 'Phred+33'

    elif min_quality >= 64 and max_qulity > 73 and max_qulity <= 104:
        return 'Phred+64'

    elif min_quality >= 59 and min_quality < 64 and max_qulity > 73:
        return 'Solexa+64'

    elif min_quality >= 33 and max_qulity <= 104:
        return 'PacBio'

    else:
        return 'Unknown score encoding!'
