# run the .py file with flag -h in terminal for help (<python parsing_report.py -h>)
import argparse
import fastqc


def main(input_fastq, output_dir = './Results/'):
    parsed_file = fastqc.read_file(input_fastq)

    sequence_length_distribution_result = fastqc.sequence_length_distribution(parsed_file)
    overrepresented_sequences_result = fastqc.overrepresented_sequences(parsed_file)
    adapter_content_result = fastqc.adapter_content(parsed_file)


if __name__ == '__main__':
    """
    Show help message: -h, -- help
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", help="path to fastq file")
    parser.add_argument("-o", "--outdir", help="path to the folder to which the result will be written")
    args = parser.parse_args()

    main(input_fastq=args.input,
         output_dir=args.outdir)
