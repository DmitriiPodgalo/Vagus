import logging
import datetime
import os
import shutil
import re
import pandas as pd

from pathlib import Path
import typer
import jinja2

import fastqc
import FastQC_functions
import FastQC_G
import FastQC_B


app = typer.Typer()


@jinja2.pass_context
def format_datetime(context, value):
    return value.strftime('%d/%m/%y  %H:%M:%S')


def directory_datetime():
    '''
    Return NOW string for directory name.
    '''
    return datetime.datetime.now().strftime('%Y%b%d_%H_%M/')


def raise_error(text):
    raise Exception(text)


def prepare_outdir(outdir):
    '''
    Check if outdir exists, if not - create it.
    Also check if images for status chesk exists, if not - create it.
    '''
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if not os.path.exists(outdir + 'check_img/error.jpg'):
        shutil.copytree('./Report_templates/check_img/', outdir + 'check_img/')


def prepair_data(input, outdir):
    '''
    Parsed file and create quality checks and plots.
    Save plots to "outdir_DATE_TIME/".
    Print log to console.
    '''
    parsed_file = FastQC_functions.read_file(input)
    prepare_outdir(outdir)
    logging.info('outdir generated')

    # 1
    sequence_length_distribution_result = fastqc.sequence_length_distribution(parsed_file, outdir)
    logging.info('sequence length distribution result generated')
    # 2
    overrepresented_sequences_result = fastqc.overrepresented_sequences(parsed_file, outdir)
    logging.info('overrepresented sequences result generated')
    # 3
    adapter_content_result = fastqc.adapter_content(parsed_file, outdir)
    logging.info('adapter content result generated')
    # 4
    per_base_seq_quality_result = FastQC_functions.plot_per_base_seq_quality(parsed_file, outdir)
    logging.info('per base sequence quality result generated')
    # 5
    per_seq_quality_scores_result = FastQC_functions.plot_per_seq_quality_scores(parsed_file, outdir)
    logging.info('per sequence quality scores result generated')
    # 6
    per_base_seq_content_result = FastQC_functions.plot_per_base_seq_content(parsed_file, outdir)
    logging.info('per base sequence content result generated')
    # 7
    gc_content_result = FastQC_G.draw_gc_content(parsed_file, outdir)
    logging.info('GC content result generated')
    # 8
    N_content_result = FastQC_G.draw_N_content(parsed_file, outdir)
    logging.info('N content result generated')
    # 9
    deduplicated_result = FastQC_G.draw_deduplicated(parsed_file, outdir)
    logging.info('deduplicated generated')

    # Basic statistics
    Encoding = FastQC_B.encoding_detector(parsed_file)
    total_sequences = len(parsed_file)
    sequence_length = FastQC_B.sequence_length(parsed_file)
    mean_GC = FastQC_B.count_all_GC(parsed_file)
    logging.info('basic statusctics generated')

    if os.path.exists(outdir+'or_seq.csv'):
        overrepresented_sequences_table = pd.read_csv(outdir+'or_seq.csv',
                                                      names=["Sequence", "Count", "Percentage"])
        overrepresented_sequences_table = overrepresented_sequences_table.to_html(index=False,
                                                                                  justify='center',
                                                                                  classes='table_dupl')
    else:
        overrepresented_sequences_table = 'No overrepresented sequences'

    if sequence_length[0] == sequence_length[1]:
        seq_length = sequence_length[0]
    else:
        seq_length = str(sequence_length[0])+'-'+str(sequence_length[1])

    input_file_short = re.search(r'\w*\.fastq$', str(input)).group(0)

    # context for html report
    context = {'now': datetime.datetime.utcnow(),
               'file': input_file_short,
               'outdir': outdir,
               'Encoding': Encoding,
               'total_sequences': total_sequences,
               'sequence_length': seq_length,
               'GC': mean_GC,

               'per_base_seq_quality_result': per_base_seq_quality_result,
               'per_seq_quality_scores_result': per_seq_quality_scores_result,
               'per_base_seq_content_result': per_base_seq_content_result,
               'gc_content_result': gc_content_result,
               'N_content_result': N_content_result,
               'sequence_length_distribution_result': sequence_length_distribution_result,
               'deduplicated_result': deduplicated_result,
               'overrepresented_sequences_result': overrepresented_sequences_result,
               'adapter_content_result': adapter_content_result,
               'overrepresented_sequences_table': overrepresented_sequences_table}

    return context


def render_report(context, template, outdir):
    '''
    Create html report.
    Save it to "outdit/Report.html".
    '''
    environment = jinja2.Environment(loader=jinja2.FileSystemLoader('./'),
                                     autoescape=False,
                                     undefined=jinja2.StrictUndefined,
                                     extensions=['jinja2.ext.loopcontrols'])

    environment.filters |= {'format_datetime': format_datetime}
    environment.globals |= {"raise_error": raise_error}

    loaded_template = environment.get_template(template)

    rendered_template = loaded_template.render(context)

    with open(outdir+'Report.html', 'w') as f:
        f.write(rendered_template)


def check_outdir(outdir, now_time):
    '''
    Check the last simbol in outdir.
    If in not the '/' - added it to path.
    '''
    pattern_end = re.compile(r".*/$")
    pattern_start = re.compile(r"\./.*/$")

    if pattern_end.match(outdir):
        outdir = outdir[:-1]

    outdir = outdir + now_time

    if not pattern_start.match(outdir):
        outdir = './' + outdir

    return outdir


now_time = directory_datetime()
DEFAULT_TEMPLATE = './Report_templates/report.html.j2'
DEFAULT_OUTPUT_DIR = './Report_data' + now_time


@app.command()
def generate(input: Path = typer.Option(...,
                                        "--input", "-i",
                                        help="Path to input file from Vagus repository",
                                        exists=True,
                                        file_okay=True,
                                        dir_okay=False,
                                        writable=False,
                                        readable=True,
                                        resolve_path=True),
             outdir: str = typer.Option(DEFAULT_OUTPUT_DIR,
                                        "--outdir", "-o",
                                        help="Path to analysis output directory from Vagus repository"),
             template: str = DEFAULT_TEMPLATE,
             log_level: str = 'info'):

    logging.basicConfig(level=getattr(logging, log_level.upper()))

    outdir = check_outdir(outdir, now_time)
    context = prepair_data(input, outdir)

    logging.info('report template: %s', template)
    logging.info('generate report')
    render_report(context, template, outdir)
    logging.info('report created')


def main():
    app()


if __name__ == '__main__':
    main()
