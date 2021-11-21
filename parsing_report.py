import logging
import datetime
import os, shutil
import re
import pandas as pd

from types import new_class
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
        shutil.copytree('./Report_data/check_img/', outdir + 'check_img/', dirs_exist_ok=True)


def prepair_data(input, outdir):
    '''
    Parsed file and create quality checks and plots.
    Save plots to "outdit/Report_data/".
    '''
    parsed_file = fastqc.read_file(input)
    prepare_outdir(outdir)
    
    # 1-3
    sequence_length_distribution_result = fastqc.sequence_length_distribution(parsed_file, outdir)
    overrepresented_sequences_result = fastqc.overrepresented_sequences(parsed_file, outdir)
    adapter_content_result = fastqc.adapter_content(parsed_file, outdir)
    # 4
    parsed_file_for_func = FastQC_functions.process_file(input)

    qualities_per_base = FastQC_functions.calculate_quality_per_base(parsed_file_for_func)
    dict_mean_qual = FastQC_functions.calculate_mean_quality_per_base(qualities_per_base)
    FastQC_functions.plot_per_base_seq_quality(qualities_per_base, dict_mean_qual, outdir)
    # 5
    d = FastQC_functions.per_sequence_quality(parsed_file_for_func)
    FastQC_functions.plot_per_seq_quality_scores(d, outdir)
    # 6
    lst_proportions = FastQC_functions.per_base_nucleotides_proportion(parsed_file_for_func)
    FastQC_functions.plot_per_base_seq_content(lst_proportions, outdir)
    # 7
    gc_content_result = FastQC_G.draw_gc_content(parsed_file, outdir)
    # 8
    N_content_result = FastQC_G.draw_N_content(parsed_file, outdir)
    #6
    deduplicated_result = FastQC_G.draw_deduplicated(parsed_file, outdir)

    # Basic statistics
    Encoding = FastQC_B.encoding_detector(parsed_file)
    total_sequences = len(parsed_file)
    sequence_length = [len(min(parsed_file_for_func[1])), len(max(parsed_file_for_func[1]))]
    mean_GC = FastQC_B.count_all_GC(parsed_file)

    if os.path.exists(outdir+'or_seq.csv'):
        overrepresented_sequences_table = pd.read_csv(outdir+'or_seq.csv',
                                                      names=["Sequence", "Count", "Percentage"])
        overrepresented_sequences_table = overrepresented_sequences_table.to_html(index = False,
                                                                                  justify = 'center',
                                                                                  classes = 'table_dupl')
    else:
        overrepresented_sequences_table = 'No overrepresented sequences'
    
    if sequence_length[0] == sequence_length[1]:
        seq_length = sequence_length[0]
    else:
        seq_length = str(sequence_length[0], '-', sequence_length[1])

    # context for html report
    context = {'now': datetime.datetime.utcnow(),
               'file': input,
               'outdir': outdir,
               'Encoding' : Encoding,
               'total_sequences': total_sequences,
               'sequence_length': seq_length,
               'GC' : mean_GC,
               'gc_content_result' : gc_content_result,
               'N_content_result' : N_content_result,
               'sequence_length_distribution_result': sequence_length_distribution_result,
               'deduplicated_result' : deduplicated_result,
               'overrepresented_sequences_result' : overrepresented_sequences_result,
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

    logging.info('report created')


def check_outdir(outdir):
    '''
    Check the last simbol in outdir.
    If in not the '/' - added it to path.
    '''
    pattern_end = re.compile(r".*/$")
    pattern_start = re.compile(r"\./.*/$")

    if not pattern_end.match(outdir):
        outdir = outdir + '/'
    
    if not pattern_start.match(outdir):
        outdir = './' + outdir
    
    return outdir


DEFAULT_TEMPLATE = './Report_templates/report.html.j2'
DEFAULT_OUTPUT_DIR = './Report_data/'


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

    outdir = check_outdir(outdir)

    context = prepair_data(input, outdir)
    render_report(context, template, outdir)

    logging.basicConfig(level=getattr(logging, log_level.upper()))
    logging.info('generate report')
    logging.info('report template: %s', template)


def main():
    app()


if __name__ == '__main__':
    main()