import logging
import datetime
from types import new_class
import os

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


def prepair_data(input, outdir):
    '''
    Parsed file and create quality checks and plots.
    Save plots to "outdit/Report_data/".
    '''
    parsed_file = fastqc.read_file(input)
    
    if not os.path.exists(outdir):
            os.makedirs(outdir)
    
    # 1-3
    sequence_length_distribution_result = fastqc.sequence_length_distribution(parsed_file, outdir)
    overrepresented_sequences_result = fastqc.overrepresented_sequences(parsed_file, outdir)
    adapter_content_result = fastqc.adapter_content(parsed_file, outdir)
    # 4
    parsed_file_for_func = FastQC_functions.process_file(input)

    qualities_per_base = FastQC_functions.calculate_quality_per_base(parsed_file_for_func)
    dict_mean_qual = FastQC_functions.calculate_mean_quality_per_base(qualities_per_base)
    FastQC_functions.plot_per_base_seq_quality(qualities_per_base, dict_mean_qual)
    # 5
    d = FastQC_functions.per_sequence_quality(parsed_file_for_func)
    FastQC_functions.plot_per_seq_quality_scores(d)
    # 6
    lst_proportions = FastQC_functions.per_base_nucleotides_proportion(parsed_file_for_func)
    FastQC_functions.plot_per_base_seq_content(lst_proportions)
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
               'adapter_content_result': adapter_content_result}

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


DEFAULT_TEMPLATE = './Report_templates/report.html.j2'
DEFAULT_OUTPUT_DIR = './Report_data/'


@app.command()
def generate(input: str = typer.Option(...),
             template: str = DEFAULT_TEMPLATE,
             outdir: str = DEFAULT_OUTPUT_DIR,
             log_level: str = 'info'):

    context = prepair_data(input, outdir)
    render_report(context, template, outdir)

    logging.basicConfig(level=getattr(logging, log_level.upper()))
    logging.info('generate report')
    logging.info('report template: %s', template)


def main():
    app()


if __name__ == '__main__':
    main()