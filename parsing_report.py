import logging
import datetime
from types import new_class

import typer
import jinja2

import fastqc
import FastQC_functions
import FastQC_G


app = typer.Typer()


@jinja2.pass_context
def format_datetime(context, value):
    return value.strftime('%d/%m/%y  %H:%M:%S')


def prepair_data(input, outdir):
    '''
    Parsed file and create quality checks and plots.
    Save plots to "outdit/Report_data/".
    '''
    parsed_file = fastqc.read_file(input)

    # 1-3
    sequence_length_distribution_result = fastqc.sequence_length_distribution(parsed_file, outdir)
    overrepresented_sequences_result = fastqc.overrepresented_sequences(parsed_file, outdir)
    adapter_content_result = fastqc.adapter_content(parsed_file, outdir)
    # 4
    qualities_per_base = FastQC_functions.calculate_quality_per_base(parsed_file)
    dict_mean_qual = FastQC_functions.calculate_mean_quality_per_base(qualities_per_base)
    FastQC_functions.plot_per_base_seq_quality(qualities_per_base, dict_mean_qual, outdir)
    # 5
    d = FastQC_functions.per_sequence_quality(parsed_file)
    FastQC_functions.plot_per_seq_quality_scores(d, outdir)
    # 6
    lst_proportions = FastQC_functions.per_base_nucleotides_proportion(parsed_file)
    FastQC_functions.plot_per_base_seq_content(lst_proportions, outdir)
    # 7
    gc_content_result = FastQC_G.draw_gc_content(parsed_file, outdir)
    # 8
    N_content_result = FastQC_G.draw_N_content(parsed_file, outdir)
    #6
    deduplicated_result = FastQC_G.draw_deduplicated(parsed_file, outdir)

    context = {'now': datetime.datetime.utcnow(),
               'file': input,
               'name': 'Yulia',
               'captions': ['intro', 'main part', 'results']}

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
    environment.globals |= {}

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