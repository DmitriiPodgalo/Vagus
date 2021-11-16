# run the .py file with flag -h in terminal for help (<python parsing_report.py -h>)
import logging
import datetime
from types import new_class

import typer
import jinja2

import fastqc
import FastQC_gorbonos


app = typer.Typer()


@jinja2.pass_context
def format_datetime(context, value):
    return value.strftime('%A')


def prepair_data(input):
    parsed_file = fastqc.read_file(input)

    sequence_length_distribution_result = fastqc.sequence_length_distribution(parsed_file)
    overrepresented_sequences_result = fastqc.overrepresented_sequences(parsed_file)
    adapter_content_result = fastqc.adapter_content(parsed_file)

    context = {'now': datetime.datetime.utcnow(),
               'file': input,
               'name': 'Yulia',
               'captions': ['intro', 'main part', 'results', 
               sequence_length_distribution_result, overrepresented_sequences_result, adapter_content_result]}

    return context


def render_report(context, template, outdir):
    
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

    context = prepair_data(input)

    render_report(context, template, outdir)

    logging.basicConfig(level=getattr(logging, log_level.upper()))
    logging.info('generate report')
    logging.info('report template: %s', template)


def main():
    app()


if __name__ == '__main__':
    main()