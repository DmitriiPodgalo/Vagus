# Vagus. FASTQ quality analyzer tool

The tool allows one to analyze the quality of sequencing reads.

### The Vagus report contains the following data:
0. **Basic statistics** - file encoding, number of reads, min and max read lenght, GC% in all reads.
1. **Per base sequence quality** - plots out the range of quality values across all bases at each position in the FastQ file. 
    - Warning - if the lower quartile for any base is less than 10, or if the median for any base is less than 25.
    - Failure - if the lower quartile for any base is less than 5 or if the median for any base is less than 20.
2. **Per sequence quality scores** - plots out the quality scor distribution over all sequences. 
    - Warning - if the most frequently observed mean quality is below 27.
    - Failure - if the most frequently observed mean quality is below 20.
3. **Per base sequence content** - plots out the proportion of each base position in a .fastq for which each of the four normal DNA bases has been called.
    - Warning - if the difference between A and T, or G and C is greater than 10% in any position.
    - Failure - if the difference between A and T, or G and C is greater than 20% in any position.
4. **Per sequence GC content** - plots out the GC content across the whole length of each sequence in a file and compares it to a modelled normal distribution of GC content.
    - Warning - if the sum of the deviations from the normal distribution represents more than 15% of the reads.
    - Failure - if the sum of the deviations from the normal distribution represents more than 30% of the reads.
5. **Per base N content** - plots out the percentage of base calls at each position for which an N was called.
    - Warning - if any position shows an N content of >5%.
    - Failure - if any position shows an N content of >20%.
6. **Sequence length distribution** - generate plot showing the distribution of fragment sizes in the file which was analysed.
    - Warning - if all sequences are not the same length.
    - Failure - if any of the sequences have zero length.
7. **Sequence duplication levels** -  counts the degree of duplication for every sequence in a library and creates a plot showing the relative number of sequences with different degrees of duplication.
    - Warning - if non-unique sequences make up more than 20% of the total.
    - Failure - if non-unique sequences make up more than 50% of the total.
8. **Overrepresented sequences** - lists all of the sequence which make up more than 0.1% of the total and write it to .csv file.
    - Warning - if any sequence is found to represent more than 0.1% of the total.
    - Failure - if any sequence is found to represent more than 1% of the total.
9. **Adapter content** - plots a cumulative percentage count of the proportion of library which has seen each of the adapter sequences at each position. 
    - Warning - if any sequence is present in more than 5% of all reads.
    - Failure - if any sequence is present in more than 10% of all reads.

# Try it!
The program analyzes only sequencing reads in the **.fastq** format, not .gz archives.

In addition to the required .fastq file, a directory path for saving analysis results can be transferred to the program input. Otherwise, the default prefix will be used.

**Important!** The time at which the report was generated will be added to output directory prefix.

As a quality check result of the program's work, a **.html** file will be generated. The **.html** report will work correctly only together with the generated .png and .csv files from the output folder.

You can test the Vagus with example .fastq file in `./Test_data` directory. The Example report placed in `./Example_report`.

See below for a description of how to install using pip or poetry and how to get started.

While the program is running, the progress of the work is displayed in the console.

# Requiered dependencies
Required Python version [3.9-3.10].

All dependencies can be installed via pip or poetry.
``` console
Jinja2 = "3.0.2"
typer = "^0"
matplotlib = "3.5.0"
scipy="1.7.2"
numpy="1.21.4"
seaborn="0.11.2"
pandas="1.3.4"
click=="8.0.3"
```

## Install and run with pip (Ubuntu)
```console
git clone https://github.com/DmitriiPodgalo/Vagus.git
cd Vagus/

# Create and activate your virtual environment

# Install venv
sudo apt install python3-venv

# create virtual environment (with python 3.9 or 3.10)
python3 -m venv ./venv

# activate virtual environment
source ./venv/bin/activate

# required by pip to build wheels
pip install wheel==0.37.0 

# Install requirements
pip install -r ./requirements.txt

# Run file (with python 3.9 or 3.10)
python3 parsing_report.py --help
python3 parsing_report.py -i ./Test_data/amp_res_1.fastq -o ./results_dir/

# Exit fron venv:
deactivate
```

## Install and run with poetry (Ubuntu)
```console
# install poetry
# for details look for https://python-poetry.org/docs/
sudo apt-get install curl
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python3.10 -

# poetry will be accessible in current session
source $HOME/.poetry/env

# prepare project
git clone https://github.com/DmitriiPodgalo/Vagus.git
cd Vagus/

# Install requirements
poetry env use python3
poetry install 

# Run file (with python 3.9 or 3.10)
poetry run python parsing_report.py --help
poetry run python parsing_report.py -i ./Test_data/amp_res_1.fastq -o ./results_dir/
```

## Install and run with pip (Windows)
```console
git clone https://github.com/DmitriiPodgalo/Vagus.git
cd Vagus/

# Chech python version with 
python --version

# Create and activate your virtual environment

# create virtual environment (with python 3.9 or 3.10)
python -m venv ./venv

# activate virtual environment
.\venv\Scripts\activate

# If the command above does not work, please enable scripting and resubmit the command above after it.
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser

# required by pip to build wheels
pip install wheel==0.37.0 

# Install requirements
pip install -r ./requirements.txt

# Run file (with python 3.9 or 3.10)
python parsing_report.py --help
python parsing_report.py -i ./Test_data/amp_res_1.fastq -o ./results_dir/

# Exit fron venv:
deactivate
```

## Install and run with poetry (Windows)
```console
## run PowerShell as administrator!

# install poetry
# for details look for https://python-poetry.org/docs/
(Invoke-WebRequest -Uri https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py -UseBasicParsing).Content | python -

# Exit the PowerShell administrator mode.
# Restart command prompt as user.

# prepare project
git clone https://github.com/DmitriiPodgalo/Vagus.git

# Install requirements
poetry env use python
poetry install

# Run
poetry run python parsing_report.py --help
poetry run python parsing_report.py -i ./Test_data/amp_res_1.fastq -o ./results_dir/
```

# Authors 
1. [**Y. Burankova**](https://github.com/Freddsle):
    - README;
    - console input parser via typer;
    - html report creation functions;
    - logging;
    - requiremets and .toml files;
    - correction of the plot design.
    - tested Vagus on Ubuntu 20.04 LTS, Python 3.9.5. and Windows 10 Pro 64x 20H2, Python 3.9.4.
2. [**A. Gorbonos**](https://github.com/IlonaGA):
    - draw_gc_content function checker and plot;
    - draw_N_content function checker and plot;
    - draw_deduplicated function checker and plot;
    - tested on MacOS Big Sur v.11.2.3, Python 3.9.8;
3. [**A. Tokareva**](https://github.com/staceyso):
    - plot_per_base_seq_quality function checker and plot;
    - plot_per_seq_quality_scores function checker and plot;
    - plot_per_base_seq_content checker and plot;
    - tested Vagus on WSL2, Ubuntu 20.04 LTS, Python 3.9.5.
3. [**D. Podgalo**](https://github.com/DmitriiPodgalo):
    - tool naming;
    - fastq file pasing function;
    - sequence_length_distribution function checker and plot;
    - overrepresented_sequences function checker and plot;
    - adapter_content function checker and plot;
    - correction of the plot design.
    - tested Vagus on MacOS v11.6(20G165), Python 3.9.7.

# Tested on
Ubuntu 20.04 LTS, Python 3.9.5.\
Windows 10 Pro 64x 20H2, Python 3.9.4.\
MacOS v11.6 (20G165), Python 3.9.7.\
MacOS v11.2.3 (20D91), Python 3.9.8.\
WSL2, Ubuntu 20.04 LTS, Python 3.9.5.
