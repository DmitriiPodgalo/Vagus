## Vagus. FASTQ quality analyzer tool

Tool is able to:
1. ...
2. ...



# Install and run with pip
## Installation

```console
git clone https://github.com/DmitriiPodgalo/Vagus.git

# Create and activate your virtual environment

# create virtual environment
python3.10 -m venv ./venv

# activate virtual environment
source ./venv/bin/activate

# if you install it not from main or master, change branch
git checkout branch_name

# required by pip to build wheels
pip install wheel==0.37.0 

# Install requirements
pip install -r ./requirements.txt
```

## Run file
```console
python3.10 parsing_report.py --input "input_file.fastq"
```

# Install and run with poetry (Ubuntu)
```console
# install poetry
# for details look for https://python-poetry.org/docs/
sudo apt-get install curl
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python3.10 -

# poetry will be accessible in current session
source $HOME/.poetry/env

# prepare project
git clone https://github.com/DmitriiPodgalo/Vagus.git

# if you install it not from main or master, change branch
git checkout branch_name

# Install requirements
poetry env use python3.10
poetry install

# Run
poetry run python parsing_report.py --input "input_file.fastq"
```


# Install and run with poetry (Windows)
```console
## run PowerShell as administrator!

# install poetry (Poetry 1.1.11)
# for details look for https://python-poetry.org/docs/
(Invoke-WebRequest -Uri https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py -UseBasicParsing).Content | python -

## close the PowerShell administrator mode.

# prepare project
git clone https://github.com/DmitriiPodgalo/Vagus.git

# if you install it not from main or master, change branch
git checkout branch_name

# Install requirements
poetry env use python
poetry install

# Run
poetry run python parsing_report.py --input "input_file.fastq"