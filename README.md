## Vagus. FASTQ quality analyzer tool

Tool is able to:
1. ...
2. ...

Required Python version [3.9-3.10].

# Requiered dependencies
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
All requrements are installed via pip or poetry.


# Install and run with pip (Ubuntu)
```console
git clone https://github.com/DmitriiPodgalo/Vagus.git
cd Vagus/

# Create and activate your virtual environment

# create virtual environment (with python 3.9 or 3.10)
python3 -m venv ./venv

# activate virtual environment
source ./venv/bin/activate

# if you install it not from main or master, change branch
git checkout branch_name

# required by pip to build wheels
pip install wheel==0.37.0 

# Install requirements
pip install -r ./requirements.txt

# Run file (with python 3.9 or 3.10)
python3 parsing_report.py --help
python3 parsing_report.py --i "input_file.fastq" -o ./out_dir/
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
cd Vagus/

# if you install it not from main or master, change branch
git checkout branch_name

# Install requirements
poetry env use python3
poetry install 

# Run file (with python 3.9 or 3.10)
poetry run parsing_report.py --help
poetry run parsing_report.py --i "input_file.fastq" -o ./out_dir/
```

# Install and run with pip (Windows)
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

# if you install it not from main or master, change branch
git checkout branch_name

# required by pip to build wheels
pip install wheel==0.37.0 

# Install requirements
pip install -r ./requirements.txt

# Run file (with python 3.9 or 3.10)
python parsing_report.py --help
python parsing_report.py --i "input_file.fastq" -o ./out_dir/
```

# Install and run with poetry (Windows)
```console
## run PowerShell as administrator!

# install poetry
# for details look for https://python-poetry.org/docs/
(Invoke-WebRequest -Uri https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py -UseBasicParsing).Content | python -

# Exit the PowerShell administrator mode.
# Restart command prompt as user.

# prepare project
git clone https://github.com/DmitriiPodgalo/Vagus.git

# if you install it not from main or master, change branch
git checkout branch_name

# Install requirements
poetry env use python
poetry install

# Run
poetry run parsing_report.py --help
poetry run parsing_report.py --i "input_file.fastq" -o ./out_dir/
```

# Tested on
Ubuntu 20.04 LTS, Python 3.9.5.
Windows 10 Pro 64x 20H2, Python 3.9.4.