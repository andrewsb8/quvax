# QuVax
### mRNA design guided by folding potential

## Getting Started

Create Conda environment

### Mac

```bash conda env create --file requirements.txt```

### Linux

Create and activate environment

```
$ conda env create -n quvax python=3.11
$ conda activate quvax
```

Install dependencies

```
$ python -m pip install -r requirements.txt

OR

$ conda install --file requirements.txt -c conda-forge
```

*The latter currently does not work.

### After Setup

Example execution:

```
$ python design.py -i examples/spike_trim.fasta
````

For help:

```$ python design.py -h```

Testing:

```$ coverage run -m pytest``` or ```$ pytest```

For running a specific set of tests:

```pytest -k "test file name without .py"``` (no path required)

### List of Dependencies in requirements.txt

```
argparse==1.4.0
biopython==1.81
dwave-neal==0.6.0
flake8==6.1.0
pandas==2.0.3
pip==23.2.1
pytest-cov==4.1.0
python-codon-tables==0.1.12
tensorflow==2.13.0
tensorflow-probability==0.21.0
```
