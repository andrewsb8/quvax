# QuVax
### mRNA design guided by folding potential

## Getting Started

Install your preferred version of conda (Anaconda or Miniconda). Create Conda environment:

```
$ conda create -n quvax python=3.11
$ conda activate quvax
```

Install dependencies:

```
$ python -m pip install -r requirements.txt

OR

$ conda install --file requirements.txt -c conda-forge
```

*The latter currently does not work.

### After Setup

QuVax has two primary functions: codon sequence and folding optimization and a small analysis suite. Optimization is executed with ```design.py``` and analyses with ```analyze.py```.

Example execution of ```design.py```:

```
$ python design.py -i examples/spike_trim.fasta
```

```analyze.py``` then takes in the output of ```design.py```, which is simply a pickled dictionary with default extension ```.qu``` (this can be changed on command line if desired). Example execution of ```analyze.py```:

```
$ python analyze.py -i quvax.qu
```

For help (each will produce unique output):

```$ python design.py -h``` or ```$ python analyze.py -h```

Testing:

```$ coverage run -m pytest``` or ```$ pytest```

Test python files are in ```tests/test_cases```. For running a specific set of tests (no path required):

```pytest -k "test file name without .py"```

### List of Dependencies in requirements.md

```
argparse==1.4.0
biopython==1.81
black==23.12.1
dwave-neal==0.6.0
flake8==6.1.0
pandas==2.0.3
pytest-cov==4.1.0
python-codon-tables==0.1.12
tensorflow==2.13.0
tensorflow-probability==0.21.0
```
