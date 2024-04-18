# QuVax
### mRNA design guided by folding potential

## Getting Started

Clone this repository. Install your preferred version of conda (Anaconda or Miniconda). Create Conda environment:

```
$ conda create -n quvax python=3.11
$ conda activate quvax
```

Install dependencies:

```
$ python -m pip install -r requirements.md

OR

$ conda install --file requirements.md -c conda-forge
```

*The latter currently does not work.

## Using QuVax

QuVax has two primary functions: codon sequence and folding optimization and a small analysis suite. Optimization is executed with ```design.py``` and analyses with ```analyze.py```.

### Optimization of mRNA sequence for a given protein sequence with ```design.py```

QuVax treats this problem as a bilevel optimization problem or a nested optimization problem. A population of mRNA sequences are generated randomly from an input protein sequence. The folding energies are then determined according to a Hamiltonian [1]. Changes are then proposed to the mRNA sequences in the population and the folding energies are recalculated. This process is repeated for a user-defined number of iterations. Example execution of ```design.py```:

```
$ python design.py -i examples/spike_trim.fasta
```

### Continuing an Optimization

An option to continue from the end of a previous execution of ```design.py``` is available. To do this, the input needs to be a valid SQLite database and the ```-resume``` option must be specified. You can still specify other options for ```design.py```, but they will be overwritten by option values read by the input database. The input database file will be used as the output file.

```
$ python design.py -i quvax.db --resume
```

### Analyzing Optimizations

```analyze.py``` then takes in the output of ```design.py```, which is simply a SQLite3 database with default extension ```.db``` (this can be changed on command line if desired). The data in the database can be viewed with any database editor (Heidi, DBeaver, etc.) and can be read by the analysis modules of QuVax. Example execution of ```analyze.py```:

```
$ python analyze.py -i quvax.db -at trajectory
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

## Known Issues

[Tensorflow does not support Mac](https://github.com/tensorflow/tensorflow/issues/61382) so you will have to use ```tensorflow==2.11``` and ```tensorflow-probability==0.19``` in ```requirements.md``` if installing via pip. All tests currently pass with these versions of the packages. Alternatively, you can try to install these packages manually via wheel from the PyPi repositories to have access to the newer versions.

## References

[1] Fox, D. M., MacDermaid, C. M., Schreij, A. M. A., Zwierzyna, M., & Walker, R. C. (04 2022). RNA folding using quantum computers. PLOS Computational Biology, 18(4), 1–17. doi:10.1371/journal.pcbi.1010032

```
@article{10.1371/journal.pcbi.1010032,
    doi = {10.1371/journal.pcbi.1010032},
    author = {Fox, Dillion M. AND MacDermaid, Christopher M. AND Schreij, Andrea M. A. AND Zwierzyna, Magdalena AND Walker, Ross C.},
    journal = {PLOS Computational Biology},
    publisher = {Public Library of Science},
    title = {RNA folding using quantum computers},
    year = {2022},
    month = {04},
    volume = {18},
    url = {https://doi.org/10.1371/journal.pcbi.1010032},
    pages = {1-17},
}
```
