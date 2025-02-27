# QuVax
### mRNA design guided by folding potential

QuVax is a python command line application designed to sample various mRNA sequences and predict their structures in a bilevel optimization formalism for use in hybrid classical and quantum computing workflows. The tool is inspired by the rise of mRNA vaccines and aspires to be a high-throughput tool to aid in their design.

## How It Works

A population of mRNA sequences are generated randomly from an input protein sequence in FASTA file format using the genetic code. In each iteration of the optimization process, changes are made to the sequences of the population and each iteration of changes is referred to as a generation. The folding energies and secondary structures are then determined according to a Hamiltonian [1] for each sequence in the population. This process is repeated for a user-defined number of iterations. Here is a basic example of the process in pseudocode:

```
for codon in codon_iterations:
  for sequence in population:
    Change mRNA or codon sequence
    for fold in folding_iterations:
      Search for minimum of fold space for given codon sequence
    if energy < minimum_energy:
      Set new minimum energy and record sequence and secondary structure
```

There are multiple optimizers that can be used to generate changes in the mRNA sequences or explore the fold space of a single mRNA sequences implemented in this tool. More details about the tool are outlined below.

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

QuVax has three functions: bilevel codon sequence and folding optimization, folding optimization for a single RNA sequence, and a small analysis suite. Bilevel optimization is executed using ```design.py```, folding with ```fold.py```, and analyses with ```analyze.py```.

### Optimization of mRNA sequence and fold for a given protein sequence with ```design.py```

Example execution of ```design.py```:

```
$ python design.py -i examples/spike_trim.fasta
```

You can see command line options and their defaults in the log file or through the help option: ```$ python design.py -h```

### Database options

QuVax has options to integrate with SQLite and Postgres databases. SQLite is the default option and the output option (```-o```) is just the name of the file containing the database. Connecting to a Postgres database requires some more information. ```-o``` should now list the database you intend to write to in your Postgres instance. An ```.ini``` file, supplied to the ```-in``` option, is required in the following format:

```
[postgresql]
user=(user)
password=(password)
```

where parentheticals, ```()```, are replaced by your values. Any other fields included in this file will be ignored! The user should have permissions to create databases in the Postgres instance, although you could create one through the client of your choice and connect to it via QuVax.

### Storing Multiple Optimizations in One Database

If you were to run the above execution of ```design.py``` twice, data from both executions will be stored within the same database (default output is SQLite database with name ```quvax.db```). Each optimization will be associated with a hash value that can be used to identify the sequence. The hash value will be used to resume optimizations and analyze optimizations.

### Continuing an Optimization

An option to continue from the end of a previous execution of ```design.py``` is available. To do this, the input needs to be a valid SQLite database and the ```--resume``` option must be specified.

```
$ python design.py -i quvax.db --resume
```

Most options for ```design.py``` are read from the database. The input database file will be used as the output file. You can use ```$ python design.py --resume -h``` to see the available options a user can specify when resuming an optimization.

### Analyzing Optimizations

The data in the database can be viewed with any database editor (Heidi, DBeaver, etc.). Custom python or SQL scripts can be written to analyze an optimization process since all of the information is stored in a standard SQLite database format. However, QuVax includes some analysis modules in ```analyze.py```. Example execution of ```analyze.py```:

```
$ python analyze.py command [OPTIONS]
$ python analyze.py fe_trajectory -i quvax.db
```

For analysis options: ```$ python analyze.py -h```
For analysis command options: ```$ python analyze.py [command] -h```

### Folding Individual Codon Sequences

If you have a codon sequence in mind for folding optimization alone, you can use ```fold.py``` to do this.

```python fold.py -i GGGAAACUGGAAGGCGGGGCGAGCUGCAGCCCCAGUGAAUCAAAUGCAGC```

Currently, the simulated annealer (```-s SA```), the Monte Carlo (```-s MC```), and an Exact Solver (```-s ES```) are available for RNA folding and both can be used here. The output of ```fold.py``` is a dot-bracket file (default: quvax.dot) which has the following format

```
\> Folded energy: -81.0
GGGAAACUGGAAGGCGGGGCGAGCUGCAGCCCCAGUGAAUCAAAUGCAGC
.....(((((............[[[[[[...)))))........]]]]]]
```

Where the final line is the folded secondary structure in dot-bracket notation. For more folding options: ```$ python fold.py -h```

### Analyzing RNA structures

After a structure is generated, you can analyze that structure in a variety of ways (see ```python analyze.py -h``` for list). Here's an example that will compare two structures based on their connectivity tables:

```python analyze.py compare_ct -i structure.ct -r ref_structure.ct```

where you would replace the filenames for ```-i``` and ```-r``` with valid filenames. Note: Most analyses only currently work with connectivity tables.

## Testing

To run tests:

```$ coverage run -m pytest``` or ```$ pytest```

Test python files are in ```tests/test_cases```. For running a specific set of tests (no path required):

```pytest -k "test file name without .py"```

## List of Dependencies in requirements.md

```
argparse==1.4.0
biopython==1.81
dwave-samplers==1.2.0
pandas==2.0.3
psycopg2-binary==2.9.9
pytest-cov==4.1.0
python-codon-tables==0.1.12
tensorflow==2.16.1
tensorflow-probability==0.24.0
tf_keras==2.16.0
```

## Known Issues

Installing tensorflow on Macs can have machine-dependent issues (x86 vs ARM, etc). Instead of pip, you may want to try installing via conda. Or, if you have an Intel (x86) Mac, try installing tensorflow 2.12. As a last resort, you can download the whl for the desired version from PyPi.

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
