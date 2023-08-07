# QuVax
### mRNA design guided by folding potential

## Getting Started

Create Conda environment

### Mac

```bash conda env create --file requirements.txt```

### Linux

Create and activate environment

```$ conda env create -n quvax python=3.11
   $ conda activate quvax
```

Install dependencies

```$ python -m pip install -f requirements.txt

OR

$ conda install --file requirements.txt -c conda-forge```

*The latter currently does not work.

Example execution

```$ python design.py```

### List of Dependencies in requirements.txt
```
biopython==1.81
dwave-neal==0.6.0
pandas==2.0.3
pip==23.2.1
pipdeptree==2.12.0
python-codon-tables==0.1.12
tensorflow==2.13.0
tensorflow-probability==0.21.0
```
