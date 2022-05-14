# Bond Energy/Bond Order Population (Beta Version)

The bond energy/bond order population ([BEBOP]()) program is a computational chemistry algorithm that computes accurate molecular energies at equilibrium and bond energies using well-conditioned Hartree-Fock orbital populations and bond orders.

This is a beta version of the code. 

## Codes

All Python codes are available under the ```codes``` directory.  
## Running the code
1. Optimize your molecular structure using Gaussian. Use the following recommended keywords:

   ```# Opt B3LYP/CBSB7```

2. Run Hartree-Fock on the optimized structure in Gaussian. Use the following keywords:

   ```# SP ROHF/CBSB3 Pop=(Full) IOp(6/27=122)```

3. Execute ```bebop.py```. Execute this in the command line:
  ```bash
./bebop.py -f {name_file} --be --sort --json > {name_file}.bop
```
where ```{name_file}``` is the Hartree-Fock Gaussian output file. 

Some examples of Hartree-Fock output files and BEBOP results are found in the ```tests``` and ```output``` directories, respectively. 

## Usage
Some details of the argument parser used in ```bebop.py``` source code.

```bash
usage: bebop.py [-h] -f F [--be] [--sort] [--json]
compute BEBOP atomization energies and bond energies (i.e., gross and net)
optional arguments:
  -h, --help  show this help message and exit
  -f F        name of the Gaussian Hartree-Fock output file
  --be        compute BEBOP bond energies (net and gross bond energies)
  --sort      sort the net BEBOP bond energies (from lowest to highest in
              energy)
  --json      save the job output into JSON
```
## Directory Information
* ```codes```- contains the BEBOP program codes
* ```tests```- contains some Gaussian Hartree-Fock output files
* ```output```- contains BEBOP output file examples

## Important Documents
For information regarding the modules present in the BEBOP program, visit:
* ```codes/code_explanation.md```

For infromation regarding the data printed by the BEBOP program, visit:
* ```output/out_explanation.md```



## License
[MIT](https://choosealicense.com/licenses/mit/)
