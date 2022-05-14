# BEBOP Program

Below are brief general dcoumentation on the BEBOP codes used by BEBOP.  In general, the BEBOP program **consists** of:
* ```bebop.py```
* ```bebop1_equation.py```
* ```bebop_g1.py```
* ```parameters.py```
* ```read_output.py```
* ```roothan.py```
* ```spatial_geom.py```

## ```bebop.py```

The main executable file for the BEBOP program. This program parses the arguments the user inputs in the command line to retrive the BEBOP data. Moreover, it prints in a human readable file that prints the BEBOP total energies and bond energies.   

## ```bebop1_equation.py```

This API computes contains prototypes that computes any type of bond energies (i.e., net and gross) and total energies from the BEBOP equation.

## ```bebop_g1.py```

The core module of the BEBOP program. This API retrives data from:
*  the Hartree-Fock output via ```read_output.py``` ,
* the calculated sigma and pi bond orders via ```roothan.py``` ,
* and the projection angles and distance matrix via ```spatial_geom.py```
to compute BEBOP energies via  ```bebop1_equation.py```'s prototypes. 

Users can import this module in a Jupyter Notebook to generate data present in ```bebop1_equation.py```. Documentation in the near future will be present for users who desire to use Jupyter Notebook. 

## ```parameters.py```

Module contains all parameters for the first row elements that needs to be retrieve to compute BEBOP energies via ```bebop1_equation.py``` . These parameters were obtained from non-linear regression for 36 reference species by:
* Fitting BEBOP potential energy data to accurate energies 
* Constraining the classical turning point and the equilibrium point

More information on this will be available soon.  

## ```read_output.py```

Module retrieves the necessary data from the Hartree-Fock output file to compute BEBOP energies (i.e., atoms, xyz coordinates, etc.)

## ```roothan.py```

This module computes sigma and pi bond orders by feeding the distance matrix, projection angles, atomic populations, and the total  bond order. It then computes the bond orders using two equations found in this [paper](https://doi.org/10.1063/1.1748100).


## ```spatial_geom.py```

This module uses the xyz coordinates retrieved from ```read_output.py``` to compute the distance matrix, a lower diagonal matrix consisting of distances between pairs of atoms, and the projection angles. 
