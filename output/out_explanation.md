# BEBOP Output File
This is an overview of the data printed by the BEBOP program.

## ```HARTREE-FOCK OUTPUT:```
Contains the path where the Hartree-Fock output file is located. 

## ```BEBOP ATOMIZATION ENERGY (0 K) = ```
Total BEBOP atomization energy computed with zero-point vibration energy. 
This energy (E_0) is computed as:

E_0 = E_cov + E_rep + E_hyb

where E_cov is the energy of a covalent bond (from approximate quantum chemistry),
E_rep is the nuclear-nuclear repulsive interaction that occurs at small distances 
where bonding is formed, and E_hyb is the orbital hybridization energy from exciting 
an electron from the 2s to 2p orbital.

For more information on BEBOP theory, please read the [BEBOP manuscript](https://chemrxiv.org/engage/chemrxiv/article-details/624dff79855ee54b39e40518).

## ```GROSS SIGMA```-, ```GROSS PI```-, and ```GROSS TOTAL```- ```BOND ENERGIES INCLUDING REPULSION``` Data

These three tables depicts the gross bond energies. The tables have similar characteristics which include:
1. the positive integer numbers presented in the column and rows represents the center number (or the number representing the position of the 
atom relative to the top atom introduced in the matrix)
2. lower diagonal square gross bond energy matrices all in kcal/mol (the dimensions are equal to the amount of atoms present in the molecule)
3. the sum of all elements within the array gives the total gross energy value (depicted as ```TOTAL GROSS...```)

The total bond energies are computed similar to the previous equation, except no E_hyb. Sigma bond energies and pi bond energies were computed
using that equation, but the bond orders were computed from the BEBOP program under the ```roothan.py``` module. The bond orders used for the total 
gross bond energy, however, were from the Hartree-Fock output (or the sum of the sigma and pi bond orders). 

## ```NET SIGMA```-, ```NET PI```-, and ```NET TOTAL```- ```BOND ENERGIES INCLUDING HYBRIDIZATION``` Data

These tables are structured exactly like the gross bond tables. The difference here is that the bond energies include hybridization.
The net total energy should equal to E_0.

## ```COMPOSITE TABLE```
This table are similar to the previous tables, but:
* Net total bond energies quantities from the previous table are the lower diagonal elements,
* the diagonal elements are the hybridization energy,
* and gross total bond energies quantities from the previous table are the upper diagonal elements

At the bottom of this table, the ```TOTAL HYBRIDIZATION ENERGY``` value is introduced which is equal to 
the sum of the diagonal elements of the composite table (i.e., matrix trace). 

## ```SORTED NET BONDING ENERGIES (LOWEST TO HIGHEST)```
Tabulated values for sorted net bond energies excluding anti-bonding interactions. 
The sorted bond identities (with the center number) are shown on the left and the value
is present on the right.
