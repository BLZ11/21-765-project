# BEBOP Output File
We will give an overview of the data printed by the BEBOP program.

## ```HARTREE-FOCK OUTPUT:```
Contains the path where the Hartree-Fock output file is located. 

## ```BEBOP ATOMIZATION ENERGY (0 K) = ```
Total BEBOP atomization energy computed with zero-point vibration energy. 
This energy (E_0) computed as:

E_0 = E_cov + E_rep + E_hyb

where E_cov is the energy of a covalent bond (from approximate quantum chemistry),
E_rep is the nuclear-nuclear repulsive interaction that occurs at small distances,
and E_hyb is the orbital hybridization energy from exciting an electron from the 2s
to 2p orbital.

For information on BEBOP theory please read the [BEBOP manuscript](https://chemrxiv.org/engage/chemrxiv/article-details/624dff79855ee54b39e40518)
