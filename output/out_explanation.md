# BEBOP Output File
We will give an overview of the data printed by the BEBOP program.

## ```HARTREE-FOCK OUTPUT:```
Contains the path where the Hartree-Fock output file is located. 

## ```BEBOP ATOMIZATION ENERGY (0 K) = ```
Total BEBOP atomization energy computed with zero-point vibration energy. 
This energy ($E_{\text{AT},0}$) computed as:

$$\Delta E_{\text{AT},0} = E_{\text{cov}} + E_{\text{rep}} + E_{\text{hyb}}$$

where $E_{\text{cov}}$ is the energy of a covalent bond (from approximate quantum chemistry),
$E_{\text{rep}}$ is the nuclear-nuclear repulsive interaction that occurs at small distances,
and $E_{\text{hyb}}$ is the orbital hybridization energy from exciting an electron from the $2s$
to $2p$. 
