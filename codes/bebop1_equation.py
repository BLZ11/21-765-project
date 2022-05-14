"""Computes the BEBOPG1 energy functional and/or bond energy contributions
   BEBOP1 :: total BEBOPG1 atomization energy including zero-point vibrational energies (ZPVE) in kcal/mol
   BEBOPBondEnergies :: BEBOPG1 pi and sigma bond energy contributions and total atomization energy (including ZPE) all in kcal/mol
   resonance :: calculate the resonance energy (kcal/mol)
   strain :: calculate the strain energy (kcal/mol)"""

import numpy as np
import parameters as par

def BEBOP1(DistanceMatrix,AtomSym,MolOcc2s,MullikenPop):
    """Compute the total BEBOPG1 energy (SCF+ZPVE) at 0 K.
    
    Parameters
    ----------
    PairDis: :obj:'np.ndarray'
        Distance matrix in \AA
    AtomSymb: :obj:'np.ndarray'
        Array of atoms present in the molecular system
    MolOcc2s: :obj:'np.ndarray' or 'float'
        Molecular 2s occupation number   
    MullikenPop: :obj:'np.ndarray'
        Mulliken MBS population matrix condensed to atoms
    
    Returns
    -------
    BEBOP: :obj:'np.float64'
        Total BEBOPG1 atomization energy including ZPVE at 0K. 
    """
    
    BEBOP = 0
    SIZE = AtomSym.shape[0]
    for l in range(1,SIZE):
        for n in range(l):
            betaAB, zetaAB, R_sigma, D_AB = par.BEBOP_Pair(AtomSym[l],AtomSym[n]) # BEBOP atom-paired parameters
            BEBOP += -2 * betaAB * MullikenPop[l][n] # extended Hückel covalent bond correction
            BEBOP += D_AB * np.exp(zetaAB * (DistanceMatrix[l][n] - R_sigma)) # short-range repulsion correction
                
    # Calculate the hybridization energy
    n = 0
    for foo in AtomSym:
        if foo in np.array(['H','He']):
            continue # skip if 'H' and/or 'He' is present (no hybridization!)
        else:
            n2s,E2s2p = par.BEBOP_Atom(foo) # Hybridization parameters
            if type(MolOcc2s) != np.ndarray:
                BEBOP += (n2s - MolOcc2s) * E2s2p
            else:
                BEBOP += (n2s - MolOcc2s[n]) * E2s2p
                n += 1
    return BEBOP # in kcal/mol

def BEBOPBondEnergies(SigmaBondOdr, PiBondOdr, MullikenPop, DistanceMatrix, AtomSym, MolOcc2s, EROCBSQB3=None):
    """Calculate the individual bond energy contributions (and scaled bond energies with respect to CBS-QB3(0K) energy if optional)
    
       Parameters
       ----------
       SigmaBondOdr: obj:'np.ndarray' 
           Sigma bond order matrix 
       PiBondOdr: obj:'np.ndarray' 
           Pi bond order matrix 
       MullikenPop: obj:'np.ndarray'
           Mulliken MBS population matrix condensed to atoms 
       DistanceMatrix: obj:'np.ndarray'
           Distance matrix in \AA  
       AtomSym: obj:'np.ndarray'
           Array of atoms present in the molecular system  
       EROCBSQB3: obj:'float', optional
           Total atomization energy of EROCBSQB3 (optional)
       
       Returns
       -------
       Esig: obj:'np.ndarray'
           Gross sigma bond energies
       Epi: obj:'np.ndarray' 
           Gross pi bond energies
       Ecov: obj:'np.ndarray' 
           Gross covalent bond energies (i.e., Ecov = Esig + Epi)
       Enet_sig: obj:'np.ndarray'
           Net sigma bond energies 
       Enet_pi: obj:'np.ndarray'
           Net pi bond energies
       Enet: obj:'np.ndarray'
           Net bond energies (i.e., Enet = Enetsig + Enetpi)
       CompositeTable: obj:'np.ndarray'
           Matrix showing hybridization energy(diagonal elements),
           gross covalent bond energies(upper diagonal elements), 
           and net bond energies (lower diagonal elements)
       Enet_RN: obj:'np.ndarray'
           Renormalized net bond energies
       BEBOP: obj:'np.float64', 
           Total BEBOP atomization energy (with ZPVE) at 0 K
       """
    
    BEBOP = 0
    SIZE = AtomSym.shape[0]
    Esig = np.zeros((SIZE,SIZE)) # the gross sigma bond energy
    Epi = np.zeros((SIZE,SIZE)) # the gross pi bond energy
    Ecov = np.zeros((SIZE,SIZE)) # the gross covalent bond energy
    Rep = np.zeros((SIZE,SIZE)) # the short-range repulsion matrix
    
    # Calculate the BEBOP pair-wise energy
    for l in range(1,SIZE):
        for n in range(l):
                betaAB, zetaAB, R_sigma, D_AB = par.BEBOP_Pair(AtomSym[l],AtomSym[n]) # get parameters
                Esig[l][n] = -2 * betaAB * SigmaBondOdr[l][n] # calculate sigma gross bond energies
                Epi[l][n] = -2 * betaAB * PiBondOdr[l][n] # calculate pi gross sigma bond energies
                Ecov[l][n] = -2 * betaAB * MullikenPop[l][n] # calculate the covalent bond energies
                BEBOP += -2 * betaAB * MullikenPop[l][n] # added to the BEBOP energy
                Rep[l][n] = D_AB * np.exp(zetaAB * (DistanceMatrix[l][n] - R_sigma)) # short-range term
                BEBOP += Rep[l][n] 
                
                # Correct BEBOP covalent energies for short-range repulsion
                if np.abs(Ecov[l][n]) < 0.01:
                    Ecov[l][n] += Rep[l][n] # covalent bond energies
                else:
                    Esig[l][n] += Rep[l][n] * Esig[l][n] / Ecov[l][n] # sigma bond energy elements
                    Epi[l][n] += Rep[l][n] * Epi[l][n] / Ecov[l][n] # pi bond energy elements
                    Ecov[l][n] += Rep[l][n] # repulsive energy 
                    
  
    TBE = np.zeros(SIZE, dtype = np.float64) # total bond energy dummy array
    Enet = np.zeros((SIZE,SIZE), dtype = np.float64) # net bond energy 
    Ehybrid = np.zeros(SIZE, dtype = np.float64) # hybridization energy
    
    # Calculate the hybridization energy
    n = 0
    for foo in range(SIZE):
        if AtomSym[foo] in np.array(['H','He']):
            continue
        else:
            n2s,E2s2p = par.BEBOP_Atom(AtomSym[foo]) # Hybridization parameters
            if type(MolOcc2s) != np.ndarray:
                Ehybrid[foo] = (n2s - MolOcc2s) * E2s2p
                BEBOP += (n2s - MolOcc2s) * E2s2p
            else:
                Ehybrid[foo] = (n2s - MolOcc2s[n]) * E2s2p
                BEBOP += (n2s - MolOcc2s[n]) * E2s2p
                n += 1
                
    # Calculate total gross bond energy
    for i in range(1,SIZE):
        for j in range(i):
            TBE[i] += Ecov[i][j]
            TBE[j] += Ecov[i][j]
    
    Enet_sig = np.zeros((SIZE,SIZE), dtype = np.float64) # calculate sigma net bond energies
    Enet_pi = np.zeros((SIZE,SIZE), dtype = np.float64) # calculate pi net bond energies
    for i in range(1,SIZE):
        for j in range(i):
            # scale the bond energies
            Factor = 1.0 + Ehybrid[i] / TBE[i] + Ehybrid[j] / TBE[j] 
            Enet[i][j] = Ecov[i][j] * Factor 
            Enet_sig[i][j] = Esig[i][j] * Factor 
            Enet_pi[i][j] = Epi[i][j] * Factor 
    CompositeTable = np.zeros((SIZE,SIZE), dtype = np.float64)     # Composite Table: Eii (hybrid),Eji (gross),
                                                                   #                  Eij (net)   ,Ejj (hybrid) 
    for i in range(SIZE):
        for j in range(l):
            CompositeTable[j][i] = Ecov[i][j] # upper diagonal elements are filled with Ecov elements
            CompositeTable[i][j] = Enet[i][j] # lower diagonal elements are filled with Enet elements
    np.fill_diagonal(CompositeTable, Ehybrid) # Fill the diagonal elements of the CompositeTable with the Ehybrid contributions
    if EROCBSQB3 == None:
        return (Esig, Epi, Ecov, Enet_sig, Enet_pi, Enet, CompositeTable, BEBOP)
    else:
        Enet_RN = Enet * EROCBSQB3 / BEBOP
        return (Esig, Epi, Ecov, Enet_sig, Enet_pi, Enet, CompositeTable, Enet_RN, BEBOP)
    
def resonance(Data):
    """Compute the resonance energy of the aromatic system (in kcal/mol)
    
        Parameters
        ----------
        Data: obj:'dict'
            A dictionary containing the following inside the tuple: (total number of bond orders, 
            total number of unique bonds, reference bond orders, atoms arrays) for each reference
            bond
        
        Return
        ------
        res_E: obj:'np.float'
            Return the total resonance energy change
    """
    
    res_E = 0
    for l in Data.keys():
            BO_tot, n, ref_BO, atoms = Data[l] 
            betaAB = par.BEBOP_Pair(atoms[0],atoms[1], res_strain = True) 
            res_E += (BO_tot - np.int(n) * ref_BO) * betaAB  
    return res_E
    
def strain(Data):
    """Compute the ring strain energy for a ring (in kcal/mol)
    
    Parameters
    ----------
    Data: obj:'dict'
        A dictionary containing the following inside the tuple: (total number of bond orders, 
        total number of unique bonds, reference bond orders, atoms arrays) for each reference
        bond
        
    Return
    ------
    res_E: obj:'np.float'
        Return the total resonance energy change
    """
    
    strain_E = 0
    betaAB = par.BEBOP_Pair('C','C', res_strain = True)
    for l in Data.keys():
        BO_tot, n, ref_BO = Data[l] 
        strain_E += (np.int(n) * ref_BO - BO_tot) * betaAB 
    return strain_E