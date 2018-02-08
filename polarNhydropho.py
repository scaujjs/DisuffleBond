import numpy as np

def normalizeHydro(x):
    y=(x+-0.09649999999999999)/(1.85+2.02)
    return y

def twoMetrix(residue):
    metrix=np.zeros((1),np.float32)
    dict_Hydro={'A':-0.17481,'R':-0.23424,'N':-0.14122,'D':0.26705,'C':-0.05853,'Q':0.11977,
                  'E':0.03450,'G':0.45310,'H':-0.13346,'I':-0.06111,'L':-0.28075,'K':0.05517,
                  'M':0.03708,'F':-0.06886,'P':0.21796,'S':-0.06886,'T':-0.02752,'W':-0.04302,
                  'Y':-0.54690,'V':-0.34277}
    # Normalisation based on solvent accessibility data for
    # graphics-built Ala-X-Ala tri-peptide models given in Table
    # II of Samanta et al. (2002) obtained using NACCESSsoftware. C a
    #  considered as main-chain atom.
    metrix[0]=dict_Hydro[residue]
    return (metrix)

#x=0
#for key in dict_residue:
#    a=normalizeHydro(dict_residue[key])
#    print(a)
