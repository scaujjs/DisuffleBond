import numpy as np

def normalizepolar(x):
    y=(x-0.16684)/(1.36+0.84)
    return y

def normalizerHydro(x):
    y=(x+0.0965)/(1.85+2.02)
    return y


def twoMetrix(residue):
    metrix=np.zeros((2),np.float32)
    dict_Hydronorm={'T':-0.01124, 'V': 0.00685, 'R': -0.18437, 'C': 0.08695, 'P': -0.09134, 'F': 0.31693,
                    'D': -0.29289, 'Q': -0.12494, 'W': 0.50297, 'K': -0.23088, 'Y': 0.26783,'M':0.08437,
                    'E':-0.49703,'A':-0.01899,'I':0.10504,'S':-0.00866,'N':-0.08359,'H':-0.01899,'L':0.16964,
                    'G':0.02235,'X':0}
    dict_polarnorm = {'Q':-0.40766,'R':-0.45766,'H':0.14689,'M':0.50144,'A':-0.10311,'F':0.50144,
                  'I':0.51962,'V':0.41962,'E':-0.42584,'L':0.47416,'C':0.54235,'N':-0.29402,
                  'Y':0.07416,'D':-0.43947,'W':0.32416,'G':-0.26220,'K':-0.61220,'T':-0.19856,
                  'S':-0.30311,'P':0,'X':0}
    # Normalisation based on solvent accessibility data for
    # graphics-built Ala-X-Ala tri-peptide models given in Table
    # II of Samanta et al. (2002) obtained using NACCESSsoftware. C a
    #  considered as main-chain atom.
    metrix[0]=dict_Hydronorm[residue]
    metrix[1]=dict_polarnorm[residue]
    return (metrix)


dict_Hydro={'A':-0.17,'R':-0.81,'N':-0.42,'D':-1.23,'C':0.24,'Q':-0.58,'E':-2.02,'G':-0.01,'H':-0.17,
            'I':0.31,'L':0.56,'K':-0.99,'M':0.23,'F':1.13,'P':-0.45,'S':-0.13,'T':-0.14,'W':1.85,
            'Y':0.94,'V':-0.07}

dict_polar={'L':1.21,'I':1.31,'V':1.09,'F':1.27,'M':1.27,'W':0.88,'A':-0.06,'C':1.36,'G':-0.41,
            'Y':0.33,'T':-0.27,'S':-0.50,'H':0.49,'Q':-0.73,'K':-1.18,'N':-0.48,'E':-0.77,'D':-0.80,
            'R':-0.84}


if(0):
    x=0
    for key in dict_polar:
       x+=dict_polar[key]
       print(x)
    # x=3.17

    aver=3.17/19
    print(aver)
    #aver==0.1668421052631579~ 0.16684

if(0):
    for key in dict_polar:
        a=normalizepolar(dict_polar[key])
        print(str(key)+':'+str(a))

if(0):
    x = 0
    for key in dict_Hydro:
        x += dict_Hydro[key]
    print(x)  # x=-1.93

    aver=-1.93/20
    print(aver)  # aver== -0.0965

if(0):
    numOfKey=0
    for key in dict_Hydro:
        a=normalizerHydro(dict_Hydro[key])
        numOfKey+=1
        print(str(key)+':'+str(a))
    print(numOfKey)