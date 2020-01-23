"""
Programma mirato a trovare similitudini tra gruppi di sequenze di RNA
Semplice idea: date N sequenze di eguale lunghezza, costruisco una "matrice delle distanze"
ossia vedo il grado di similitudine delle stringhe dando valore al numero di residui (basi
azotate in comune)

es GUUACAC e CAAUGUG di lunghezza 7 hanno similitudine 
0 se si conta la posizione del residuo
6 se si conta in assoluto la presenza di una base

"""

#%%
import numpy as np

def String_Similarity(s1, s2):

    if (len(s1) != len(s2)):
        raise ValueError ("Strings must be of equal size")
    
    simil    =   np.zeros(len(s1))  

    for i in range(len(s1)):
        if (s1[i]==s2[i]):
            simil[i] =  1 
        else :simil[i] = 0
    
    return np.sum(simil)/len(s1)

#%%
#GRUPPO LUNGO

N = 25

A   =   "GCUGGGGUGGGGCGGAUCGGUGUUG"
B   =   "AGACGAGGCCGGGCUUGUCCCCGGC"
C   =   "GGAGCUGGCCCUGCGGGGCCCGGCG"

Gen_1 = (A, B, C)

# costruisco matrice somiglianze
dist   =  np.zeros((3,3))
i = 0


for a in Gen_1:
    j=0
    for b in Gen_1:
        dist[i,j]   =   String_Similarity(a,b)
        j+=1
    i+=1

# sommo su un indice e ottengo una stima della somiglianza

sameness = np.sum(dist, axis=0)

print("La matrice delle somiglianze della Generazione1  è \n")
print(dist)
print("Il centroide è dato dalla sequenza %s" %str((np.argmax(sameness)+1)))
print("Con il valore di %s" %str(np.max(sameness)))

# %%

D   =   "CGGUGUUGCU"
E   =   "CGCUGUGGUC"
F   =   "GUGGUCCCCG"
G   =   "GGGGUGGGGC"
H   =   "CGAGGCCGGG"
F   =   "GCGGGGCCCG"
Dneg=   "AGCAACACCG"

Gen_2   = (D, E, F, G, H, F, Dneg)

# costruisco matrice somiglianze
dist   =  np.zeros((7,7))
i = 0

for a in Gen_2:
    j=0
    for b in Gen_2:
        dist[i,j]   =   String_Similarity(a,b)
        j+=1
    i+=1

# sommo su un indice e ottengo una stima della somiglianza

sameness = np.sum(dist, axis=0)

print("La matrice delle somiglianze della Generazione 2 è \n :")
print(dist)
print("Il centroide è dato dalla sequenza %s" %str((np.argmax(sameness)+1)))
print("Con il valore di %s" %str(np.max(sameness)))

# %%


# %%
