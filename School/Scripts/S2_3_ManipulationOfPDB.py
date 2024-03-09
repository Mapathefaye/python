#!/usr/bin/env python3

import string, sys



import S2_S3_structureToolsM1BIBS as structureToolsM1BIBS


# gestion des arguments
#=======================

#infile = sys.argv[1]

infile = "./School/1brs.pdb"


#=======================
#      MAIN
#=======================


# creation du dico dPDB contenant les coords 3D stockes dans infile
"""on appelle la fonction parser_PDB a partir du module PDBTools
que nous avons cree"""
dPDB = structureToolsM1BIBS.PDB_parser(infile)



#####  1ERE PARTIE #####
#----------------------#


# get number of chains
#-----------------------
print ("Nb of chains:",len(dPDB["chains"]))


# get nb res per chain + total nb of residues
#---------------------------------------------
total = 0
for chaini in dPDB["chains"] :
    print ("Chain %s has %d residues:"%(chaini, len(dPDB[chaini]["reslist"])))
    total+=len(dPDB[chaini]["reslist"])
    # eq a total = total + len(dPDB[chaini]["reslist"])
    # de facon generale : i+=x <---> i = i + x

print("Total nb of residues", total)



# get nb LYS/GLN
#-----------------
totLYS = 0 # initiation de nb LYS (total)
totGLN = 0 # idem pr GLN

for chaini in dPDB["chains"] :
    cmtLYS = 0 # initiation de nb LYS (per chain)
    cmtGLN = 0 # idem pr GLN

    for resi in dPDB[chaini]["reslist"] :
        # compte nb LYS/GLN par chaine
        if dPDB[chaini][resi]["resname"] == "LYS" :
            cmtLYS+=1
        elif dPDB[chaini][resi]["resname"] == "GLN" :
            cmtGLN+=1
            
    # affiche nb LYS/GLN pour la chaine i
    print("Chains {} has {} Lysines\n".format(chaini, cmtLYS))
    print("Chains {} has {} Glutamines\n".format(chaini, cmtGLN))
 
    
    # compte le nb total de GLN/LYS sur tte la prot
    totLYS+= cmtLYS
    totGLN += cmtGLN

# affiche totaux            
print("There are {} LYS".format(totLYS))
print("There are {} GLN".format(totGLN))



# get nb aa per chain
#---------------------
"""la fct est stockee dans structureToolsM1BIBS.py car cette fct est assez generique et peut etre
   amenee a etre appelee ds d'autres programmes donc on veut y avoir acces facilement"""

# ex avec nb de PHE

d_nbPHE = structureToolsM1BIBS.getNbAA(dPDB, "PHE")

for chaini in dPDB["chains"] :
    print("Nb of PHE {} in chain {}".format(d_nbPHE[chaini], chaini))



#####  2NDE PARTIE #####
#----------------------#

# les fonctions sont toutes dans structureToolsM1BIBS.py. On les reutilisera probablement.

# il y a bien sur plein d'autres possibilites, ceci est l'une des solutions simples.

# ex de calcul de distance entre res 12 chaine A et res 8 chaine A

d_res8 = dPDB["A"]["8"]
d_res12 = dPDB["A"]["12"]
print(d_res12)
dist_8_12_minval = structureToolsM1BIBS.compDistance(d_res8, d_res12, mode = "atom")
print("dist between res 8 and res 12 (minval): ",dist_8_12_minval)

dist_8_12_centroids = structureToolsM1BIBS.compDistance(d_res8, d_res12, mode = "centroid")
print("dist between res 8 and res 12 (cent): ", dist_8_12_centroids)


####  !!!!! PROFITEZ DE PYMOL POUR VERIFIER VOTRE PROGRAMME !!!! #####
