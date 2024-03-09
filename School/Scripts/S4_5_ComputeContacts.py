#!/usr/bin/env python3

import random, math, numpy, string, sys, os # modules de base
import argparse # gestion ds args
import S4_S5_structureToolsM1BIBS as structureToolsM1BIBS # notre module maison
import matplotlib.pyplot as plt  # pour le plot final

########################################################################################################
#
#                   FUNCTIONS
#
########################################################################################################



def computeContactMatrix(d_coords, chain, mode, writematrix = False) :
    """ input : list of lists which contains the coords of every atoms
    output : distance matrix
    """
    nbres = len(d_coords[chain]["reslist"])
    distmat = numpy.zeros((nbres,nbres))

    if writematrix:
        fout = open("mat.txt", "w")

    for i in range(nbres) :
        j = i + 1
        d_coordi = d_coords[chain][d_coords[chain]["reslist"][i]]
        while j < nbres :
            d_coordj = d_coords[chain][d_coords[chain]["reslist"][j]]
            dij = structureToolsM1BIBS.compDistance(d_coordi, d_coordj, mode)
            distmat[i,j], distmat[j,i] = dij, dij
            j += 1
            if writematrix:
                fout.write("%s\n"%(dij))

    if writematrix :
        fout.close()

    return distmat            


def getResid(contactsList, dPDB) :
    """input: contact pairs extracted from the matlist (indexes from 0 to Nb of residues -1)
       purpose:   indexes do not necessarily correspond to the PDB numbers, this function
       based on the indexes extracted from the matrix, search for the corresponding residues 
       in the dict
       output: contact pairs with PDB numbers"""

    resPairs = []
    
    for pair in contactsList :
        resPairs.append([dPDB[chain]["reslist"][pair[0]], dPDB[chain]["reslist"][pair[1]]])        

    return resPairs



#####################################################################################
#
#                       MAIN
#
#####################################################################################




# Get Arguments
#===============

parser = argparse.ArgumentParser( prog = 'ComputeContact.py',
                    description = 'Calculates the contact map of a given chain in a PDB file',
                    epilog = """if you use this program please cite abcd et al, 2023 """)

# add args
parser.add_argument("--pdb", help="pdb to treat", metavar = "input", type=str, required = True)
parser.add_argument("--chain", help="chain to treat", metavar = "chain", type=str, default = "first")
parser.add_argument("--mode", help="""if mode = atom, compute the distance between all the atoms of the two 
                 residues and return the smallest distance.
                 If  mode = 'centroid', compute the distance between the two centers of mass
                 of the two residues and return it.""", metavar = "comp_mode", type=str, default = "atom")
parser.add_argument("--oplot", help="name of the contact map output file", metavar = "output_map",type=str, default = "distance_matrix.eps")
parser.add_argument("--opairs", help="name of the text file containing all distance pairs", metavar = "output_dist", type=str, default = "contactsPairs.txt")
parser.add_argument("--seuil", help="threshold to define a contact (in Angstrom)", metavar = "5.0",type=float, default = 5.0)
parser.add_argument("--writemat", help="write the matrix of pairwise distances", action='store_true')

# parse args and store the corresponding values in var
args = parser.parse_args()

infile = args.pdb
chain = args.chain
mode = args.mode
outplot = args.oplot
outname = args.opairs
seuil = args.seuil
writemat = args.writemat


# Computes distances between every pair of residues, stores them in distmat (array) and returns a contact matrix in eps format
#=======================================================================================

# parses the pdb file
dPDB = structureToolsM1BIBS.PDB_parser(infile)

if chain == "first" : # by default, computes the map for the first chain
    chain = dPDB["chains"][0]

# computes the distances
print("computing the distances according to the", mode)
matdist = computeContactMatrix(dPDB, chain, mode, writematrix = writemat)

# plots
fig = plt.figure()# use param figsize=(8,8) to specify a given size
plt.pcolormesh(matdist)
plt.title('map distance of {} chain {}'.format(infile, chain), fontweight ="bold")
plt.colorbar()
#plt.show() --> if you do not want to show the plot
plt.savefig(outplot,dpi=100)  # if you want to save the plot

# extracts the contacts
contacts = structureToolsM1BIBS.extractContactResidues(matdist, seuil)
resPairs = getResid(contacts, dPDB)

# writes the pairs in contact
fout = open(outname, "w")

for pairs in resPairs :
    fout.write("%s\t%s\n"%(pairs[0],pairs[1]))

fout.close()



