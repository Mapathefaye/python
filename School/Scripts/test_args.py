import sys, argparse


print(sys.argv)



parser = argparse.ArgumentParser()

parser.add_argument("-f", "--file", help="nom du fichier a lire", type=str, required = True)
parser.add_argument("-m", "--mode", help="mode de calcul de la distance inter-residus", type=str, required = False, default = "atom")

#parser.add_argument("-t", "--type", help="type de molecule", type=str, default = "prot", nargs = "?", const="nucl")
parser.add_argument("-w", "--weigth", action="store_true", help="pondere les termesde mafct")

parser.add_argument("-p", "--peigth", help="ponderation pou non, si pondere, val par defaut = 0.5", type=float, default = False, nargs = "?", const=0.5)

args = parser.parse_args()

infile = args.file
weigth_fun = args.weigth
mode_calcul = args.mode
print(infile)
print("le calcul se fera sur :", mode_calcul)

if weigth_fun :
    print("ok je pondere !")
else:
    print("pas de ponderation pour l'operation")


peigth_fun = args.peigth

print(peigth_fun)
