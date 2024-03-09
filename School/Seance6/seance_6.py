import os, glob

# 1) créer un répertoire "Seance6" (je suis déjà au bon endroit)
#os.mkdir("Seance6")
# 2) téléchargement et extraction des fichiers fasta (fait manuellement)
# 3) Créer un rép start dans Senace6
# je me place d'abord dans Seance6
#Seance6 = "C:\\Users\\Etudiant FST\\OneDrive - UPEC\\Bureau\\Master1 BIBS\\S2BIBS\\Python\\School\\Seance6"
os.chdir("Seance6")
##os.mkdir("Start")
# 4) 5)
# Affiche le chemin de start
#Start = "C:\\Users\\Etudiant FST\\OneDrive - UPEC\\Bureau\\Master1 BIBS\\S2BIBS\\Python\\School\\Seance6\\Start"
os.chdir("/School/Start") # me place dans le rep
os.getcwd() # affiche le chemin où je suis
os.listdir() # liste rep et file du rep
#files = "C:\\Users\\Etudiant FST\\OneDrive - UPEC\\Bureau\\Master1 BIBS\\S2BIBS\\Python\\School\\Seance6\\Start\\files"
os.chdir("files") # me place dans le rep files
os.listdir()   # list le contenu de files
# 6) identifier les items correspondant à des fichiers, à des répertoires et compter le nombre d'items dans chaque catégorie
fileslist = glob.glob("../Start/*fasta") # stock tous les files avec l'extension .fasta
print(fileslist)
list = os.listdir() # contient la liste des fichier dans le rep files
# Comptage du nombre d'items dans chaque catégorie
nb_files = 0
nb_dirs = 0
nb_fasta = 0

for item in list:
    if os.path.isfile(item):
        nb_files += 1
        if item.endswith(".fasta"):
            nb_fasta += 1
    elif os.path.isdir(item):
        nb_dirs += 1

# 7. Générer le fichier "report_fasta.txt"
with open("report_fasta.txt", "w") as report_file:
    # Écriture du résumé
    report_file.write("#summary\n")
    report_file.write(f"nb_items = {len(list)}\n")
    report_file.write(f"nb_files = {nb_files}\n")
    report_file.write(f"nb_fasta = {nb_fasta}\n")
    report_file.write(f"nb_dir = {nb_dirs}\n")
    report_file.write("\n")
    # Détails pour chaque fichier fasta
    report_file.write("#details for each fasta file\n")
    report_file.write("# filename # size #nb_lines #nb_seq\n")
    for fasta_file in list:
        extensions_fasta = [".fasta", ".fsa", ".fa"]
        if any(fasta_file.endswith(ext) for ext in extensions_fasta):
            with open(fasta_file, "r") as f:
                lines = f.readlines()
                num_lines = len(lines)
                num_seqs = sum(1 for line in lines if line.startswith(">"))
                file_size = os.path.getsize(fasta_file)
                report_file.write(f"{fasta_file} {file_size}B {num_lines} {num_seqs}\n")


