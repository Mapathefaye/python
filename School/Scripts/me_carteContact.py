import math

# Importer les fonctions du module de parsing PDB
from S2_S3_structureToolsM1BIBS import PDB_parser

# Définir les fonctions de calcul des distances et de génération de la carte de contact

def distance_points(point1, point2):
    """ Calcule la distance entre deux points dans l'espace """
    x_diff = point1[0] - point2[0]
    y_diff = point1[1] - point2[1]
    z_diff = point1[2] - point2[2]
    return math.sqrt(x_diff**2 + y_diff**2 + z_diff**2)


def compute_contact_map(dPDB, mode="atom"):
    """ Calcule la carte de contact résidu-résidu """
    contact_map = {}
    # Loop through all residues
    for chain1 in dPDB["chains"]:
        for res1 in dPDB[chain1]["reslist"]:
            for chain2 in dPDB["chains"]:
                for res2 in dPDB[chain2]["reslist"]:
                    # Calculate distance between residues
                    distance = compute_distance(dPDB[chain1][res1], dPDB[chain2][res2], mode)
                    # Add distance to contact map
                    contact_map[(chain1, res1, chain2, res2)] = distance
    return contact_map

def compute_distance(res1, res2, mode="atom"):
    """ Calcule la distance entre deux résidus """
    # Mode atom: calculer la distance entre les atomes les plus proches
    if mode == "atom":
        distances = []
        for atom1 in res1["atomlist"]:
            for atom2 in res2["atomlist"]:
                dist = distance_points((res1[atom1]["x"], res1[atom1]["y"], res1[atom1]["z"]),
                                       (res2[atom2]["x"], res2[atom2]["y"], res2[atom2]["z"]))
                distances.append(dist)
        return min(distances)
    

# Fonction principale
def main():
    # Demander à l'utilisateur le nom du fichier PDB et les paramètres de calcul
    infile = input("Entrez le nom du fichier PDB : ")
    mode = input("Choisissez le mode de calcul de distance (atom/centroid) : ")
    seuil = float(input("Entrez le seuil de distance pour la carte de contact : "))

    # Parser le fichier PDB
    infile = "./School/1brs.pdb"
    dPDB = PDB_parser(infile)

    # Calculer la carte de contact
    contact_map = compute_contact_map(dPDB, mode)

    # Générer et enregistrer la carte de contact dans un fichier
    with open("contact_map.txt", "w") as f:
        for key, value in contact_map.items():
            if value < seuil:
                f.write(f"Résidus en contact : {key} - Distance : {value}\n")

    print("La carte de contact a été générée avec succès !")

# Appeler la fonction principale
if __name__ == "__main__":
    main()
import numpy as np
import matplotlib.pyplot as plt

# Supposons que vous ayez déjà généré votre carte de contact sous forme de matrice nommée "contact_matrix"

def plot_contact_matrix(contact_map):
    # Créer une figure
    plt.figure(figsize=(8, 6))
    # Définir le seuil pour les distances faibles (en rouge)
    red_threshold = 3.0  
    # Utiliser pcolor pour afficher la matrice avec un colormap
    plt.pcolor(contact_map, cmap='coolwarm', vmin=0, vmax=10)  # Réglez vmin et vmax selon vos besoins
    # Définir la couleur rouge pour les valeurs inférieures au seuil
    cmap = plt.cm.get_cmap('coolwarm')
    cmap.set_under('red')
    # Ajouter une barre de couleur
    plt.colorbar(label='Distance')
    # Ajouter des titres et des étiquettes d'axe
    plt.title('Carte de contact résidu-résidu')
    plt.xlabel('Résidu j')
    plt.ylabel('Résidu i')
    # Afficher la figure
    plt.show()

# Exemple d'utilisation avec une matrice de contact fictive
contact_map = np.random.rand(10, 10)  # Matrice de contact aléatoire pour l'exemple
plot_contact_matrix(contact_map)


