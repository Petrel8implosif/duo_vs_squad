import numpy as np

def lire_matrice(fichier):
    with open(fichier, 'r') as f:
        lignes = f.readlines()[0:5]  
    matrice = np.array([list(map(float, ligne.split())) for ligne in lignes])
    return matrice

def calculer_valeurs_propres(matrice):
    valeurs_propres, _ = np.linalg.eig(matrice)
    return valeurs_propres

# Exemple d'utilisation
if __name__ == "__main__":

    fichier_A = 'A_origin.txt'
    matrice_A = lire_matrice(fichier_A)
    valeurs_propres_A = calculer_valeurs_propres(matrice_A)
    print("Les valeurs propres de la matrice A_origin sont:", valeurs_propres_A)