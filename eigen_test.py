import numpy as np

def lire_matrice(fichier):
    with open(fichier, 'r') as f:
        lignes = f.readlines()
    matrice = np.array([list(map(float, ligne.split())) for ligne in lignes])
    return matrice

def calculer_valeurs_propres(matrice):
    valeurs_propres, _ = np.linalg.eig(matrice)
    return valeurs_propres

# Exemple d'utilisation
if __name__ == "__main__":

    fichier_A = 'A_devoir.txt'
    matrice_A = lire_matrice(fichier_A)
    valeurs_propres_A = calculer_valeurs_propres(matrice_A)
    print("Les valeurs propres de la matrice A_devoir sont:", valeurs_propres_A)