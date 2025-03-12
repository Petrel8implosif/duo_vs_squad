import numpy as np

def lire_matrice(fichier):
    with open(fichier, 'r') as f:
        lignes = f.readlines()[1:6]  # Ignorer la première ligne
    matrice = np.array([list(map(float, ligne.split())) for ligne in lignes])
    return matrice

def calculer_valeurs_propres(matrice):
    valeurs_propres, _ = np.linalg.eig(matrice)
    return valeurs_propres

# Exemple d'utilisation
if __name__ == "__main__":
    fichier = 'matrix_Q.txt'
    matrice = lire_matrice(fichier)
    valeurs_propres = calculer_valeurs_propres(matrice)
    print("Les valeurs propres de la matrice sont:", valeurs_propres)