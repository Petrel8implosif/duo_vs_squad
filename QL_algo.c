#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPSILON 1e-14

// Algorithme TQL2 pour matrice symétrique tridiagonale (adapté de Numerical Recipes in C)
// d[] : éléments diagonaux
// e[] : éléments sous-diagonaux (e[0] n'est pas utilisé)
// Le résultat (les valeurs propres) est écrit dans d[]
void tql2(int n, double d[], double e[]) {
    int l, m, i, iter;
    double s, r, p, g, f, dd, c, b;
    // Décaler e vers la gauche : on ne considère pas e[0]
    for (i = 1; i < n; i++) {
        e[i-1] = e[i];
    }
    e[n-1] = 0.0;
    
    for (l = 0; l < n; l++) {
        iter = 0;
        do {
            for (m = l; m < n-1; m++) {
                dd = fabs(d[m]) + fabs(d[m+1]);
                if (fabs(e[m]) + dd == dd)
                    break;
            }
            if (m != l) {
                if (iter++ == 30) {
                    printf("Too many iterations in tql2\n");
                    break;
                }
                g = (d[l+1] - d[l]) / (2.0 * e[l]);
                r = sqrt(g * g + 1.0);
                // Calculer le décalage implicite
                g = d[m] - d[l] + e[l] / (g + (g >= 0 ? fabs(r) : -fabs(r)));
                s = 1.0;
                c = 1.0;
                p = 0.0;
                for (i = m - 1; i >= l; i--) {
                    f = s * e[i];
                    b = c * e[i];
                    r = sqrt(f * f + g * g);
                    e[i+1] = r;
                    if (r == 0.0) {
                        d[i+1] -= p;
                        e[m] = 0.0;
                        break;
                    }
                    s = f / r;
                    c = g / r;
                    g = d[i+1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    p = s * r;
                    d[i+1] = g + p;
                    g = c * r - b;
                }
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        } while (m != l);
    }
}

// Fonction simple de tri (tri à bulles) pour ranger les valeurs propres
void sort_eigenvalues(int n, double d[]) {
    for (int i = 0; i < n-1; i++) {
        for (int j = i+1; j < n; j++) {
            if (d[i] > d[j]) {
                double tmp = d[i];
                d[i] = d[j];
                d[j] = tmp;
            }
        }
    }
}


void analytical_values(int n, double d[]) {
    double dx = 2.0 / (n + 1);
    for (int i = n; i > 0; i--) {
        d[n-i] = -4.0 /(dx * dx)  * (sin((i) * M_PI / (2 * (n + 1))) * sin((i) * M_PI / (2 * (n + 1))));
    }
}

int main() {
    int n = 2; // Nombre de noeuds internes
    double dx = 2.0 / (n + 1); // Pas de discrétisation sur [-1, 1]
    // On construit la matrice tridiagonale du Laplacien discrétisé
    // Pour le Laplacien, sur chaque noeud, la diagonale vaut -2/(dx^2)
    // et la sous-diagonale vaut 1/(dx^2)
    double *d = (double *)calloc(n, sizeof(double)); // diagonale
    double *e = (double *)calloc(n, sizeof(double)); // sous-diagonale (e[0] non utilisé)
    for (int i = 0; i < n; i++) {
        d[i] = -2.0 / (dx * dx);
        if (i < n - 1) {
            e[i+1] = 1.0 / (dx * dx);
        }
    }
    
    // Calcul des valeurs propres avec l'algorithme TQL2
    tql2(n, d, e);
    
    // Tri des valeurs propres (optionnel, pour comparer avec la formule analytique rangée par ordre croissant)
    sort_eigenvalues(n, d);
    
    printf("Eigenvalues computed by TQL2:\n");
    for (int i = 0; i < n; i++) {
        printf("%f\n", d[i]);
    }
    
    // Calcul des valeurs propres analytiques
    double *d_analytical = (double *)calloc(n, sizeof(double));
    analytical_values(n, d_analytical);

    printf("\nAnalytical eigenvalues:\n");
    for (int i = 0; i < n; i++) {
        printf("%f\n", d_analytical[i]);
    }

    // Différence entre les valeurs propres calculées et analytiques
    printf("\nDifference between analytical and computed eigenvalues:\n");
    for (int i = 0; i < n; i++) {
        printf("%f\n", d[i] - d_analytical[i]);
    }

    
    free(d_analytical);
    free(d);
    free(e);
    return 0;
}
