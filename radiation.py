import numpy as np 
import random 
import matplotlib as plt
from tools import *


def sender(Map_hauteur, i, j, Map_energie, parameters, debug=True, dprint=True):
    """
   @input:  Matrix Map_hauteur | Matrice des hauteur de flamme
            Matrix Map_energie | Matrice des energies a émettre
            int i | Indice en ligne de la cellule emettrice (les cellules etant numerotes de 0 a N-1)
            int j | Indice en colonne de la cellule emettrice
            float R | Rayon de la cellule emettrice
            float H | Hauteur de la cellule emettrice
            float dteta | Angle elementaire considere pour les calculs
            float dh | Hauteur elementaire considere pour les calculs
            int nphoton | Nombre de photons emit par surface dS
            float hauteur_cellule_combustible | Hauteur combustible

   @output: Matrix Matrice des energies recues en provenance de (i,j)

   @limite: On suppose que toutes les cellules de la grille sont de rayon R
            On suppose que toute les cellules combustible ont la meme hauteur

   @schema: Voici la situation

                              /                          ^ey
                             / rand_alpha                |
                            /                            |i    j
                           /  ..    ..    [i*R,j*R]             --> ex
                       ..  \            ..
                            \
                   ..        \              ..
                              \   teta_ds
                 ..            \              ..
                                O - - - - - - -
                 ..                           ..

                  ..                         ..


   @utilite: Fonction utilise par reception_matrix
   """
    l = parameters.Taille_grille

    # On récupère les infos sur la case qui va emettre
    energie = Map_hauteur[i, j]
    nombre_subdivision_h = int(Map_hauteur[i, j]/parameters.dh)

    i = l-i-1  # On modife le i pour qu'il corresponde au schéma ci dessus

    # On calcule des valeurs pour notre ds et l'énergie élémentaire d'un photon
    nombre_subdivision_teta = int((2*np.pi)/parameters.dteta)
    photon_energie = energie/(nombre_subdivision_h*nombre_subdivision_teta*parameters.nphoton)
    result_matrix = np.zeros(shape=(l, l))
    if dprint:
        print("------------------------------")
        print("Nombre de surfaces elementaires considerees {}".format(str(nombre_subdivision_h*nombre_subdivision_teta)))


    for u in range(nombre_subdivision_h):
        if u == int(nombre_subdivision_h/4) and dprint:
            print("25% ..")
        if u == int(nombre_subdivision_h/2) and dprint:
            print("50% ..")
        if u == int(3*nombre_subdivision_h/4) and dprint:
            print("75% ..")
        if u == 0 and debug:
            xx = []
            yy = []
        for v in range(nombre_subdivision_teta):
            hauteur_envoi = parameters.topography[i][j]+dh*(u+1)  # TODO Crochets pas sur du tout
            tetads = dteta*v
            depart_x = (j*2*R + R) + R*np.cos(tetads)
            depart_y = (i*2*R + R) + R*np.sin(tetads)
            depart_z = hauteur_envoi
            if u == 0 and debug:
                xx.append(depart_x)
                yy.append(depart_y)
            uu = []
            vv = []

            for counter in range(parameters.nphoton):

                # On genere deux angles randoms/ rand_alpha doit etre genere de maniere particuliere
                rand_phi = 2*np.arcsin(np.sqrt(random.random()))
                rand_alpha = random.random()*np.pi

                #  On regarde dans quelle direction part notre photon
                vecteur_directeur_photon_x = np.cos(rand_alpha)*np.sin(tetads)+np.sin(rand_alpha)*np.cos(tetads)
                vecteur_directeur_photon_y = (-1)*np.cos(rand_alpha)*np.cos(tetads)+np.sin(rand_alpha)*np.sin(tetads)
                vecteur_directeur_photon_z = np.cos(rand_phi)

                if u == 0 and debug:
                    uu.append(vecteur_directeur_photon_x)
                    vv.append(vecteur_directeur_photon_y)

                # On regarde a quel angle cette direction fait avec axe des x (axe des j)
                teta_par_raport_horizontale = np.arctan2(vecteur_directeur_photon_y, vecteur_directeur_photon_x)
                point_depart_photon = (depart_x, depart_y, depart_z)
                vecteur_photon = (vecteur_directeur_photon_x, vecteur_directeur_photon_y, vecteur_directeur_photon_z)

                # On determine dans quel quart de plan le photon se dirige
                (incr_i, incr_j) = 1, 1  # Haut Droite
                if teta_par_raport_horizontale > np.pi/2 and teta_par_raport_horizontale < np.pi:
                    (incr_i, incr_j) = 1, -1  # Haut Gauche
                if teta_par_raport_horizontale < 0 and teta_par_raport_horizontale > -np.pi/2:
                    (incr_i, incr_j) = -1, 1  # Bas Droite
                if teta_par_raport_horizontale > -np.pi and teta_par_raport_horizontale < -np.pi/2:
                    (incr_i, incr_j) = -1, -1  # Bas Droite
                if teta_par_raport_horizontale == 0:
                    (incr_i, incr_j) = 0, 1
                if teta_par_raport_horizontale == np.pi or teta_par_raport_horizontale == -np.pi:
                    (incr_i, incr_j) = 0, -1
                if teta_par_raport_horizontale == np.pi/2:
                    (incr_i, incr_j) = 1, 0
                if teta_par_raport_horizontale == -np.pi/2:
                    (incr_i, incr_j) = -1, 0

                # Ici on a deux ensembles tx et ty qui correspondent aux parametres a entrer dans les equations
                # parametrique pour tomber sur une intersections avec le quadrillage selon x ou selon y.
                tx = []
                ty = []
                if vecteur_directeur_photon_x != 0:
                    tx = [((n*2*R-depart_x)/vecteur_directeur_photon_x) for n in range(0, l+1) if ((n*2*R-depart_x)/vecteur_directeur_photon_x) >= 0]
                if vecteur_directeur_photon_y !=0:
                    ty = [((n*2*R-depart_y)/vecteur_directeur_photon_y) for n in range(0, l+1) if ((n*2*R-depart_y)/vecteur_directeur_photon_y) >= 0]
                tx = sorted(tx, reverse=True)
                ty = sorted(ty, reverse=True)
                
                # On trie tx et ty et pour avancer de case en case on regarde si on a intersection avec un axe vertical ou horizontal
                # Pour cela  on compare le min des deux ensembles ce qui revient a comparer le premier element de chaque liste
                # puisque les listes sont tries.
                iii, jjj = i, j
                distance_parcourue_jusque_intersection = -1
                while (tx != [] or ty != []):
                    if tx == []:
                        ty.pop()
                        iii += incr_i
                    elif ty == []:
                        tx.pop()
                        jjj += incr_j
                    else:
                        if tx[-1] > ty[-1]:
                            ty.pop()
                            iii += incr_i
                        elif tx[-1] < ty[-1]:
                            tx.pop()
                            jjj += incr_j
                        else:
                            ty.pop()
                            tx.pop()
                            jjj += incr_j
                            iii += incr_i
                    if not((iii <= l-1) and (iii >= 0) and (jjj <= l-1) and (jjj >= 0)):
                        break

                    # Maintenant que on connait la nouvelle case que va traverser le rayon on regarde si il intersecte avec le cylindre dessus
                    # On va pouvoir utiliser la fonction intersection
                    coordonne_cylindre = ((jjj*2*R +R), (iii*2*R+R))
                    rayon_cylindre = R
                    hauteur_cylindre = Map_hauteur[l-1-iii, jjj]
                    iscombustible = False
                    if hauteur_cylindre == 0:
                        iscombustible = True
                        hauteur_cylindre = hauteur_cellule_combustible
                    intersec, distance = intersection(point_depart_photon, vecteur_photon, coordonne_cylindre, rayon_cylindre, hauteur_cylindre,parameters.topography)
                    if intersec:
                        distance_parcourue_jusque_intersection = distance
                        if not(iscombustible):
                            break
                        else:
                            # CALCUL PROBA POUR VOIR SI LE PHOTON A PAS ETE ABSORBE PAR AIR
                            # on utilise la loi approchée (modele SNB) proba = a+ b*(l^c) a=1.213 , b =-0.253 , c=0.17 pour 25% d'humidité
                            proba_ab = 1.213 -0.253*(distance**(0.17))
                            if random.random() < proba_ab:
                                # CAS NON ABSORBE PAR L'AIR
                                # ON AJOUTE ENERGIE TRANSPORTE PAR PHOTON A LA CELLULE III JJJ
                                result_matrix[l-1-iii, jjj] = result_matrix[l-1-iii, jjj] + photon_energie
                            break
                        
                if distance_parcourue_jusque_intersection != -1:  # AU CAS OU LE TRAJET DU PHOTON INTERSECTE AVEC AUCUN CYLINDRE
                    continue
 
    return result_matrix



def reception_matrix(para,rad, Map_energie):
    """
   @input:  Object para | Objet contenant les parametres
            Matrix Map_energie | Matrice des energies emises

   @output: Matrix Matrice des energies recues pour toutes les cellules

   @utilite: Est utilise pour calculer les radiations
   
   @valide du calcul: energie emise ne varie pas pendant la vie d'une flamme

   """
    
    result_matrix = np.zeros(shape=(para.Taille_grille, para.Taille_grille))
    for i in range(0, l):
        for j in range(0, l):
                mat = rad[i][j]*Map_energie[i][j])
                result_matrix = np.add(result_matrix,mat)

    return result_matrix


def test_sender():
    """
    @utilite:
    Cette fonction a pour seul but de verifier si la fonction sender
    fonctionne correctement

    """
    Map_hauteur = np.zeros(shape=(5, 5))

    Map_hauteur[1, 0] = 2
    a = sender(Map_hauteur, 1, 0, 100000, 2, 0.5*10**-1, 0.4*10**-1, 20, 1)
    plt.matshow(a)
    plt.show()


def test_receptionmatrix(printvaleur=True):
    Map_hauteur = np.zeros(shape=(50, 50))
    Map_energie = np.zeros(shape=(50, 50))
    Map_hauteur[25, 25] = 5
    Map_energie[25, 25] = 100000000
    # Map_hauteur[5,6] =1
    # Map_energie[5, 6] =100000
    R = 2
    dteta = 5*10**-2
    dh = 0.5*10**-1
    nphoton = 20
    mat = sender(Map_hauteur, 25, 25, Map_energie, R, dteta, dh, nphoton, 1)

    def f(x):
        return np.log10(x)

    f = np.vectorize(f)
    np.set_printoptions(precision=3)
    print(mat)
    fig, ax = plt.subplots()
    ax.matshow(f(mat), cmap=plt.get_cmap("afmhot"))
    if printvaleur:
        for (i, j), z in np.ndenumerate(mat):
            ax.text(j, i, '{:0.01f}'.format(z), ha='center', va='center')
    plt.savefig('uniform_high_res.png', format='svg', dpi=2000)
    fig.show()
