import numpy as np
from param import *
from tools import *


def reception_matrix_conv(Map_longeur, Map_temperature,para= Param()):
    """
   @input: Matrix Map_longeur | Matrice des longeurs de flamme
           Matrix Map_temperature | Matrice des temperatures
          float | dt temps élementaire avec lequel on travaille

   @output: Matrix Matrice des energies recues pour toutes les cellules

   @utilite: Est utilise pour calculer les convection

   """
    l = len(Map_longeur)
    result_matrix = np.zeros(shape=(l, l))
    for i in range(0, l):
        for j in range(0, l):
            if Map_longeur[i, j] != 0:
                result_matrix = np.add(sender_conv(i, j, Map_temperature, Map_longeur[i, j],para), result_matrix)
    return result_matrix


def sender_conv(i, j, Map_temperature, longeur_flamme,para= Param()):
    """
   @input:  Matrix Map_temperature | Matrice des temperatures de flamme
            int i | Indice en ligne de la cellule emettrice
            int j | Indice en colonne de la cellule emettrice
            float | temperature Rayon de la cellule emettrice
            float | longeur_flamme de la cellule emettrice
            float | dt temps élementaire avec lequel on travaille


   @output: Matrix Matrice des energies recues en prevenance de (i,j)

   @utilite: Fonction utilise par reception_matrix_conv
   """
    l = len(Map_temperature)
    R = para.D/2
    result_matrix = np.zeros(shape=(l, l))
    i = l-1-i
    depart_x = (j*2*R + R)
    depart_y = (i*2*R + R)
    vecteur_directeur_u_x = np.cos(para.Teta_horizontal_vent)*para.U
    vecteur_directeur_u_y = np.sin(para.Teta_horizontal_vent)*para.U
    teta_par_raport_horizontale = para.Teta_horizontal_vent
    point_depart_u = (depart_x, depart_y, 0.5)
    vecteur_u = (vecteur_directeur_u_x, vecteur_directeur_u_y, 0)

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
    if vecteur_directeur_u_x != 0:
        tx = [((n*2*R-depart_x)/vecteur_directeur_u_x) for n in range(0,l+1) if ((n*2*R-depart_x)/vecteur_directeur_u_x) >= 0]
    if vecteur_directeur_u_y !=0:
        ty = [((n*2*R-depart_y)/vecteur_directeur_u_y) for n in range(0,l+1) if ((n*2*R-depart_y)/vecteur_directeur_u_y) >= 0]
    tx = sorted(tx, reverse=True)
    ty = sorted(ty, reverse=True)

    # On trie tx et ty et pour avancer de case en case on regarde si on a intersection avec un axe vertical ou horizontal
    # Pour cela  on compare le min des deux ensembles ce qui revient a comparer le premier element de chaque liste
    # puisque les listes sont tries.
    iii, jjj = i, j 
    distance_parcourue_jusque_intersection = -1
    while (tx != [] or ty != []) :
        if tx == []:
            ty.pop()
            iii += incr_i
        elif  ty == []:
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
        hauteur_cylindre = 1
        intersec,distance = intersection(point_depart_u, vecteur_u, coordonne_cylindre, rayon_cylindre, hauteur_cylindre)
        if intersec:
            distance_parcourue_jusque_intersection = distance
            # ON AJOUTE ENERGIE TRANSPORTE PAR CONVECTION A CELLULE III JJJ
            duv = distance_parcourue_jusque_intersection
            v = para.cinematic_viscosity_air((Map_temperature[l-1-i, j]+para.Tambiant)/2)
            k = para.thermal_conductivity_air((Map_temperature[l-1-i, j]+para.Tambiant)/2)
            Re = para.U * duv/v
            h = 0.037*k*(Re**(0.8))*(para.Pr**(1/3))/duv
            conv_puissance = (h/para.delta)*(Map_temperature[l-1-i, j]-Map_temperature[l-iii-1, jjj])*np.exp((-1)*0.3*duv/longeur_flamme)
            result_matrix[l-iii-1, jjj] = result_matrix[l-iii-1, jjj] + conv_puissance

    return result_matrix

