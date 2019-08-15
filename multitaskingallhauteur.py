# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 19:23:23 2019

"""
import multiprocessing 
from multiprocessing import Process, current_process, Queue, Array
import numpy as np
import random
import datetime
from tqdm import tqdm
import sys
import traceback


def intersection(point_depart_photon,vecteur_photon,coordonne_cylindre,rayon_cylindre,hauteur_cylindre):
    """
   @input:  (float,float,float) point_depart_photon | Point de depart du photon
            (float,float,float) vecteur_photon | Vecteur associe au photon
            (float,float) coordonne_cylindre | Coordonnes du centre du cylindre
            float rayon_cylindre | Rayon du cylindre
            float hauteur_cylindre | Hauteur du cylindre


   @output: (bool,float) | Nous indique si intersection et la distance

   @utilite: Fonction utilise pour le calcul des radiations
   """

    px,py,pz = point_depart_photon
    xd,yd,zd = vecteur_photon
    xc,yc = coordonne_cylindre

    a = (xd**2)+(yd**2)
    b= (2*(px-xc)*xd)+(2*(py-yc)*yd)
    c = ((px-xc)**2) + ((py-yc)**2) - (rayon_cylindre**2)
    racines = np.roots([a,b,c])
    
    racine = [x for x in racines if x >= 0 and np.isreal(x)]
    if len(racine) == 0:
        return (False,-1)
    if len(racine) ==1 :
        return (False,-1)
         #cas impossible

    dist = (lambda x,y,z : np.sqrt((x**2)+(y**2)+(z**2)))
    if len(racine)==2:
        z1 = pz+racine[0]*zd
        z2 = pz+racine[1]*zd
        t1 = racine[0]
        t2 = racine[1]

        if (0 < z1 and z1 < hauteur_cylindre):
            if (0 < z2 and z2 < hauteur_cylindre):
                mm = min(t1,t2)
                return (True,dist((px+mm*xd-px),(py+mm*yd-py),(pz+mm*zd-pz)))
            else:
                if z2 >= hauteur_cylindre:
                    t2 = (hauteur_cylindre-pz)/zd
                else:
                    t2 = (0-pz)/zd
        else:
             if (0 < z2 and z2 < hauteur_cylindre):
                if z1 >= hauteur_cylindre:
                    t1 = (hauteur_cylindre-pz)/zd
                else:
                    t1 = (0-pz)/zd
             else:
                if z1 < hauteur_cylindre and z2 < hauteur_cylindre:
                    return (False,-1)
                if z1 > hauteur_cylindre and z2 > hauteur_cylindre:
                    return (False,-1)
                pass#TODO cas possible
        mm = min(t1,t2)
        return (True,dist((px+mm*xd-px),(py+mm*yd-py),(pz+mm*zd-pz)))
    
    
def sender(Map_hauteur, i, j, Map_energie, R, dteta, dh, nphoton, hauteur_cellule_combustible,debug=True,dprint=True):
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
    l = len(Map_hauteur)
 
    # On récupère les infos sur la case qui va emettre
    energie = Map_hauteur[i, j]
    nombre_subdivision_h = int(Map_hauteur[i, j]/dh)
    
    i = l-i-1  # On modife le i pour qu'il corresponde au schéma ci dessus
    
    # On calcule des valeurs pour notre ds et l'énergie élémentaire d'un photon
    nombre_subdivision_teta = int((2*np.pi)/dteta)
    photon_energie = energie/(nombre_subdivision_h*nombre_subdivision_teta*nphoton)
    result_matrix = np.zeros(shape=(l,l))
    if dprint:
        print("------------------------------")
        print("Nombre de surfaces elementaires considerees {}".format(str(nombre_subdivision_h*nombre_subdivision_teta)))
    for  u in range(nombre_subdivision_h):
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
            hauteur_envoi = dh*(u+1)
            tetads = dteta*v
            depart_x = (j*2*R + R) + R*np.cos(tetads)
            depart_y = (i*2*R + R) + R*np.sin(tetads)
            depart_z = hauteur_envoi
            if u == 0 and debug:
                xx.append(depart_x)
                yy.append(depart_y)
            uu = []
            vv = []

            for counter in range(nphoton):

                # On genere deux angles randoms/ rand_alpha doit etre genere de maniere particuliere
#                rand_alpha = np.arccos(random.random()*2 -1)
#                rand_phi = (random.random())*np.pi -(np.pi/2)
                rand_phi = 2*np.arcsin(np.sqrt(random.random()))
                rand_alpha = random.random()*np.pi
                #  On regarde dans quelle direction part notre photon
                vecteur_directeur_photon_x = np.cos(rand_alpha)*np.sin(tetads)+np.sin(rand_alpha)*np.cos(tetads)
                vecteur_directeur_photon_y = (-1)*np.cos(rand_alpha)*np.cos(tetads)+np.sin(rand_alpha)*np.sin(tetads)
                vecteur_directeur_photon_z = np.cos(rand_phi)

                if u == 0 and debug:
                    uu.append(vecteur_directeur_photon_x)
                    vv.append(vecteur_directeur_photon_y)
                    
                #  On regarde a quel angle cette direction fait avec axe des x (axe des j)
                teta_par_raport_horizontale = np.arctan2(vecteur_directeur_photon_y, vecteur_directeur_photon_x)
                # VERIF DE CET ANGLE EFFECTUEE: donne des valeurs entre -np.pi et np.pi
                
                
                point_depart_photon = (depart_x, depart_y, depart_z)
                vecteur_photon = (vecteur_directeur_photon_x, vecteur_directeur_photon_y, vecteur_directeur_photon_z)
               

                #  On determine dans quel quart de plan le photon se dirige

                (incr_i, incr_j)= 1, 1 #  Haut Droite
                if teta_par_raport_horizontale > np.pi/2 and teta_par_raport_horizontale < np.pi:
                    (incr_i, incr_j)= 1, -1  #  Haut Gauche
                if teta_par_raport_horizontale < 0 and teta_par_raport_horizontale > -np.pi/2:
                    (incr_i, incr_j)= -1, 1  #  Bas Droite
                if teta_par_raport_horizontale > -np.pi and teta_par_raport_horizontale < -np.pi/2:
                    (incr_i, incr_j)= -1, -1  #  Bas Droite
                if teta_par_raport_horizontale == 0:
                    (incr_i, incr_j) = 0, 1
                if teta_par_raport_horizontale == np.pi or teta_par_raport_horizontale == -np.pi:
                    (incr_i, incr_j) = 0, -1
                if teta_par_raport_horizontale == np.pi/2:
                    (incr_i, incr_j) = 1, 0
                if teta_par_raport_horizontale == -np.pi/2:
                    (incr_i, incr_j) = -1, 0



                #  Ici on a deux ensembles tx et ty qui correspondent aux parametres a entrer dans les equations
                #  parametrique pour tomber sur une intersections avec le quadrillage selon x ou selon y.
                tx = []
                ty = []
                if vecteur_directeur_photon_x != 0:
                    tx = [((n*2*R-depart_x)/vecteur_directeur_photon_x) for n in range(0,l+1) if ((n*2*R-depart_x)/vecteur_directeur_photon_x) >= 0]
                if vecteur_directeur_photon_y !=0:
                    ty = [((n*2*R-depart_y)/vecteur_directeur_photon_y) for n in range(0,l+1) if ((n*2*R-depart_y)/vecteur_directeur_photon_y) >= 0]
                tx = sorted(tx, reverse=True)
                ty = sorted(ty, reverse=True)
                
                #  On trie tx et ty et pour avancer de case en case on regarde si on a intersection avec un axe vertical ou horizontal
                #  Pour cela  on compare le min des deux ensembles ce qui revient a comparer le premier element de chaque liste
                #  puisque les listes sont tries.
               
                iii, jjj= i, j 
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
                    if not((iii <= l-1) and (iii >= 0) and (jjj <= l-1) and (jjj >= 0)and (abs(i-iii)<10)and (abs(j-jjj)<10)):
                        break
                    
                    #Maintenant que on connait la nouvelle case que va traverser le rayon on regarde si il intersecte avec le cylindre dessus 
                    #On va pouvoir utiliser la fonction intersection
                    coordonne_cylindre = ((jjj*2*R +R), (iii*2*R+R))
                    rayon_cylindre = R
                    
                    
                    hauteur_cylindre = Map_hauteur[l-1-iii, jjj]
                    iscombustible = False
                    if hauteur_cylindre == 0:
                        iscombustible = True
                        hauteur_cylindre = hauteur_cellule_combustible
                        
                    
                    intersec,distance = intersection(point_depart_photon,vecteur_photon,coordonne_cylindre,rayon_cylindre,hauteur_cylindre)
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
                                # ON AJOUTE ENERGIE TRANSPORTE PAR PHOTON A CELLULE III JJJ
                            
                                      
                                result_matrix[l-1-iii,jjj] = result_matrix[l-1-iii,jjj] + photon_energie
                            break
                        

                if distance_parcourue_jusque_intersection != -1: #AU CAS OU LE TRAJET DU PHOTON INTERSECTE AVEC AUCUN CYLINDRE
                    continue
            
        
    return result_matrix


dteta =  2*10**-2
dh = 0.1
nphoton = 25
R = 2.54/2
h = 0.51 #Hauteur des autres cellules //IMPORTANT //
taille =50
premierehauteur = 0.51
espacement = 0.2

def worker(part):

    dteta =  2*10**-2
    dh = 0.1
    nphoton = 20
    R = 2.54/2
    h = 0.51 #Hauteur des autres cellules //IMPORTANT //
    taille =50
    premierehauteur = 0.51
    espacement = 0.2

    hauteur_flamme= round(premierehauteur+espacement*part,1)
    print("Worker {} , started ".format(hauteur_flamme))

    Map_energie = np.zeros(shape=(taille,taille))+1

  
    Map_hauteur = np.zeros(shape=(taille, taille))
    Map_hauteur[taille//2][taille//2] = hauteur_flamme
    d= {}
    if h !=0:
        w=sender(Map_hauteur, taille//2, taille//2, Map_energie, R, dteta, dh, nphoton,h,debug=False,dprint=False)
        for i in range(taille):
            for j in range(taille):
                em_ij = np.zeros(shape=(taille, taille))
                for u in range(taille):
                    for v in range(taille):
                        if taille//2+abs(i-u) >0 and  taille//2+abs(i-u) < taille and taille//2+abs(j-v) >0 and taille//2+abs(j-v) <taille:
                            em_ij[u,v]=w[taille//2+abs(i-u),taille//2+abs(j-v)]
                d[(i,j)] = em_ij
    else:
        d={}
	print("\n")
    print("Worker {} , finished \n".format(hauteur_flamme)) 
	print("\n")	
    np.save('{}_50x50.npy'.format(part),d)         

def listener(q,nbr):
    pbar = tqdm(total = nbr)
    sys.stdout.flush()
    for item in iter(q.get, None):
        pbar.update()
    sys.stdout.flush()



if __name__ == '__main__':
    print("-------------------------------")
    print(" Précalcul des radiations ")
    print("-------------------------------")
    print("Taille de la grille: {}x{} ".format(taille,taille))
    print("Rayon cylindres: {} ".format(R))
    print("Hauteur combustible: {} ".format(h))
    print("dh: {} ".format(dh))
    print("dteta: {} ".format(dteta))
    print("nphoton: {} ".format(nphoton))
    print("-------------------------------")
    print("Heure début: {}".format(datetime.datetime.now()))
    print("-------------------------------")
    manager = multiprocessing.Manager()
    jobs = []
    y = datetime.datetime.now()

    with multiprocessing.Pool(processes=6) as p:
        max_ = 50
        with tqdm(total=max_) as pbar:
            for i, _ in tqdm(enumerate(p.imap_unordered(worker, range(0, max_)))):
                print("\n")
                pbar.update()
                print("\n")

    print(datetime.datetime.now()-y)
    print("Terminé")
    input("Close2")
    
    





    
