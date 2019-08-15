# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 13:02:55 2018

"""

from enum import Enum
import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib.colors import colorConverter
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import datetime
import os
"""
06/04/19 Ammélioration de la drastique du temps de la calcul

"""




"""
Importer les radiations précalculées
Mettre les fichiers de la forme: entier_nom.npy dans le meme dossier
que ce fichier
"""
Taille_grille = 50
Nomrad = "50x50"
Init_h = 0.51
Espacement_H = 0.2
Nombre_hauteur_differentes = 50
Precal_path = "ray051_02_2/" 


"""
Creation des differents objets
"""

class Param(object):
    def  __init__(self):
        """ Constantes & Variables que l'on peut ajuster """
        self.Vchar = 0.3  # Contenu en résidu charbonneux
        self.Kf = 0.4  # Coefficient d'extinction de la flamme,  utile pour calculer l'émissivité de celle ci
        self.Xr = 0.35  # Fraction de Qpoint libérée par rayonnement (entre 0.2-0.4 depend végétaux) | sans unité
        self.Tambiant = 307  # Température de l'air | K
        self.Tpyro = 500  # Température de pyrolyse | K
        self.sigma_boltzmann = 5.67*(10**-8)  # Stefan-Boltzman constante | kg.s^-3 .K^-4
        self.sigma_k = 12240  # Rapport surface volume des aiguilles (entre 4000 et 12000 selon espece)| m^-1
        self.alpha_k = 0.0012  # Fraction volumique de la phase solide | sans unité
        self.delta = 4/(self.sigma_k*self.alpha_k)
        self.Efb = 0.9  # émissivité de du combustible
        self.U = 10  # Vitesse du vent | m.s^-1
        self.Teta_horizontal_vent = np.pi/3
        self.g = 9.81
        self.D = 2.54  # Diametre d'une cellule | m
        self.Tc = 5 # Temps d'extinction de la flamme | s
        self.mprimeprimeDFF = 0.313  # Charge de combustible | kg.m^-2
        self.deltahc = 15.6*(10**6)  # Chaleur de combustion |" J.kg^-1
        self.Rho_wff = 612  # Masse volumique du combustible végétal humide  | kg.m^-3
        self.Rho_dff = 720  # Masse volumique du combustible végétal sec  | kg.m^-3
        self.Cdff = 1470  # Chaleur spécifique combustible végétal sec  | J.kg^-1 .K^-1
        self.Cpeau = 4181  # Capacité de l'eau | J.kg^-1 .K^-1
        self.FMC = 0.058  # Teneur en eau initiale du combustible végétal | sans unité
        self.Lpyro = 15.6*(10**6)  # Pin 1.77 +- 0.31 MJ. kg^-1 | J.kg^-1
        self.Lvap = 2.25*(10**6)  # Chaleur latente de vaporisation de l’eau à 373K | J.kg^-1
        self.Pr = 0.707  # Nombre de Prandtl (0.707 pour air a 25°C) | sans unité
        self.hcom = 0.51 # Hauteur de la strate végétale supossé de hauteur constante
        self.absorptivity = 0.9
        self.delta = min(self.hcom,self.delta)
        self.thermal_conductivity_air = (lambda T: 1.5207*(10**(-22))*(T**(3)) - 4.8574*(10**(-8))*(T**(2))+ 1.0184*(10**(-4))*T -3.9333*(10**-4))  # WIKI
        self.cinematic_viscosity_air = (lambda T: (-1.363528*(10**-14)*T**3)+(1.00881778*(10**-10) *(T**2))+(3.452139*(10**-8)*T)-(3.400747*(10**-6)))  # WIKI
        
        """
        Constantes sources:
        
        Conductivité thermique de l'air:
        http://bouteloup.pierre.free.fr/lica/phythe/don/air/air_k_plot.pdf
        
        Lpyro table:
        https://lib.dr.iastate.edu/cgi/viewcontent.cgi?referer=https://www.google.fr/&httpsredir=1&article=1156&context=me_pubs
        """
    
    
    
class CaseType(Enum):
    """
    On enumere tout les types de cases possibles
    """
    TERRE = (1, "peru")
    VEG = (2, "green")

    def __init__(self, tid, couleur):
        self.id = tid
        self.couleur = couleur


class Case(object):

    """
    Contient toute l'information d'une case

    @input: bool feu | Definit si la case est en feu
            float temperature | Temperature de la case
            CaseType ctype | Nous indique de quel type est cette case

    @output: Case / Nouvelle case
    """

    def  __init__(self, feu=False,isdry=False,para=Param(), temperature=373, ctype=CaseType.VEG, hascombustible=True):
        mass_dry= (1/4)*para.Rho_dff*np.pi*para.D**2*para.hcom
        self.ctype = ctype
        self.feu = feu
        self.temperature = para.Tambiant
        self.isdry = isdry
        self.hascombustible = hascombustible
        self.mass_dry = mass_dry
        self.mass_water = para.FMC*mass_dry
        self.produitpyrolyse = 1 - para.Vchar
        self.cFMC = para.FMC


class Grille2D(object):

    """
    Notre grille est une matrice d'objets Case

    @input: int taille | Taille de la grille

    @output: Grille | Nouvelle grille

    @fonction: update(float) -> ()
                    Utilite: Avance le temps de dt
                    Limites: si dt trop grand non cohérent
    """
    def __init__(self, taille=50,msg=True,para=Param()):
        self.taille = taille
        self.grille = []
        self.radiation_enreception = np.zeros(shape=(self.taille, self.taille))
        self.convection_enreception = np.zeros(shape=(self.taille, self.taille))
        self.time = 0
        self.rad = []
        self.loadrad()
        self.count= 0
        self.fires=[]
        self.radgrids= {}
        self.pa= para
        if msg:
            self.message_recap()
        
    def message_recap(self):
        print("--------------------------------------------")
        print("     Simulation ")
        print("--------------------------------------------")
        print("Vchar: {}".format(self.pa.Vchar))
        print("Kf: {}".format(self.pa.Kf))
        print("Xr: {}".format(self.pa.Xr))
        print("Tambiant: {} Kelvin".format(self.pa.Tambiant))
        print("Tpyro: {} Kelvin".format(self.pa.Tpyro))
        print("sigma_k: {} m^-1".format(self.pa.sigma_k))
        print("alpha_k: {}".format(self.pa.alpha_k))
        print("delta: {}".format(self.pa.delta))
        print("Efb: {}".format(self.pa.Efb))
        print("U: {} m.s^-1".format(self.pa.U))
        print("Teta_horizontal_vent: {}".format(self.pa.Teta_horizontal_vent))
        print("D: {} m".format(self.pa.D))
        print("Tc: {} s".format(self.pa.Tc))
        print("mprimeprimeDFF: {} kg.m^-2".format(self.pa.mprimeprimeDFF))
        print("Rho_wff: {}  kg.m^-3".format(self.pa.Rho_wff))
        print("Rho_dff: {}  kg.m^-3".format(self.pa.Rho_dff))
        print("Cdff: {} J.kg^-1 .K^-1".format(self.pa.Cdff))
        print("Cpeau: {} J.kg^-1 .K^-1".format(self.pa.Cpeau))
        print("FMC: {}".format(self.pa.FMC))
        print("hauteur_combustible: {}".format(self.pa.hcom))
        print("--------------------------------------------")
        print("Avant utilisation:")
        print("--------------------------------------------")
        print("Utilisez multitaskinghauteur.py pour ")
        print("générer les fichiers radiation, faire correspondre hcom avec")
        print("avec le h ligne 281 de multitaskinghauteur.py ")
        print("penser à modifier Nomrad egalement.")
        print("--------------------------------------------")
    def loadrad(self):
        for o in range(Nombre_hauteur_differentes):
            self.rad.append((np.load(Precal_path+str(o)+'_{}.npy'.format(Nomrad))))

    def image_temperature(self, Map_Temperature, fig1):
        def f(x):
            return np.log10(x)

        f = np.vectorize(f)
        ax = fig1.add_subplot(2, 2, 4)
        ax.matshow(f(Map_Temperature), cmap=plt.get_cmap("afmhot"))
        ax.xaxis.tick_bottom()
        ax.set_title('Image thermique')

    def image_situation_3D(self, Map_hauteur, fig1):
        ax1 = fig1.add_subplot(222, projection="3d")
        xpos = []
        ypos = []
        zpos = []
        dz = []
        colorl = []
        width = depth = 2
        for i in range(self.taille):
            for j in range(self.taille):
                case = self.grille[i][j]
                radius = self.pa.D/2
                height = Map_hauteur[i, j]
                if height == 0:
                    height = self.pa.hcom
                x_center = (j*2*radius + radius)
                y_center = (i*2*radius + radius)
                xpos.append(x_center)
                ypos.append(y_center)
                zpos.append(0)
                dz.append(height)
                colorl.append(case.temperature)
        cmap = cm.get_cmap('jet')
        max_height = max(colorl)  
        min_height = min(colorl)
        rgba = [cmap((k-min_height)/max_height) for k in colorl] 
        ax1.bar3d(xpos, ypos, zpos, width, depth, dz, color=rgba)
        ax1.set_title('Vue 3D')

    def image_situation(self, fig1):
        print("------------------------------")
        print("Temps : {} ".format(str(datetime.timedelta(seconds=self.time))))
        print((self.grille[23][6]).temperature)
        z = []
        for i in range(self.taille):
            o = []
            for j in range(self.taille):
                case = self.grille[i][j]
                if case.feu:
                    o.append(colorConverter.to_rgb("red"))
                elif case.mass_dry == 0:
                    o.append(colorConverter.to_rgb("black"))
                elif case.temperature == 373:
                    o.append(colorConverter.to_rgb("orange"))
                elif case.temperature > 373:
                    o.append(colorConverter.to_rgb("brown"))
                elif (case.temperature < 373) and (case.temperature > self.pa.Tambiant):
                    redco = 200*(1-(case.temperature-self.pa.Tambiant)/(373-self.pa.Tambiant))
                    c = (round(redco, 0), 100.0, 255.0)
                    o.append(c)
                else:
                    o.append(colorConverter.to_rgb(case.ctype.couleur))
            z.append(o) 
        ax1 = fig1.add_subplot(1, 2, 1)
        ax1.imshow(z)
        ax1.set_title('Image situation')
        


    def xupdate(self, times=10, timage=100000, stockage="stockage/"):
        self.count = 0
        for _ in range(times-1):
            if self.count % timage == 0:
                fig1 = plt.figure()
                self.image_situation(fig1)
                self.update(imagesituation3D=True,figs=fig1)
                path = stockage+str(self.count)+'_u.png'
                if os.path.isfile(path):
                    os.remove(path)
                fig1.savefig(path, format='png')
                
            else:
                self.update(imagesituation3D=False)
            self.count += 1

    def update(self, dt=0.1, debug=False, imagesituation3D=False, figs=None):
        """
        Protocole:
            But: Calculer nouvelles températures, nouvelles masses en eau
                 + préparer certains calculs pour la prochaine update
            Etapes:
                Pour chaque case on fait un bilan énergétique, conformément
                à celui proposé par l'article.
            Resultat:
                - self.grille est mis à jour
                - self.radiation_enreception & self.convection_enreception
                  préparés pour la prochaine update.

        """
        self.time += dt
        hauteurflamme_matrix = np.zeros(shape=(self.taille, self.taille))
        longeurflamme_matrix = np.zeros(shape=(self.taille, self.taille))
        temperature_matrix = np.zeros(shape=(self.taille, self.taille))
        energie_matrix = np.zeros(shape=(self.taille, self.taille))
        
        for i in range(self.taille):
            for j in range(self.taille):
                caseij = self.grille[i][j]
                if not(caseij.hascombustible):
                    continue
                
                #######################################################
                temperature_matrix[i, j] = caseij.temperature
                sum_qrad_entrant = self.radiation_enreception[i, j]
                sum_qconv_entrant = self.convection_enreception[i, j]

                # Calcul de qrad- 
                qradmoins = self.pa.Efb*self.pa.sigma_boltzmann*((caseij.temperature**4) - (self.pa.Tambiant**4))/self.pa.delta
                
                # Bilan de puissance
                Volume = np.pi*min(self.pa.hcom,self.pa.delta)*((self.pa.D/2)**2)
                Puissance_rechauf = (sum_qrad_entrant+sum_qconv_entrant-qradmoins)*Volume
         
                # Avant de calculer la new temperature on regarde l energie degage par les flammes
                # Calcul préparatoire de l'update suivante
                if caseij.feu:
                    Qpoint = (self.pa.mprimeprimeDFF/self.pa.Tc)*self.pa.deltahc*(np.pi*(self.pa.D**2))/4
                    hflamme_no_vent = 0.0148 * (Qpoint)**(2/5) + (-1)*1.02*self.pa.D
                    hflamme = hflamme_no_vent * ((1+4*((self.pa.U**2)/(self.pa.g*hflamme_no_vent)))**((-1)*0.5))
                    hauteurflamme_matrix[i, j] = hflamme
                    longeur_flamme = hflamme_no_vent
                    longeurflamme_matrix[i, j] = longeur_flamme
                    Sf = np.pi * self.pa.D*longeur_flamme
                    Pfprime = self.pa.Xr*Qpoint/Sf
                    energie_matrix[i, j] = self.pa.absorptivity*(Pfprime/Volume)*Sf
                    caseij.temperature = (Pfprime/(self.pa.sigma_boltzmann*(1-np.exp(-self.pa.Kf*longeur_flamme))))**(1/4)
         
                    temperature_matrix[i, j] = caseij.temperature
                #######################################################
                # On calcule nos nouvelles temperatures et nouvelles masses

                Cwff = (1/(1+self.pa.FMC))*self.pa.Cdff + self.pa.FMC/(1+self.pa.FMC)*self.pa.Cpeau
                if not(caseij.isdry):
                    if caseij.temperature <373:
                        deltaTemperature = (Puissance_rechauf) * (1 / (self.pa.Rho_wff * self.pa.alpha_k * Cwff)) * dt
                        caseij.temperature = caseij.temperature + deltaTemperature
                        if caseij.temperature > 373:
                            caseij.temperature = 373
                    else:
                        dFMC = (Puissance_rechauf) * ((-1)/(self.pa.Rho_dff*self.pa.alpha_k*self.pa.Lvap)) * dt
                        caseij.cFMC = caseij.cFMC + dFMC
                        if caseij.cFMC < 0:
                            caseij.mass_water = 0
                            caseij.cFMC = 0
                            caseij.isdry = True
                else:
                    if caseij.temperature < self.pa.Tpyro and 373 <= caseij.temperature:
                        deltaTemperature = (Puissance_rechauf)*(1/(self.pa.Rho_dff*self.pa.alpha_k*self.pa.Cdff))*dt
                        caseij.temperature = caseij.temperature+deltaTemperature
                        if caseij.temperature >= self.pa.Tpyro:
                            caseij.temperature = self.pa.Tpyro
                    else:
                        caseij.feu = True
                        deltaproduitpyrolyse = (Puissance_rechauf)*((-1)/(self.pa.Rho_dff*self.pa.alpha_k*self.pa.Lpyro))*dt
                        caseij.produitpyrolyse = caseij.produitpyrolyse + deltaproduitpyrolyse
                        if caseij.produitpyrolyse >= 1 or caseij.produitpyrolyse <0:
                            caseij.mass_dry = 0
                            caseij.hascombustible = False
                            caseij.feu = False

        self.radiation_enreception,self.fires,self.radgrids = reception_matrix(self.rad, hauteurflamme_matrix, energie_matrix,self.radiation_enreception,self.fires,self.radgrids)
        self.convection_enreception = reception_matrix_conv(longeurflamme_matrix, temperature_matrix)
        if debug:
            print(temperature_matrix)
        if imagesituation3D:
            self.image_temperature(temperature_matrix, figs)
            self.image_situation_3D(hauteurflamme_matrix, figs)
            


"""
Quelques fonctions...
"""


def sender(Map_hauteur, i, j, Map_energie, R, dteta, dh, nphoton, hauteur_cellule_combustible, debug=True, dprint=True):
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
                    intersec, distance = intersection(point_depart_photon, vecteur_photon, coordonne_cylindre, rayon_cylindre, hauteur_cylindre)
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

def create_mat(i,j,taille,w):
    em_ij = np.zeros(shape=(taille, taille))
    for u in range(taille):
        for v in range(taille):
            if taille//2+abs(i-u) >0 and  taille//2+abs(i-u) < taille and taille//2+abs(j-v) >0 and taille//2+abs(j-v) <taille and abs(i-u)<10 and abs(j-v) <10:
                em_ij[u,v]=w[taille//2+abs(i-u),taille//2+abs(j-v)]
    return em_ij

def reception_matrix(rad, Map_hauteur, Map_energie,last,fires,radgrids):
    """
   @input: Matrix Map_hauteur | Matrice des hauteurs de flamme
           Matrix Map_energie | Matrice des energies emises

   @output: Matrix Matrice des energies recues pour toutes les cellules

   @utilite: Est utilise pour calculer les radiations
   
   @valide du calcul: energie emise ne varie pas pendant la vie d'une flamme

   """
    
    l = len(Map_hauteur)
    traite =[]
    result_matrix = np.zeros(shape=(l, l))
    for i in range(0, l):
        for j in range(0, l):
            if Map_hauteur[i, j] != 0:
                traite.append((i,j))
                if (i,j) in fires:
                    pass
               
                k = int((Map_hauteur[i, j]-Init_h)//Espacement_H)
                mat = create_mat(i,j,l, ((rad[k][()])[(l//2, l//2)])*Map_energie[i][j])
                radgrids[i,j] = mat
                result_matrix = np.add(result_matrix,mat)
    for (i,j) in fires:
        if not((i,j) in traite):
            result_matrix = np.subtract(result_matrix,radgrids[i,j])
            fires.remove((i,j))
    return result_matrix,traite,radgrids


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


def intersection(point_depart_photon, vecteur_photon, coordonne_cylindre, rayon_cylindre, hauteur_cylindre):
    """
   @input:  (float,float,float) point_depart_photon | Point de depart du photon
            (float,float,float) vecteur_photon | Vecteur associe au photon
            (float,float) coordonne_cylindre | Coordonnes du centre du cylindre
            float rayon_cylindre | Rayon du cylindre
            float hauteur_cylindre | Hauteur du cylindre


   @output: (bool,float) | Nous indique si intersection et la distance

   @utilite: Fonction utilise pour le calcul des radiations
   """

    px, py, pz = point_depart_photon
    xd, yd, zd = vecteur_photon
    xc, yc = coordonne_cylindre

    a = (xd**2)+(yd**2)
    b = (2*(px-xc)*xd)+(2*(py-yc)*yd)
    c = ((px-xc)**2) + ((py-yc)**2) - (rayon_cylindre**2)
    racines = np.roots([a, b, c])
    racine = [x for x in racines if x >= 0 and np.isreal(x)]
    if len(racine) == 0:
        return (False, -1)
    if len(racine) == 1:
        return (False, -1)
    dist = (lambda x, y, z: np.sqrt((x**2)+(y**2)+(z**2)))
    if len(racine) == 2:
        z1 = pz+racine[0]*zd
        z2 = pz+racine[1]*zd
        t1 = racine[0]
        t2 = racine[1]
        if (0 < z1 and z1 < hauteur_cylindre):
            if (0 < z2 and z2 < hauteur_cylindre):
                mm = min(t1, t2)
                return (True, dist((px+mm*xd-px), (py+mm*yd-py), (pz+mm*zd-pz)))
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
                    return (False, -1)
                if z1 > hauteur_cylindre and z2 > hauteur_cylindre:
                    return (False, -1)
                pass 
        mm = min(t1, t2)
        return (True, dist((px+mm*xd-px), (py+mm*yd-py), (pz+mm*zd-pz)))


def test_intersection():
    """
    @utilite:
    Cette fonction a pour seul but de verifier si la fonction intersection
    fonctionne correctement

    """

    assert intersection((1,0,1),(1,0,0),(3,0),1,3)[0] == True, "Problème 1"
    assert intersection((1,0,5),(1,0,0),(3,0),1,3)[0] == False, "Problème 2"
    assert intersection((1,0,5),(1,0,-2),(3,0),1,3)[0] == True, "Problème 3"
    assert intersection((1,0,5),(1,-5,-2),(3,0),1,3)[0] == False, "Problème 4"
    assert intersection((1,0,1),(1,0,2),(3,0),1,4)[0] == True, "Problème 5"
    assert intersection((1,0,1),(1,0,0),(3,0),1,4)[1] == 1,"Problème 6"
    assert intersection((1,0,1),(1,1,0),(3,0),1,4)[0] == False,"Problème 7"
    assert intersection((0,1,0.0000001),(1,0,0),(1.5,1.5),0.5,2)[0] == True, "Problème 8 Bord extreme"
    assert intersection((0,1.5,0.0000001),(1,0,0),(1.5,1.5),0.5,2)[0] == True, "Problème 9 Milieu"
    assert intersection((0,1.5,1),(1,0,1+0.0000001),(1.5,1.5),0.5,2)[0] == False, "Problème 10 Bord haut milieu tir"
    assert intersection((0,1.5,1),(1,0,1-0.0000001),(1.5,1.5),0.5,2)[0] == True, "Problème 10 Bord haut milieu tir bis"
    assert intersection((0,1.5,1),(1,0,1-0.0000001),(1.5,1.5),0.5,2)[0] == True, "Problème 10 Bord haut milieu tir bis"
    assert intersection((0,1,0.5+0.0000001),(1,0,1),(1.5,1.5),0.5,2)[0] == False, "Problème 10 Bord haut milieu tir extreme"
    assert intersection((0,1,0.5-0.0000001),(1,0,1),(1.5,1.5),0.5,2)[0] == True, "Problème 10 Bord haut milieu tir extreme bis"
    assert intersection((0.5,1.5,3.5-0.0000001),(1,0,-1),(1.5,1.5),0.5,2)[0] == True,"Problème 11 Couvercle haut extreme"
    assert intersection((0.5,1.5,3.5+0.0000001),(1,0,-1),(1.5,1.5),0.5,2)[0] == False,"Problème 11 Couvercle haut extreme bis"
 

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


def test_grille():
    obj_grille = Grille2D(taille=Taille_grille)  # On créer l'objet grille
    obj_grille.grille = [[Case() for j in range(Taille_grille)] for i in range(Taille_grille)]  # On rempli l'objet grille avec des cases végétales
    for k in range(1):
        obj_grille.grille[25+k][25] = Case(feu=True, isdry=True, temperature=500)  # On place une source de feu
    dt = 0.1
    obj_grille.xupdate(int((1/dt))*20, timage=5)  # On utilise une fonction qui met à jour la grille et nous stocke une image de la situation regulierement
