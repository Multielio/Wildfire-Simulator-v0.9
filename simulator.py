# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 13:02:55 2018

"""

from enum import Enum
import numpy as np
import random
import matplotlib.pyplot as plt
import os
from param import *
from radiation import *
from convection import *
from imager import *

"""
Creation des differents objets
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

class Extractor(object):
    def __init__(self,grille):
        self.states = []
        self.current = grille
        self.count = 0
    def add(self,Grille):
        self.states.append(Grille)
    def generate(self,dt,secondsToAdd,storeOneOf_States,liveProduce=True,storage="store/"):
        for i in range(int(secondsToAdd//dt)):
            self.current.update(dt, debug=False, figs=None)
            if i % storeOneOf_States ==0:
                self.add(self.current)
                if liveProduce:
                    fig1 = plt.figure()
                    image_situation(self.current,fig1)
                    image_temperature(self.current,fig1)
                    image_situation_3D(self.current,fig1)
                    path = storage+str(self.count)+'_u.png'
                    self.count +=1
                    if os.path.isfile(path):
                        os.remove(path)
                    fig1.savefig(path, format='png')

    

class Grille2D(object):

    """
    Notre grille est une matrice d'objets Case

    @input: int taille | Taille de la grille

    @output: Grille | Nouvelle grille

    @fonction: update(float) -> ()
                    Utilite: Avance le temps de dt
                    Limites: si dt trop grand non cohérent
    """
    def __init__(self,msg=True,para=Param()):
        self.taille = para.Taille_grille
        self.grille = []
        self.radiation_enreception = np.zeros(shape=(self.taille, self.taille))
        self.convection_enreception = np.zeros(shape=(self.taille, self.taille))
        self.time = 0
        self.rad = []
        self.fires=[]
        self.radgrids= {}
        self.pa= para
        self.loadrad()
        if msg:
            self.pa.printdata()
    
    def loadrad(self):
        for o in range(self.pa.Nombre_hauteur_differentes):
            self.rad.append((np.load(self.pa.Precal_path+str(o)+'_{}.npy'.format(self.pa.Nomrad))))

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
        self.hauteurflamme_matrix = np.zeros(shape=(self.taille, self.taille))
        self.longeurflamme_matrix = np.zeros(shape=(self.taille, self.taille))
        self.temperature_matrix = np.zeros(shape=(self.taille, self.taille))
        self.energie_matrix = np.zeros(shape=(self.taille, self.taille))
        
        for i in range(self.taille):
            for j in range(self.taille):
                caseij = self.grille[i][j]
                if not(caseij.hascombustible):
                    continue
                
                #######################################################
                self.temperature_matrix[i, j] = caseij.temperature
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
                    self.hauteurflamme_matrix[i, j] = hflamme
                    longeur_flamme = hflamme_no_vent
                    self.longeurflamme_matrix[i, j] = longeur_flamme
                    Sf = np.pi * self.pa.D*longeur_flamme
                    Pfprime = self.pa.Xr*Qpoint/Sf
                    self.energie_matrix[i, j] = self.pa.absorptivity*(Pfprime/Volume)*Sf
                    caseij.temperature = (Pfprime/(self.pa.sigma_boltzmann*(1-np.exp(-self.pa.Kf*longeur_flamme))))**(1/4)
         
                    self.temperature_matrix[i, j] = caseij.temperature
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

        self.radiation_enreception,self.fires,self.radgrids = reception_matrix(self.pa,self.rad, self.hauteurflamme_matrix, self.energie_matrix,self.radiation_enreception,self.fires,self.radgrids)
        self.convection_enreception = reception_matrix_conv(self.longeurflamme_matrix, self.temperature_matrix)
  
        
            







