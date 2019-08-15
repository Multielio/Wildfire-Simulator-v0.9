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
from radiation import *
from param import *


def worker(part):
    parameters= Param()
    dteta =  parameters.dteta
    dh = parameters.dh
    nphoton = parameters.nphoton
    R = parameters.D/2
    h = parameters.hcom
    taille =parameters.Taille_grille
    premierehauteur = parameters.Init_h
    espacement = parameters.Espacement_H

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


parameters= Param()
if __name__ == '__main__':
    print("-------------------------------")
    print(" Précalcul des radiations ")
    print("-------------------------------")
    print("Données chargées depuis config.yml")
    print("-------------------------------")
    print("Heure début: {}".format(datetime.datetime.now()))
    print("-------------------------------")
    manager = multiprocessing.Manager()
    jobs = []
    y = datetime.datetime.now()

    with multiprocessing.Pool(processes=parameters.nbr_cores) as p:
        max_ = 50
        with tqdm(total=max_) as pbar:
            for i, _ in tqdm(enumerate(p.imap_unordered(worker, range(0, max_)))):
                print("\n")
                pbar.update()
                print("\n")

    print(datetime.datetime.now()-y)
    print("Terminé")
    input("Close2")
    
    





    
