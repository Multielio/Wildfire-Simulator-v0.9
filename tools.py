import numpy as np

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
 