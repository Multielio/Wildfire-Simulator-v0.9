from simulator import *



def test_grille():
    parameters = Param()
    obj_grille = Grille2D(msg=True,para=parameters)  # On créer l'objet grille
    obj_grille.grille = [[Case() for j in range(parameters.Taille_grille)] for i in range(parameters.Taille_grille)]  # On rempli l'objet grille avec des cases végétales
    for k in range(1):
        obj_grille.grille[25+k][25] = Case(feu=True, isdry=True, temperature=500)  # On place une source de feu
    extract = Extractor(obj_grille)
    extract.generate(0.1,10,2,True)
    

test_grille()