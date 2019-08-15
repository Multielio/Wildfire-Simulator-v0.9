import numpy as np

class Param(object):
    def  __init__(self):
        self.Taille_grille = 50
        self.Nomrad = "50x50"
        self.Init_h = 0.51
        self.Espacement_H = 0.2
        self.Nombre_hauteur_differentes = 50
        self.Precal_path = "ray051_02_2/" 

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
        
    
    def printdata(self):
        print("Vchar: {}".format(self.Vchar))
        print("Kf: {}".format(self.Kf))
        print("Xr: {}".format(self.Xr))
        print("Tambiant: {} Kelvin".format(self.Tambiant))
        print("Tpyro: {} Kelvin".format(self.Tpyro))
        print("sigma_k: {} m^-1".format(self.sigma_k))
        print("alpha_k: {}".format(self.alpha_k))
        print("delta: {}".format(self.delta))
        print("Efb: {}".format(self.Efb))
        print("U: {} m.s^-1".format(self.U))
        print("Teta_horizontal_vent: {}".format(self.Teta_horizontal_vent))
        print("D: {} m".format(self.D))
        print("Tc: {} s".format(self.Tc))
        print("mprimeprimeDFF: {} kg.m^-2".format(self.mprimeprimeDFF))
        print("Rho_wff: {}  kg.m^-3".format(self.Rho_wff))
        print("Rho_dff: {}  kg.m^-3".format(self.Rho_dff))
        print("Cdff: {} J.kg^-1 .K^-1".format(self.Cdff))
        print("Cpeau: {} J.kg^-1 .K^-1".format(self.Cpeau))
        print("FMC: {}".format(self.FMC))
        print("hauteur_combustible: {}".format(self.hcom))
        print("--------------------------------------------")

    def loadDefault(self):
        self.Taille_grille = 50
        self.Nomrad = "50x50"
        self.Init_h = 0.51
        self.Espacement_H = 0.2
        self.Nombre_hauteur_differentes = 50
        self.Precal_path = "ray051_02_2/" 

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

    