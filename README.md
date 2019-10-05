# Wildfire-Simulator

![alt text](https://github.com/Multielio/Wildfire-Simulator/blob/master/img.png)

Simulator based on the real physical phenomenons acting in a wildfire, deals with wind, topography (with few changes), heat capacities..an more.
It took me 1 year to build it.


**It is mainly based on the following papers**
- Mohamed Drissi / Modeling the spreading of large-scale wildland fires 
- J.K.Adou,  A.D.V.Broub, B.Porterie / Modeling wildland fire propagation using a semi-physical network model / Case Studies in Fire Safety 4 (2015) 11–18 

**Main assumptions**

- No conduction
- No screening effect
- Flames are cylinders

**Required libraries**
- Matplolib
- Numpy (pip install numpy==1.16.2)
- Scipy
- tqdm


**How it works**
- Fill config.yml with your parameters
- Launch precalculus.py 
- Launch test.py
- The simulator (simulator.py) loads the precalculated data generated by multitaskingallhateurs.py
- The extractor object supervise the simulation by updating the grille2D object and puts the images of the simulation into the folder storage/ (you can tweak the extractor to get more images)

- At each update of Grille2D object {
1. The simulator makes a energetic balance for each cell 
2. The simulator determines the next state of each cell based on the energetic balance.
3. The simulator prepare data for the next update of the simulator.
}

**Quickstart**
- Launch test.py (it will use the default config.yml parameters and the default precalculus data)
- Wait 
- Images generated will be in the "store/" folder.

