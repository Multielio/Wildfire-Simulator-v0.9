
import matplotlib.pyplot as plt
from matplotlib.colors import colorConverter
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.patches as mpatches
import numpy as np
import datetime

def image_temperature(grille2D, fig1):
        def f(x):
            return np.log10(x)

        f = np.vectorize(f)
        ax = fig1.add_subplot(2, 2, 4)
        ma = ax.matshow(f(grille2D.temperature_matrix), cmap=plt.get_cmap("afmhot"))
        ax.xaxis.tick_bottom()
        ax.set_title('Image thermique')
        cl= fig1.colorbar(ma, ticks=[-1,0,1], orientation='horizontal')
        cl.set_label("Température")

def image_situation_3D(grille2D, fig1):
    ax1 = fig1.add_subplot(222, projection="3d")
    xpos = []
    ypos = []
    zpos = []
    dz = []
    colorl = []
    width = depth = 2
    for i in range(grille2D.taille):
        for j in range(grille2D.taille):
            case = grille2D.grille[i][j]
            radius = grille2D.pa.D/2
            height = grille2D.hauteurflamme_matrix[i, j]
            if height == 0:
                height = grille2D.pa.hcom
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

def image_situation(grille2D, fig1):
    print("------------------------------")
    print("Temps : {} ".format(str(datetime.timedelta(seconds=grille2D.time))))
    print((grille2D.grille[23][6]).temperature)
    z = []
    for i in range(grille2D.taille):
        o = []
        for j in range(grille2D.taille):
            case = grille2D.grille[i][j]
            if case.feu:
                o.append(colorConverter.to_rgb("red"))
            elif case.mass_dry == 0:
                o.append(colorConverter.to_rgb("black"))
            elif case.temperature == 373:
                o.append(colorConverter.to_rgb("orange"))
            elif case.temperature > 373:
                o.append(colorConverter.to_rgb("brown"))
            elif (case.temperature < 373) and (case.temperature > grille2D.pa.Tambiant):
                x = 1-(case.temperature-grille2D.pa.Tambiant)/(373-grille2D.pa.Tambiant)
                c = ((1-x)*255, 80*x +(1-x)*140, 0)
                o.append(c)
            else:
                o.append(colorConverter.to_rgb(case.ctype.couleur))
        z.append(o) 
    ax1 = fig1.add_subplot(1, 2, 1)
    
    im= ax1.imshow(z)
    cmap = {1:colorConverter.to_rgb("red"),2:colorConverter.to_rgb("brown"),3:colorConverter.to_rgb("orange"),4:colorConverter.to_rgb("yellow"),5:colorConverter.to_rgb("black")}
    labels = {1:'Flamme',2:'Combustible 373<T<T_pyro',3:'Combustible T=373K ',4:"Combustible T_amb<T<373",5:'Charbon'}  
    patches =[mpatches.Patch(color=cmap[i],label=labels[i]) for i in cmap]
    ax1.legend(handles=patches,loc=3, bbox_to_anchor=(0, -0.25, 0.30, 0.30) )
    ax1.set_title('Image situation')
        
