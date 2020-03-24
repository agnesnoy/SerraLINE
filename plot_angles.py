import matplotlib.pyplot as plt
import numpy as np

#Colors
c_red = '#BC0000'
c_green = '#10A100'
c_blue = '#1000FF'
c_purple = '#C200CD'
c_cyan = '#00BBEF'
c_yellow = '#D6D600'
label_s = 7 #Label size. For ticks
legend_s = 5 #legend size
line_w = 1 #linewidth
mark_s = 2 #markersize

colours=[c_red, c_green, c_blue, c_purple, c_cyan, c_yellow]
############################################################################
#plot

l = 0
for x in range(1, 6, 1):

    if l > len(colours)-1 :
        l = 0

    macolor = colours[l]
    l = l + 1

    path = "subfragment_{}.out".format(x)
    dat = np.loadtxt(path)

    malabel = "step {}".format(x) #length label

    plt.plot( dat[:,0], dat[:,1], color=macolor, label=malabel) 
    plt.fill_between( dat[:,0], dat[:,1] - dat[:,2],  dat[:,1] + dat[:,2],
    alpha=0.5, edgecolor=macolor, facecolor=macolor, linewidth=0)


##############################################################################

plt.grid(True)
plt.xlabel('Angle position along the DNA')
plt.ylabel('Bending angle (deg)')
plt.legend()

plt.show()

