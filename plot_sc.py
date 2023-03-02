"""
https://gist.github.com/WetHat/1d6cd0f7309535311a539b42cccca89c
"""

import matplotlib.pyplot as plt

from matplotlib.text import Annotation
from matplotlib.patches import FancyArrowPatch

from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d.axes3d import Axes3D

from skspatial.objects import Sphere

import numpy as np

from itertools import product, combinations


class Annotation3D(Annotation):

    def __init__(self, text, xyz, *args, **kwargs):
        super().__init__(text, xy=(0, 0), *args, **kwargs)
        self._xyz = xyz

    def draw(self, renderer):
        x2, y2, z2 = proj_transform(*self._xyz, self.axes.M)
        self.xy = (x2, y2)
        super().draw(renderer)

def _annotate3D(ax, text, xyz, *args, **kwargs):
    '''Add anotation `text` to an `Axes3d` instance.'''

    annotation = Annotation3D(text, xyz, *args, **kwargs)
    ax.add_artist(annotation)

setattr(Axes3D, 'annotate3D', _annotate3D)

class Arrow3D(FancyArrowPatch):

    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)

    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)
        
    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

        return np.min(zs) 

def _arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    '''Add an 3d arrow to an `Axes3D` instance.'''

    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)


setattr(Axes3D, 'arrow3D', _arrow3D)

def eq_to_cart(arr_eq, normed=False):

    u, v = np.deg2rad(arr_eq[0]), np.deg2rad(90-arr_eq[1])

    r = arr_eq[2]
    if normed:
        r = 1 

    x = r*np.cos(u)*np.sin(v)
    y = r*np.sin(u)*np.sin(v)
    z = r*np.cos(v)

    return x, y, z



fig = plt.figure()
ax = fig.add_subplot(projection='3d')

#ax.set_aspect("equal")

r_earth = 6356.8/1000


# draw sphere
"""
u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:20j]
x = r_earth*np.cos(u)*np.sin(v)
y = r_earth*np.sin(u)*np.sin(v)
z = r_earth*np.cos(v)
ax.plot_wireframe(x, y, z, color="b", linewidths=0.5)
"""

sphere = Sphere([0, 0, 0], r_earth)
sphere.plot_3d(ax, alpha=0.2)

"""
dic_sc = {
    'Fermi': [58.14010, -23.94314, 6913.6/1000],
    'Swift': [15.24823, -20.53904, 6928.4/1000],
    'Mars-Odyssey': [313.343, -18.819,  201618152.0/1000],
    'Wind': [2.14702, 1.94708, 1428868.3/1000],
    'INTEGRAL':[222.735, 52.717, 132660.2/1000]
}
"""

lst_instr = """Fermi Swift AGILE Astrosat GECAM-B  SATech-01-HEBS  Insight-HXMT 
    INTEGRAL Spektr-RG  Wind 
    Mars-Odyssey  BepiColombo  Solar-Orbiter Voyager-1 Voyager-2""".split()

lst_instr_distant = """INTEGRAL  Spektr-RG  Wind  Mars-Odyssey  
    BepiColombo  Solar-Orbiter  Voyager-1 Voyager-2""".split()

dic_pos = {
    'Fermi':         '20221009 47819.988 298.61316  -4.57939  6910.4',
    'Swift':         '20221009 47819.988  59.23429  -3.62475  6924.3',
    'AGILE':         '20221009 47819.988 279.70061  -1.84728  6825.5',
    'Astrosat':      '20221009 47819.988 132.78290   4.48190  7012.9',
    'SATech-01-HEBS':'20221009 47820.050 342.07152  61.19406  6863.5',
    'Insight-HXMT':  '20221009 47820.050 272.88385 -10.92975  6911.2',
    'GECAM-B':       '20221009 47819.988  59.38016  -2.25436  6960.7',
    'INTEGRAL':      '20221009 47819.988 199.505    69.422  137648.0',
    'Spektr-RG':     '20221009 47995.000  14.65504   6.24952  1.5106176314E+06',
    'Wind':          '20221009 47821.648 181.71703  -3.62390   1294914.5',
    'Mars-Odyssey':  '20221009 47819.988  81.7184   22.8347   1.099847e+08',
    'BepiColombo':   '20221009 47819.988 219.4658  -19.0687   1.358263e+08',
    'Solar-Orbiter': '20221009 47819.988 203.5964  -10.8673   1.842560e+08',
    'Voyager-1':     '20221009 47819.988 258.28835  12.13962  2.3659203499E+10',
    'Voyager-2':     '20221008 35774.758 301.43064 -58.89899  1.9658309639E+10',
} 

# RA, Dec for 2022-10-09T13:17:00
dic_ground_instr = {
    'IceCube':   [180.757,-89.875, 6357.6/1000], 
    'LHAASO' :   [317.286, 29.104, 6377.4/1000],
    'Carpet-3':  [259.883, 43.248, 6368.1/1000],
    'Baikal-GVD':[325.363, 53.272, 6364.3/1000],
    'KM3NeT':    [233.273, 36.160, 6367.2/1000],
}

dic_sc = {}
for instr in lst_instr:

    lst_ = dic_pos[instr].split()

    dic_sc[instr] = [float(lst_[2]), float(lst_[3]), float(lst_[4])/1000]

# plot LEO instruments
for name in dic_sc.keys():

    if name in set(lst_instr_distant):
        continue

    x,y,z = eq_to_cart(dic_sc[name])
    ax.scatter([x], [y], [z], color="b", s=50)
    ax.annotate3D(name, (x, y, z), xytext=(3, 3), textcoords='offset points')

# plot ground instruments
for name in dic_ground_instr.keys():

    x,y,z = eq_to_cart(dic_ground_instr[name])
    ax.scatter([x], [y], [z], color="r", s=50)
    ax.annotate3D(name, (x, y, z), xytext=(3, 3), textcoords='offset points')

# plot distant SC
for name in lst_instr_distant:

    x,y,z = eq_to_cart(dic_sc[name], normed=True)

    f = 9
    ax.arrow3D(0,0,0,
               x*f,y*f,z*f,
               mutation_scale=20,
               arrowstyle="-|>",
               linestyle='dashed')
    
    ax.annotate3D('{:s} at {:.1f} x $10^6$ km'.format(name, dic_sc[name][2]/1e3), 
        (x*f, y*f, z*f), xytext=(3, 3), textcoords='offset points')

# plot GRB direction
dic_grb ={'GRB221009A': [288.2643, 19.7712, 1e+10]}

x,y,z = eq_to_cart(dic_grb['GRB221009A'], normed=True)

f = 9
ax.arrow3D(0,0,0,
           x*f,y*f,z*f,
           mutation_scale=20,
           arrowstyle="-|>",
           linestyle='solid', color='red')

ax.annotate3D('GRB 221009A', 
    (x*f, y*f, z*f), xytext=(3, 3), textcoords='offset points')

ax.set_xlabel('X (1000 km)')
ax.set_ylabel('Y (1000 km)')
ax.set_zlabel('Z (1000 km)')

ax.set_xlim([-8, 8])
ax.set_ylim([-8, 8])
ax.set_zlim([-8, 8])

ax.set_title('Spacecraft position at 2022-10-09 13:17:00')

fig.tight_layout()

plt.show()