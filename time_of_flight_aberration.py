import sys

import numpy as np

c = 299792.458 # km/s

from jplephem.spk import SPK

bsp = 'c:/work/IPN_triangulation/pyipn/triangulation/kernels/de405.bsp'

AE = 1.49597870e13 # cm
t_unit = 86400.0 #s

def JD(nYear, nMonth, nDay, fTime):
    """
    a, y, and m are integers
    """ 
    a = (14 - nMonth) // 12
    y = nYear + 4800 - a
    m = nMonth + 12 * a - 3
    nJD = (153 * m + 2) // 5 + y // 4 - y // 100 + y // 400

    fJD = nJD + nDay + 365.0 * y - 32045.0
    fJD += fTime / 86400.0 - 0.5

    return fJD

def MJD(nYear, nMonth, nDay, fTime):

    return JD(nYear, nMonth, nDay, fTime) - 2400000.5

def GetSolarEphemeris(JD):

    t = JD - 2451545.0  # time since epoch J2000.0
    g = 357.529 + 0.98560028 * t
    q = 280.459 + 0.98564736 * t

    g = np.deg2rad(g)

    # ecliptic longitude
    L = q + 1.915 * np.sin(g) + 0.020 * np.sin(2 * g) # degrees!

    e = 23.439 - 0.00000036 * t # degrees!

    L = np.deg2rad(L)
    e = np.deg2rad(e)

    RA = np.rad2deg(np.arctan2(np.cos(e) * np.sin(L), np.cos(L))) # degrees
    
    if (RA < 0):
        RA += 360

    Dec = np.rad2deg(np.arcsin(np.sin(e) * np.sin(L))) # degrees
    
    return RA, Dec


def GES_to_DEC(r_GES):
    """
    from GES to DEC
    """

    r = r_GES[0]
    alpha = np.deg2rad(r_GES[1])
    delta = np.deg2rad(r_GES[2])

    r_DEC = np.zeros(3)
    r_DEC[0] = r * np.cos(alpha) * np.cos(delta)
    r_DEC[1] = r * np.sin(alpha) * np.cos(delta)
    r_DEC[2] = r * np.sin(delta)

    return r_DEC

def DEC_to_GES(r_DEC):
    """
    from DEC to GES
    """
    r = np.linalg.norm(r_DEC)
    alpha = np.arctan2(r_DEC[1], r_DEC[0])

    if (alpha < 0):
        alpha += 2*np.pi
    
    delta = np.arcsin(r_DEC[2] / r)

    r_GES = np.zeros(3)
    r_GES[0] = r
    r_GES[1] = np.rad2deg(alpha)
    r_GES[2] = np.rad2deg(delta)

    return r_GES

def print_vector_G(arr):

    print(f'r: {arr[0]} RA: {arr[1]} DEC: {arr[2]}')

def print_vector_D(arr):

    print(f'x: {arr[0]} y: {arr[1]} z: {arr[2]}')

def Aberr_rel(r1_G, r2_G, V, dt):
    """
    Input:
    r1_G - position of s/c 1 (r [km], ra [deg], dec [deg])
    r2_G - position of s/c 2 (r [km], ra [deg], dec [deg])

    V - the Earth velosity relative to the Sun in units of light speed
    dt - GRB time delay

    Output:
    dr_G - s/c 1 to s/c 2 vector (annulus center in case of triangulation)
   
    Comment:
    fdt > 0 - light travel time from s/c1 to s/c2
    on s/c1 - the burst came first
    r1_G, r2_G - R, RA, Dec
    """

    #print_vector_G(r1_G)
    #print_vector_G(r2_G)

    # transforn s/c vectors from Sperical to Cartesian coordinates
    r1_D = GES_to_DEC(r1_G)
    r2_D = GES_to_DEC(r2_G)

    # convert to light seconds
    r1_D /= c
    r2_D /= c

    #print_vector_D(r1_D)
    #print_vector_D(r2_D)

    dr_D = r1_D - r2_D

    V_abs = np.linalg.norm(V)
    dr_D_abs = np.linalg.norm(dr_D)
    #print(f'{V_abs=}')
    #print(f'{dr_D_abs=}')

    # Calculate projection of V on dr_D
    # Vdr = V dot dr_D
    Vdr = np.dot(V, dr_D)
    #print(f'{Vdr=} {dt=}')

    r_corr = V * (Vdr - dt)
    r_corr_abs = np.linalg.norm(r_corr)
    #print(f'{r_corr_abs=}')

    # Aberration correction
    # see eq. 21 in http://adsabs.harvard.edu/abs/1981Ap%26SS..75..219B
    for i in range(3):
        dr_D[i] += V[i] * (Vdr - dt)

    dr_G = DEC_to_GES(dr_D)

    return dr_G, dr_D

def get_solar_system_data(fJD):
    """
    From https://rhodesmill.org/skyfield/planets.html
    0 -> 3    SOLAR SYSTEM BARYCENTER -> EARTH BARYCENTER
    3 -> 399  EARTH BARYCENTER -> EARTH
    0 -> 10   SOLAR SYSTEM BARYCENTER -> SUN

    compute_and_differentiate returns position in kilometers and velocity  kilometers per day

    see https://github.com/brandon-rhodes/python-jplephem/issues/19


    """

    kernel = SPK.open(bsp)

    #print(kernel.comments())

    rEarth_bc, vEarth_bc = kernel[0,3].compute_and_differentiate(fJD)
    #rEarth, vEarth = kernel[3,399].compute_and_differentiate(jd)
    rSun_bc, vSun_bc = kernel[0,10].compute_and_differentiate(fJD)

    # Check if the units are correct 
    #print('rEarth_bc={:.3e} vEarth_bc={:.3e}'.format(
    #    np.linalg.norm(rEarth_bc) * 1e5/AE, 
    #    np.linalg.norm(vEarth_bc) / t_unit))
    #sys.exit()
    

    return rEarth_bc, vEarth_bc / t_unit, rSun_bc, vSun_bc / t_unit

def correct_time_delay(r1_G, fT1, r2_G, fdT, nYear, nMonth, nDay):
    
    """
    Input:

    r1_G - position of s/c 1 in sperical coordinates
    fT1  - time zerof of s/c 1
    r2_G - position of s/c 2 in sperical coordinates
    fdT  - can be positive or negative
    
    Output:
 
    r12_G  - vector from s/c 1 to s/c 2 in sperical coordinates (r, alpha, delta)
    rSun_G
    vSun
    """

    fJD = JD(nYear, nMonth, nDay, fT1)

    # r in km, v in km/s
    rEarth_bc, vEarth_bc, rSun_bc, vSun_bc = get_solar_system_data(fJD)

    rSun_G = DEC_to_GES(rSun_bc - rEarth_bc)
    rSun_G[0] *= 1e5/AE
    rSun_G[1] = np.rad2deg(rSun_G[1])
    rSun_G[2] = np.rad2deg(rSun_G[2])

    #print("Sun coordinates (DE405):");
    #print("RA = {:.5f} deg Dec = {:+.5f} deg R = {:.3f} AU\n".format(rSun_G[1], rSun_G[2], rSun_G[0]))

    fSunRA, fSunDec = GetSolarEphemeris(fJD)

    #print("Sun coordinates (approximate):")
    #print("RA = {:.5f} deg Dec = {:+.5f} deg\n".format(fSunRA, fSunDec))


    # Calculate the Earth velocity relative to the solar system barycenter in the units of light speed
    
    

    #for i in range(3):
    #    pEarth[i] = -vSun_[i]/s_l #*AE/86400

    # скорость Солнца (абсолютное значение) км/сек
    #vSun = 1e-5 * s_l * (pEarth[0]*pEarth[0] + pEarth[1]*pEarth[1] + pEarth[2]*pEarth[2])**0.5 # km/sec
    vSun = np.linalg.norm(vEarth_bc)


    vEarth_c = vEarth_bc / c
    r12_G, dr12_D = Aberr_rel(r1_G, r2_G, vEarth_c, fdT)

    """
    if (fdT < 0):

        r12_G, dr12_D = Aberr_rel(r1_G, r2_G, vEarth_c, fdT)
    
    else:
    
        r12_G, dr12_D = Aberr_rel(r2_G, r1_G, vEarth_c, -fdT)
    """
        
    return r12_G, dr12_D

class SC:

    def __init__(self, *args, **kwargs):

        if len(args) == 3 and len(kwargs) ==0:
            self.date = args[0] # YYYYMMDD
            self.time = args[1] # s UTC
            self.pos = args[2]

        elif len(args) == 1 and len(kwargs) ==0:
            lst_ = args[0].split()
            self.date = lst_[0] # YYYYMMDD
            self.time = float(lst_[1]) # s UTC

            # s/c geocentric coordinates in sperical system
            # RA [deg], Dec[deg], R[km]
            self.pos = np.array( [float(lst_[4]), float(lst_[2]), float(lst_[3])]) 

    def __str__(self,):

        return "{:s} {:.3f} {:f} {:f} {:f}".format(self.date, self.time, self.pos[0], self.pos[1], self.pos[2])

class TimeOfFlight:

    def __init__(self, ra, dec, sc_1, sc_2):

        self.sc_1 = sc_1
        self.sc_2 = sc_2
        self.ra = ra
        self.dec = dec
        self.dt, self.dt_sc_2, self.dt_aber = self.get_ToF(ra, dec, sc_1, sc_2)

    def __str__(self,):

        str_ ="""\
            Info                      sc1             sc2      GRB    
            RA (deg)       {:14.3f}  {:14.3f}  {:9.3f}  
            Dec (deg)      {:14.3f}  {:14.3f}  {:9.3f}  
            R (km)         {:14.3f}  {:14.3f}           
            T0 (s)         {:14.3f}  {:14.3f}           

            R (lt-s)       {:14.3f}  {:14.3f}
            Angle (deg)    {:14.3f}  {:14.3f}
            Prop. time (s) {:14.3f}  {:14.3f}  {:9.3f}
            sc2 time (s)   {:14s}  {:14s}  {:9.3f}

        """.format(self.sc_1.pos[1], self.sc_2.pos[1], self.ra,
                   self.sc_1.pos[2], self.sc_2.pos[2], self.dec,
                   self.sc_1.pos[0], self.sc_2.pos[0], 
                   self.sc_1.time,   self.sc_2.time,
                   self.dist_1,      self.dist_2,
                   self.angle_1,     self.angle_2,
                   self.prop_1,      self.prop_2, self.dt, 
                   "",  "", self.dt_sc_2)
        str_ = "\n".join([s.lstrip() for s in str_.split('\n')])
 
        return str_
    
    
    def get_ToF(self, ra, dec, sc_1, sc_2):
    
        r_src = GES_to_DEC(np.array([1.0, ra, dec]))
        r_sc_1 = GES_to_DEC(sc_1.pos)
        r_sc_2 = GES_to_DEC(sc_2.pos)

        r_sc_1_len = np.linalg.norm(r_sc_1)
        r_sc_2_len = np.linalg.norm(r_sc_2)

        self.prop_1 = np.dot(r_sc_1, r_src) / c
        self.prop_2 = np.dot(r_sc_2, r_src) / c

        if r_sc_1_len > 0.0:
            self.angle_1 = np.rad2deg( np.arccos(np.dot(r_sc_1, r_src)/r_sc_1_len) )
        else:
            self.angle_1 = 0.0

        if r_sc_2_len > 0.0:
            self.angle_2 = np.rad2deg( np.arccos(np.dot(r_sc_2, r_src)/r_sc_2_len) )
        else:
            self.angle_2 = 0.0

        self.dist_1 = r_sc_1_len / c
        self.dist_2 = r_sc_2_len / c
    
        dt = np.dot((r_sc_1 - r_sc_2), r_src) / c


        nYear, nMonth, nDay = int(sc_1.date[:4]), int(sc_1.date[4:6]), int(sc_1.date[6:8])

        dr12_G, dr12_D = correct_time_delay(sc_1.pos, sc_1.time, sc_2.pos, dt, nYear, nMonth, nDay)

        #print(np.linalg.norm(dr12_D))
        dt_aber = np.dot(dr12_D, r_src)
        
     
        dt_sc = sc_1.time - sc_2.time 
    
        return dt, dt_sc + dt, dt_aber

def main():

    """HEBS is mouted at SATech 01 NORAD ID is mentioned at 
     https://www.orbitalfocus.uk/Diaries/Launches/GeoSS/ss-PRC.php
    """

    ra, dec = 288.2643, 19.7712  # deg GRB position

    lst_instr = """Fermi  Swift  AGILE  SATech-01-HEBS  Insight-HXMT-HE GECAM-B
        INTEGRAL-SPI-ACS L2_SRG-ARTXC  Wind-Konus 
        Mars-Odyssey-HEND  BepiColombo-MGNS  Solar-Orbiter-STIX Voyager_1 Voyager_2
        IceCube  LHAASO Carpet-3 Baikal-GVD KM3NeT
     """.split()

    #lst_instr = ['Voyager_1', 'Voyager_2']

    dic_pos = {
        'Fermi':             SC('20221009 47819.988 298.61316  -4.57939  6910.4'),
        'Swift':             SC('20221009 47819.988  59.23429  -3.62475  6924.3'),
        'AGILE':             SC('20221009 47819.988 279.70061  -1.84728  6825.5'),
        'GECAM-B':           SC('20221009 47819.988 59.38016   -2.25436  6960.7'),
        'SATech-01-HEBS':    SC('20221009 47820.050 342.07152  61.19406  6863.5'),
        'Insight-HXMT-HE':   SC('20221009 47820.050 272.88385 -10.92975  6911.2'),
        'INTEGRAL-SPI-ACS':  SC('20221009 47819.988 199.505    69.422  137648.0'),
        'L2_SRG-ARTXC':      SC('20221009 47995.000  14.65504  6.24952 1.5106176314E+06'),
        'Wind-Konus':        SC('20221009 47821.648 181.71703 -3.62390 1294914.5'),
        'Mars-Odyssey-HEND': SC('20221009 47819.988  81.7184  22.8347 1.099847e+08'),
        'BepiColombo-MGNS':  SC('20221009 47819.988 219.4658 -19.0687 1.358263e+08'),
        'Solar-Orbiter-STIX':SC('20221009 47819.988 203.5964 -10.8673 1.842560e+08'),
        'Voyager_1':         SC('20221008 65712.487 258.28582  12.14229  2.3656606778E+10'),
        'Voyager_2':         SC('20221008 35774.758 301.43064 -58.89899  1.9658309639E+10'),
        'IceCube':   SC('20221009 47819.988 180.757 -89.875 6357.6'), 
        'LHAASO' :   SC('20221009 47819.988 317.286  29.104 6377.4'),
        'Carpet-3':  SC('20221009 47819.988 259.883  43.248 6368.1'),
        'Baikal-GVD':SC('20221009 47819.988 325.363  53.272 6364.3'),
        'KM3NeT':    SC('20221009 47819.988 233.273  36.160 6367.2'),
    } 

    pos_Earth_center = SC('20221009 47819.988 0.0 0.0 0.0')
    print(f'{"Instrument":18s} {"ToF":>11s} {"ToF aberr.":>11s} {"aberr. corr.":>12s}')
    print(f'{" ":18s} {"(s)":>11s} {"(s)":>11s} {"(s)":>12s}')

    for instr in lst_instr:

        tof = TimeOfFlight(ra, dec, pos_Earth_center, dic_pos[instr])
        print(f'{instr:18s} {tof.dt:>11.3f} {tof.dt_aber:>11.3f} {tof.dt_aber-tof.dt:>12.3e}')

    #with open('ToF_SRG_Wind.txt', 'w') as f:
    #    f.write(str(tof))


main()