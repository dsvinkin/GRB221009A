import numpy as np

c = 299792.458 # km/s

class SC:

    def __init__(self, *args, **kwargs):
        if len(args) == 3 and len(kwargs) ==0:
            self.date = date # YYYYMMDD
            self.time = time # s UTC
            self.pos = arr_pos

        elif len(args) == 1 and len(kwargs) ==0:
            lst_ = args[0].split()
            self.date = lst_[0] # YYYYMMDD
            self.time = float(lst_[1]) # s UTC
            self.pos = [ float(x) for x in lst_[2:5] ]

    def __str__(self,):

        return "{:s} {:.3f} {:f} {:f} {:f}".format(self.date, self.time, self.pos[0], self.pos[1], self.pos[2])

class TimeOfFlight:

    def __init__(self, ra, dec, sc_1, sc_2):

        self.sc_1 = sc_1
        self.sc_2 = sc_2
        self.ra = ra
        self.dec = dec
        self.dt, self.dt_sc_2 = self.get_ToF(ra, dec, sc_1, sc_2)

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

        """.format(self.sc_1.pos[0], self.sc_2.pos[0], self.ra,
                   self.sc_1.pos[1], self.sc_2.pos[1], self.dec,
                   self.sc_1.pos[2], self.sc_2.pos[2], 
                   self.sc_1.time,   self.sc_2.time,
                   self.dist_1,      self.dist_2,
                   self.angle_1,     self.angle_2,
                   self.prop_1,      self.prop_2, self.dt, 
                   "",  "", self.dt_sc_2)
        str_ = "\n".join([s.lstrip() for s in str_.split('\n')])
 
        return str_
    
    def cartesian_to_spherical(self, r):
        
        Rabs = np.linalg.norm(r)
        fDec = arcsin(r[2] / Rabs)
        fRA = arctan2(r[1], r[0])
        
        if (fRA < 0): 
            fRA += 2 * np.pi
    
        return fRA, fDec
    
    def spherical_to_cartesian(self, r):
        
        r_cart = np.zeros(3)
       
        r_cart[0] = r[2] * np.cos(np.deg2rad(r[1])) * np.cos(np.deg2rad(r[0]))
        r_cart[1] = r[2] * np.cos(np.deg2rad(r[1])) * np.sin(np.deg2rad(r[0]))
        r_cart[2] = r[2] * np.sin(np.deg2rad(r[1]))

        return r_cart
    
    def get_ToF(self, ra, dec, sc_1, sc_2):
    
        r_src = self.spherical_to_cartesian(np.array([ra, dec, 1.0]))
        r_sc_1 = self.spherical_to_cartesian(sc_1.pos)
        r_sc_2 = self.spherical_to_cartesian(sc_2.pos)

        r_sc_1_len = np.linalg.norm(r_sc_1)
        r_sc_2_len = np.linalg.norm(r_sc_2)

        self.prop_1 = np.dot(r_sc_1, r_src) / c
        self.prop_2 = np.dot(r_sc_2, r_src) / c

        self.angle_1 = np.rad2deg( np.arccos(np.dot(r_sc_1, r_src)/r_sc_1_len) )
        self.angle_2 = np.rad2deg( np.arccos(np.dot(r_sc_2, r_src)/r_sc_2_len) )

        self.dist_1 = r_sc_1_len / c
        self.dist_2 = r_sc_2_len / c
    
        dt = np.dot((r_sc_1 - r_sc_2), r_src) / c
     
        dt_sc = sc_1.time - sc_2.time 
    
        return dt, dt_sc + dt

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
        'Voyager_1':         SC('20221009 47819.988 258.28835 12.13962  2.3659203499E+10'),
        'Voyager_2':         SC('20221008 35774.758 301.43064 -58.89899  1.9658309639E+10'),
        'IceCube':   SC('20221009 47819.988 180.757 -89.875 6357.6'), 
        'LHAASO' :   SC('20221009 47819.988 317.286  29.104 6377.4'),
        'Carpet-3':  SC('20221009 47819.988 259.883  43.248 6368.1'),
        'Baikal-GVD':SC('20221009 47819.988 325.363  53.272 6364.3'),
        'KM3NeT':    SC('20221009 47819.988 233.273  36.160 6367.2'),
    } 

    pos_Earth_center = SC('20221009 47819.988 0.0 0.0 0.0')
    
    for instr in lst_instr:

        tof = TimeOfFlight(ra, dec, pos_Earth_center, dic_pos[instr])
        print(f'{instr:18s}: {tof.dt:>11.3f} seconds')

    #with open('ToF_SRG_Wind.txt', 'w') as f:
    #    f.write(str(tof))


main()