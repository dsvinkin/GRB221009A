import sys
import os

from astropy.coordinates import EarthLocation
from astropy.time import Time

proxy = 'http://www-proxy:3128'

os.environ['http_proxy'] = proxy 
os.environ['HTTP_PROXY'] = proxy
os.environ['https_proxy'] = proxy
os.environ['HTTPS_PROXY'] = proxy

#EarthLocation._get_site_registry(force_download=True)

t = Time('2022-10-09T13:17:00', format='isot', scale='utc')

#lst_sites = EarthLocation.get_site_names()
#print(lst_sites)
#sys.exit()

#site = EarthLocation.of_site('ICECUBE')  
#print(site.geodetic) 
#print(site.info)
#print(site.get_gcrs(t))


lst_sites = 'IceCube LHAASO Carpet-3 Baikal-GVD KM3NeT'.split()
dic_sites = {
    'IceCube':    EarthLocation.of_site('ICECUBE'),
    'LHAASO':     EarthLocation.from_geodetic(100.1375, 29.358611, height=4410),  # 29.358611 100.1375
    'Carpet-3':   EarthLocation.from_geodetic(42.6667, 43.41667, height=0.0), #  43o25'N  42o40'E
    'Baikal-GVD': EarthLocation.from_geodetic(108.1650, 53.5587, height=0.0), # 53.5587 N, 108.1650 E
    'KM3NeT':     EarthLocation.of_site('KM3NeT ARCA'),
}

for instr in lst_sites:

    gcrs_ = dic_sites[instr].get_gcrs(t)
    print('{:12s} {:>8.3f} {:>8.3f} {:>8.1f}'.format(instr, gcrs_.ra.value, gcrs_.dec.value,  gcrs_.distance.value/1000))
    



