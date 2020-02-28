#This script takes the boss catalog and filters
#it by ra and dec
#ra:  0-45 and  317-360F
#dec -1.25 - 1.25
#then compares this subset to the FIRST catalog (with converted ra and dec values from hms to deg)
#and rejects objects(galaxies) that are within 1 arcminute (1') of a FIRST object
#This final catalog of boss galaxies is then saved.

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

print('Opening catalogs')
#open FIRST catalog:
first = pd.read_csv('../catalogs/south_clean_radec.csv')

#open boss catalog:
boss = pd.read_csv('../catalogs/boss_gal_3.csv')

print('Filtering BOSS')
#filter ra and dec
restrict_ra = boss.loc[(boss['ra'] <=45) | (boss['ra'] >=317)]
overlap = restrict_ra.loc[(restrict_ra['dec'] <= 1.25) & (restrict_ra['dec'] >= -1.25)]

filter_boss= Table.from_pandas(overlap)
filter_boss = filter_boss.to_pandas()

print('Cleaning BOSS requiring at least 1 arcmin separation from FIRST')

b_coord = SkyCoord(ra=filter_boss['ra']*u.degree, dec=filter_boss['dec']*u.degree)
f_coord = SkyCoord(ra=first['ra']*u.degree, dec=first['dec']*u.degree)

idx, d2d, d3d = b_coord.match_to_catalog_sky(f_coord)
b_mask = d2d < 1*u.arcminute

indx_b = np.where(b_mask)[0]

print('Cleaning complete!')
#drop bad indices from filtered boss catalog
clean_boss = filter_boss.drop(indx_b)

#save cleaned catalog as a new file
clean_boss.to_csv('cleaned_boss.csv')
print('Clean catalog saved.')
