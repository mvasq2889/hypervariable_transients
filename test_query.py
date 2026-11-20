#same as all sky query, but specific date and saves all oid that pass quality thresh

from filter_funcs import *
from antares_devkit.models import DevKitLocus
import itertools
import os
from astropy.time import Time

oid_sample = 0
too_targets = 0
ra_s = np.arange(0, 360, 10)
dec_s = np.arange(-90, 90, 10)
count = 0
start = 61361.5

for ra, dec in list(itertools.product(ra_s, dec_s)):
    RA_MIN = float(ra)
    RA_MAX = float(RA_MIN + 10)
    DEC_MIN = float(dec)
    DEC_MAX = float(DEC_MIN + 10)
    
    loci = query_by_chunk(RA_MIN, RA_MAX, DEC_MIN, DEC_MAX,
                          mjd_start = start, mjd_end=start+1)
    if isinstance(loci, list):
        print(f'({RA_MIN} - {RA_MAX}, {DEC_MIN} - {DEC_MAX}) : {len(loci_ids)}')
        count += len(loci)
        for locus in loci:
            oid_here = oid_get_info(locus.alerts)
            
            if isinstance(oid_here, dict) == 1:
                oid_here['id'] = locus.locus_id
                if isinstance(oid_sample, int) == 1:
                    oid_sample = pd.DataFrame(oid_here, index = [0])
                else:
                    oid_here = pd.DataFrame(oid_here, index = [-1])
                    oid_sample = pd.concat([oid_sample, oid_here]).reset_index(drop=True)

                if oid_here['delta_m'] >= 1:
                        if isinstance(too_targets, int) == 1:
                            too_targets = pd.DataFrame(get_target_info(locus, oid_here), index=[0])
                        else:
                            target = pd.DataFrame(get_target_info(locus, oid_here), index=[-1])
                            too_targets = pd.concat([too_targets, target]).reset_index(drop=True)
                   
cwd = os.getcwd()
oid_sample.to_csv(f'{cwd}/test_data/oid_sample_{start-0.5}.csv')
 
if isinstance(too_targets, int) == 0:
    too_targets.to_csv(f'{cwd}/test_data/too_targets_{start-0.5}.csv')
    
