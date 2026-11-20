from filter_funcs import *
from antares_devkit.models import DevKitLocus
import itertools
import os
from astropy.time import Time

too_targets = 0
ra_s = np.arange(0, 360, 10)
dec_s = np.arange(-90, 90, 10)
count = 0
start_query = time.time()

#query by chunk, gets list of loci w/ alerts w/in last 24 hrs and mag<21
#check loci 1 by 1 for pass quality and variability thresholds
#create target list
for ra, dec in list(itertools.product(ra_s, dec_s)):
    RA_MIN = float(ra)
    RA_MAX = float(RA_MIN + 10)
    DEC_MIN = float(dec)
    DEC_MAX = float(DEC_MIN + 10)
        
    loci = query_by_chunk(RA_MIN, RA_MAX, DEC_MIN, DEC_MAX)
    if isinstance(loci, list):
        print(f'({RA_MIN} - {RA_MAX}, {DEC_MIN} - {DEC_MAX}) : {len(loci_ids)}')
        count += len(loci)
        for locus in loci:
            oid_here = oid_get_info(locus.alerts)
                
            if isinstance(oid_here, dict) == 1 and oid_here['delta_m'] >= 1:
                if isinstance(too_targets, int) == 1:
                    too_targets = pd.DataFrame(get_target_info(locus, oid_here), index=[0])
                else:
                    target = pd.DataFrame(get_target_info(locus, oid_here), index=[-1])
                    too_targets = pd.concat([too_targets, target]).reset_index(drop=True)      

if isinstance(too_targets, int) == 0:
    cwd = os.getcwd()
    too_targets.to_csv(f'{cwd}/query_data/too_targets_{int(Time.now().mjd)}.csv')
end_query = time.time()
print(f'searched {count} loci in {end_query - start_query} s ({(end_query - start_query)/3600} min)')
