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
    
    loci_ids = query_by_chunk(RA_MIN, RA_MAX, DEC_MIN, DEC_MAX,
                              mjd_start = start, mjd_end=start+1)
    print(f'({RA_MIN} - {RA_MAX}, {DEC_MIN} - {DEC_MAX}) : {len(loci_ids)}')
    count += len(loci_ids)
    for l_id in loci_ids:
        locus = search.get_by_id(l_id)
        oid_dv = locus.to_devkit()
        oid_locus = DevKitLocus.model_validate(oid_dv)
        
        oid_here = oid_get_info(oid_locus.alerts)
        
        if isinstance(oid_here, dict) == 1:
            oid_here['id'] = l_id
            
            if isinstance(oid_sample, int) == 1:
                oid_sample = pd.DataFrame(oid_here, index = [0])
            else:
                if oid_here['delta_m'] >= 1:
                    if isinstance(too_targets, int) == 1:
                        too_targets = pd.DataFrame(get_target_info(locus, oid_here), index=[0])
                    else:
                        target = pd.DataFrame(get_target_info(locus, oid_here), index=[-1])
                        too_targets = pd.concat([too_targets, target]).reset_index(drop=True)
                oid_here = pd.DataFrame(oid_here, index = [-1])
                oid_sample = pd.concat([oid_sample, oid_here]).reset_index(drop=True)
                

if isinstance(too_targets, int) == 1:
    print(f'no hypervariables, checked {count} loci')
else:
    cwd = os.getcwd()
    too_targets.to_csv(f'{cwd}/test_data/too_targets_{start-0.5}.csv')
    oid_sample.to_csv(f'{cwd}/test_data/oid_sample_{start-0.5}.csv')
##    print(f'{too_targets}\n\nchecked {count} loci\n\n{oid_sample}')

