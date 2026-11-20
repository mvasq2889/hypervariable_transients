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

for ra, dec in list(itertools.product(ra_s, dec_s)):
    RA_MIN = float(ra)
    RA_MAX = float(RA_MIN + 10)
    DEC_MIN = float(dec)
    DEC_MAX = float(DEC_MIN + 10)
        
    loci_ids = query_by_chunk(RA_MIN, RA_MAX, DEC_MIN, DEC_MAX)
    if isinstance(loci_ids, list):
        print(f'({RA_MIN} - {RA_MAX}, {DEC_MIN} - {DEC_MAX}) : {len(loci_ids)}')
        count += len(loci_ids)
        for l_id in loci_ids:
            locus = search.get_by_id(l_id)
            oid_locus = DevKitLocus.model_validate(locus.to_devkit())
            oid_here = oid_get_info(oid_locus.alerts)
                
            if isinstance(oid_here, dict) == 1:
                oid_here['id'] = l_id
                    
                if isinstance(too_targets, int) == 1:
                    too_targets = pd.DataFrame(get_target_info(locus, oid_here), index=[0])
                else:
                    target = pd.DataFrame(get_target_info(locus, oid_here), index=[-1])
                    too_targets = pd.concat([too_targets, target]).reset_index(drop=True)      

if isinstance(too_targets, int) == 0:
    cwd = os.getcwd()
    too_targets.to_csv(f'{cwd}/query_data/too_targets_{int(Time.now().mjd)}.csv')
end_query = time.time()
print(f'uab ps 1.docxin {end_query - start_query} s ({(end_query - start_query)/3600} min)')
