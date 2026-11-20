from antares_client import search
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import pandas as pd
import math
import time
import warnings
import itertools
import pdb
from astropy.io import fits

#LSST
#quality check
def lsst_quality(alert):
    quality = True
    quality &= alert.properties['lsst_diaSource_reliability'] >= 0.6
    quality &= alert.properties['lsst_diaSource_psfFlux_flag'] == False
    quality &= alert.properties['lsst_diaSource_centroid_flag'] == False
    quality &= alert.properties['lsst_diaSource_shape_flag'] == False
    quality &= alert.properties['lsst_diaSource_isDipole'] == False
    quality &= alert.properties['lsst_diaSource_pixelFlags_saturated'] == False
    quality &= alert.properties['lsst_diaSource_pixelFlags_edge'] == False
    quality &= alert.properties['lsst_diaSource_pixelFlags_cr'] == False
    quality &= alert.properties['lsst_diaSource_pixelFlags_streak'] == False
    quality &= alert.properties['lsst_diaSource_scienceFlux'] >= 0
    quality &= alert.properties['lsst_diaSource_templateFlux'] >= 0

    if 'lsst_diaSource_ssObjectId' in alert.properties:
        quality &= alert.properties['lsst_diaSource_ssObjectId'] in (0, None, "0")
    return quality

#extract data
def lsst_alert(alert):
    obs_date = float(alert.properties['lsst_diaSource_midpointMjdTai'])
    oid_mag = float(-2.5*np.log10(alert.properties['lsst_diaSource_scienceFlux']) + 31.4)
    oid_mag_err = float(-2.5*alert.properties['lsst_diaSource_scienceFluxErr']/(alert.properties['lsst_diaSource_scienceFlux']*np.log(10)))
    temp_mag = float(-2.5*np.log10(alert.properties['lsst_diaSource_templateFlux']) + 31.4)
    temp_mag_err = float(-2.5*alert.properties['lsst_diaSource_templateFluxErr']/(alert.properties['lsst_diaSource_templateFlux']*np.log(10)))
    diff_mag = float(temp_mag - oid_mag)
    diff_mag_err = math.sqrt(oid_mag_err**2 + temp_mag_err**2)
    
    bandpass = alert.properties['lsst_diaSource_band']
    band = 'lsst_' + str(bandpass)

    if alert.properties['lsst_diaSource_isNegative'] == True:
        pos_or_neg = -1
    else:
        pos_or_neg = 1
    
    data = {'mjd': obs_date, 'band': band, 'm_sci': oid_mag, 'm_sci_err': oid_mag_err, 'm_template': temp_mag, 'm_template_err': temp_mag_err, 
                'm_diff': None, 'm_diff_err': None, 'delta_m': diff_mag, 'delta_m_err': diff_mag_err, 'pos_or_neg': pos_or_neg, 'm_sci_flag': None}
       
    return data


#ZTF
#quality check
def ztf_quality(alert):
    quality = True
    quality &= alert.properties['ztf_rb'] >= 0.65
    quality &= alert.properties['ztf_nbad'] == 0
    quality &= alert.properties['ztf_fwhm'] <= 5
    quality &= alert.properties['ztf_elong'] <= 1.2
    quality &= abs(alert.properties['ztf_magdiff']) <= 1.
    quality &= (1.0857/alert.properties['ztf_sigmapsf']) < 8
    quality &= alert.properties['ztf_distnr'] <= 1.
    return quality

#extract data
#uses formula from 260609 6/15 markdown cell
def ztf_alert(alert):
    bands = ['x', 'g', 'R', 'i']
    band_id = alert.properties['ztf_fid']
    obs_date = float(alert.properties['ant_mjd'])
    temp_mag = alert.properties['ztf_magnr']
    temp_mag_err = alert.properties['ztf_sigmagnr']
    diff_img_mag = alert.properties['ztf_magpsf']
    diff_img_mag_err = alert.properties['ztf_sigmapsf']
    n = 10.**(0.4*(temp_mag - diff_img_mag))
    flag = 0
    
    if alert.properties['ztf_isdiffpos'] == 'f' or alert.properties['ztf_isdiffpos'] == 0:
        pos_or_neg = -1
        if 1-n <= 0.:
            diff_mag = 10**(-10)
            oid_mag_err = np.sqrt( 2**2 *temp_mag_err**2 + diff_img_mag_err**2 )
            flag = 1
        else:
            diff_mag = -2.5* np.log10(1-n)
            oid_mag_err = np.sqrt( (1+n/(1-n))**2 * temp_mag_err**2 + (-n/(1-n))**2 * diff_img_mag_err**2 )
    else:
        pos_or_neg = 1
        if 1+n <= 0.:
            diff_mag = 10**(-10)
            oid_mag_err = np.sqrt( 2**2 * temp_mag_err**2 + diff_img_mag_err**2 )
            flag = 1
        else:
            diff_mag = -2.5* np.log10(1+n)
            oid_mag_err = np.sqrt( (1-n/(1+n))**2 * temp_mag_err**2 + (n/(1+n))**2 * diff_img_mag_err**2 )
    oid_mag = diff_mag + temp_mag
    diff_mag_err = np.sqrt(oid_mag_err**2 + diff_img_mag_err**2)
    
    if oid_mag <= 20.5:
        band = 'ztf_' + str(bands[band_id])
        data = {'mjd': obs_date, 'band': band, 'm_sci': oid_mag, 'm_sci_err': oid_mag_err, 'm_template': temp_mag, 'm_template_err': temp_mag_err, 
                'm_diff': diff_img_mag, 'm_diff_err': diff_img_mag_err, 'delta_m': diff_mag, 'delta_m_err': diff_mag_err, 'pos_or_neg': pos_or_neg, 'm_sci_flag': flag}
        return data


#unpacks all alerts (that pass a quality test)
def oid_get_info(alerts):
    new_alert = alerts[-1]
    oid = None
    
    if new_alert.properties['ant_survey'] == 1:
        if ztf_quality(new_alert) == True and (1.0857/new_alert.properties['ztf_sigmapsf']) < 8:
            oid = ztf_alert(new_alert)
    elif new_alert.properties['ant_survey'] == 4:
        if lsst_quality(new_alert) == True:
            oid = lsst_alert(new_alert)
    return oid


#catalog cross-matching
#provided mostly from jessie
def catalog_check(locus):
    flags = {'known_star': False, 'known_agn': False, 'known_galaxy': False, 'gaia_likely_star': False, 
                            'gaia_astro_star': False, 'ps1_pointlike': False}
    for cat in locus.catalog_objects:
        if cat == 'asassn_variable_catalog' or cat == 'asassn_variable_Catalog_v2_20190802' or cat == 'sdss_stars':
            flags['known_star'] = True
        if cat == 'veron_agn_qso' or cat == 'milliquas':
            flags['known_agn'] = True
        if cat == '2mass_xsc' or cat == 'nyu_valueadded_gals' or cat == 'sdss_gals' or cat == 'RC3':
            flags['known_galaxy'] = True
        if cat == 'gaia_dr3_gaia_source':
            if isinstance(locus.catalog_objects[cat][0]['classprob_dsc_combmod_star'], float) == 1 and locus.catalog_objects[cat][0]['classprob_dsc_combmod_star'] > 0.5:
                flags['gaia_likely_star'] = True
            if isinstance(locus.catalog_objects[cat][0]['parallax'], float) == 1:
                gdr3_check1 = (locus.catalog_objects[cat][0]['parallax'] + 0.017)/locus.catalog_objects[cat][0]['parallax_error']
                alpha = locus.catalog_objects[cat][0]['pmra']/locus.catalog_objects[cat][0]['pmra_error']
                beta = locus.catalog_objects[cat][0]['pmdec']/locus.catalog_objects[cat][0]['pmdec_error']
                gamma = locus.catalog_objects[cat][0]['pmra_pmdec_corr']
                gdr3_check2 = np.sqrt((alpha**2 + beta**2 -2*alpha*beta*gamma)/(1 - gamma**2))
                if gdr3_check1 > 5 or gdr3_check2 > 5:
                    flags['gaia_astro_star'] = True
        alert_temp = locus.alerts[-1].properties
        if alert_temp['ant_survey'] == 1 and ((alert_temp['ztf_sgscore1'] > 0.5 and alert_temp['ztf_distpsnr1'] < 1) or (alert_temp['ztf_sgscore2'] > 0.5 and alert_temp['ztf_distpsnr2'] < 1) or (alert_temp['ztf_sgscore3'] > 0.5 and alert_temp['ztf_distpsnr3'] < 1)):
            flags['ps1_pointlike'] = True
    return flags


#query ANTARES alert stream
#goes by 10x10 deg^2 chunks, retries until it happens
def query_by_chunk(ra_min, ra_max, dec_min, dec_max, mjd_start=Time.now().mjd-1, mjd_end=Time.now().mjd, 
                   faint_lim=21, min_ndet=3, max_results=10000,
                   retries=3, wait=180):
    # set the query to fetch loci with alerts matching search criteria
    query = {
        "query": {
            "bool": {
                "filter": [
                    { "range": { "properties.newest_alert_observation_time": { "gte": mjd_start, "lte": mjd_end} } },
                    { "range": { "properties.brightest_alert_magnitude": {"lte": faint_lim} } },
                    { "range": { "properties.num_mag_values" : { "gte": 3} } },
                    { "range": { "ra": { "gte": ra_min, "lte": ra_max} } },
                    { "range": { "dec": { "gte": dec_min, "lte": dec_max} } }
                      ]
            }
        }
    }

    # sometimes the ANTARES server times out so try a few times
    # before giving up and returning an empty list (to avoid crashing the entire query)
    for attempt in range(retries):
        try:
            loci = list(search.search(query))
            if len(loci) >= max_results:
                warnings.warn(
                    f'Chunk MJD {mjd_start:.1f}–{mjd_end:.1f} hit the {max_results}-result max.'
                    f'— Reduce QUERY_DT or this window is incomplete.',
                    stacklevel=2
                )
            return loci
##            return get_ids(loci)
        except Exception as e:
            if attempt < retries - 1:
                print(f'  Chunk {mjd_start}–{mjd_end} failed (attempt {attempt+1}/{retries}): {e}.')
                time.sleep(15)
            else:
                time.sleep(wait)
                # pdb.set_trace()
                attempt = 0


#gets info for ToO target generation, one oid at a time
def get_target_info(locus, oid):
    target = pd.DataFrame(columns = ['too_id','fiber_type','catalogid','sdss_id','gaia_dr3_source_id',
                                     'twomass_pts_key','sky_brightness_mode','ra','dec','pmra','pmdec',
                                     'epoch','parallax','lambda_eff','u_mag','g_mag','r_mag','i_mag','z_mag',
                                     'optical_prov','gaia_bp_mag','gaia_rp_mag','gaia_g_mag','h_mag',
                                     'delta_ra','delta_dec','can_offset','inertial','n_exposures','priority',
                                     'active','expiration_date','observed', 'added_on'])
    for cat in locus.catalog_objects:
        if cat == 'sdss_gals':
            target.loc[0,'sdss_id'] = locus.catalog_objects[cat][0]['Objid']
        elif cat == 'sdss_stars':
            target.loc[0,'sdss_id'] = locus.catalog_objects[cat][0]['Objid']
        elif cat == 'gaia_dr3_gaia_source':
            target.loc[0,'gaia_dr3_source_id'] = locus.catalog_objects[cat][0]['source_id']
            target.loc[0, 'parallax'] = locus.catalog_objects[cat][0]['parallax']
            target.loc[0, 'pmdec'] = locus.catalog_objects[cat][0]['pmdec']
            target.loc[0, 'pmra'] = locus.catalog_objects[cat][0]['pmra']
            target.loc[0, 'gaia_bp_mag'] = locus.catalog_objects[cat][0]['phot_bp_mean_mag']
            target.loc[0, 'gaia_rp_mag'] = locus.catalog_objects[cat][0]['phot_rp_mean_mag']
            target.loc[0, 'gaia_g_mag'] = locus.catalog_objects[cat][0]['phot_g_mean_mag']
            target.loc[0, 'epoch'] = locus.catalog_objects[cat][0]['ref_epoch']
        elif cat == '2mass_xsc':
            target.loc[0,'twomass_pts_key']= locus.catalog_objects[cat][0]['id']

    band = str(oid['band'])
    if band[0] == 'l' and band[-1] != 'y':
        target.loc[0, band[-1]+'_mag'] = oid['m_sci']
    if band[0] == 'z':
        if band[-1] == 'R':
            target.loc[0,'r_mag'] = oid['m_sci']
        else:
            target.loc[0, band[-1]+'_mag'] = oid['m_sci']
            
    target.loc[0,'fiber_type'] = 'BOSS'
    target.loc[0,'sky_brightness_mode'] = 'dark'
    target.loc[0, 'ra'] = locus.ra
    target.loc[0,'dec'] = locus.dec
    target.loc[0, 'added_on'] = Time.now().to_value('iso', subfmt='date')
    target.loc[0, 'can_offset'] = True
    target.loc[0, 'n_exposures'] = 3
    target.loc[0, 'priority'] = 5
    target.loc[0, 'active'] = True
    target.loc[0, 'observed'] = False
    
    return target
