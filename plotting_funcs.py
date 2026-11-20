from filter_funcs import *
from antares_client import search
import matplotlib.pyplot as plt
import math


#make lightcurve from all alerts in a locus
def lightcurve(alerts):
    lc = 0
    for alert in alerts:
        if alert.properties['ant_survey'] == 4:
            data = lsst_alert(alert)
        if alert.properties['ant_survey'] == 1:
            data = ztf_alert(alert)
        else:
            continue
        if isinstance(data,dict) == 1:
            if isinstance(lc, int) == 1:
                lc = pd.DataFrame(data, index=[0])
            else:
                h = pd.DataFrame(data, index=[-1])
                lc = pd.concat([lc, h]).reset_index(drop=True)
    return lc


#make lightcuve from only quality alerts in a locus
def lightcurve_quality(alerts):
    lc = 0
    for alert in alerts:
        if alert.properties['ant_survey'] == 4 and lsst_quality(alert) == True:
            data = lsst_alert(alert)
        elif alert.properties['ant_survey'] == 1 and ztf_quality(alert) == True:
            data = ztf_alert(alert)
        else:
            continue
        if isinstance(data,dict) == 1:
            if isinstance(lc, int) == 1:
                lc = pd.DataFrame(data, index=[0])
            else:
                h = pd.DataFrame(data, index=[-1])
                lc = pd.concat([lc, h]).reset_index(drop=True)
    return lc


#plot lightcurve using only quality alerts
def plot_lc(oid_id):
    colors = {'ztf_g': 'green', 'ztf_R': 'orangered', 'ztf_i': 'black'}
    fig, ax = plt.subplots(1,1, figsize=(9,6))
    
    locus = search.get_by_id(oid_id)
    lc = lightcurve_quality(locus.alerts)
    title = oid_id
    flags_str = 'Flags: '
    
    for band in colors:
        band_data = lc[(lc['band'] == band)]
        if len(band_data) > 2:
            ax.errorbar(band_data['mjd'], band_data['m_sci'], yerr=band_data['m_sci_err'], fmt='.', ms=2, lw=1, 
                        c = colors[band], ecolor = colors[band], label=band, alpha = 0.5)
            ax.hlines(y=np.average(band_data['m_template']), xmin = locus.properties['oldest_alert_observation_time'], 
                      xmax= locus.properties['newest_alert_observation_time'], colors = colors[band], alpha = 0.3, ls = '--')
            
    flags = catalog_check(locus)
    for flag in flags:
        if flags[flag] == True and len(flags_str) == 7:
            flags_str += flag
        elif flags[flag] == True:
            flags_str += ', ' + flag
    ax.set(ylim=(min(lc['m_sci'])-0.5, max(lc['m_sci'])+0.5))
    ax.yaxis.set_inverted(True)
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)
        
    if len(flags_str) > 7:
        fig.suptitle(title + '\n' + flags_str, fontname = 'serif')
    else:
        fig.suptitle(title + '\n' + flags_str + 'None', fontname = 'serif')
    plt.xlabel('MJD', fontname = 'serif')
    ax.set_ylabel('$mag_{AB}$', fontname = 'serif')
    ax.grid(alpha=0.25)
    plt.tight_layout()
    plt.show()

#same as plot_lc but also plots using all data
def plot_lc2(oid_id):
    colors = {'ztf_g': 'green', 'ztf_R': 'orangered', 'ztf_i': 'black'}
    fig, ax = plt.subplots(1,2, figsize=(10,4), sharey=False)
    
    locus = search.get_by_id(oid_id)
    lc = lightcurve(locus.alerts)
    lc_q = lightcurve_quality(locus.alerts)
    title = oid_id
    flags_str = 'Flags: '

    empty = 0
    empty_q = 0
    
    for band in colors:
        band_data = lc[(lc['band'] == band)]
        band_data_q = lc_q[(lc_q['band'] == band)]
        if len(band_data) > 2:
            ax[0].errorbar(band_data['mjd'], band_data['m_sci'], yerr=band_data['m_sci_err'], fmt='.', ms=2, lw=1, c = colors[band], ecolor = colors[band], label=band, alpha = 0.5)
            ax[0].hlines(y=np.average(band_data['m_template']), xmin = locus.properties['oldest_alert_observation_time'], xmax= locus.properties['newest_alert_observation_time'], 
                          colors = colors[band], alpha = 0.3, ls = '--')
        else:
            empty += 1
        if len(band_data_q) > 2:
            ax[1].errorbar(band_data_q['mjd'], band_data_q['m_sci'], yerr=band_data_q['m_sci_err'], fmt='.', ms=2, lw=1, c = colors[band], ecolor = colors[band], label=band, alpha = 0.5)
            ax[1].hlines(y=np.average(band_data_q['m_template']), xmin = locus.properties['oldest_alert_observation_time'], xmax= locus.properties['newest_alert_observation_time'], 
                          colors = colors[band], alpha = 0.3, ls = '--')
        else:
            empty_q += 1
    
    if empty == 3:
        print(f'{oid_id}: not enough data')
        plt.close()
    else:
        flags = catalog_check(locus)
        for flag in flags:
            if flags[flag] == True and len(flags_str) == 7:
                flags_str += flag
            elif flags[flag] == True:
                flags_str += ', ' + flag
        ax[0].set(ylim=(min(lc['m_sci'])-0.5, max(lc['m_sci'])+0.5))
        ax[1].set(ylim=(min(lc_q['m_sci'])-0.5, max(lc_q['m_sci'])+0.5))
        ax[0].yaxis.set_inverted(True)
        ax[1].yaxis.set_inverted(True)
        ax[1].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)
    
        if len(flags_str) > 7:
            fig.suptitle(title + '\n' + flags_str, fontname = 'serif')
        else:
            fig.suptitle(title + '\n' + flags_str + 'None', fontname = 'serif')
        plt.xlabel('MJD', fontname = 'serif')
        ax[0].set_ylabel('$mag_{AB}$', fontname = 'serif')
        ax[1].set_ylabel('$mag_{AB}$', fontname = 'serif')
        ax[0].grid(alpha=0.25)
        ax[1].grid(alpha=0.25)
        plt.tight_layout()
        plt.show()
