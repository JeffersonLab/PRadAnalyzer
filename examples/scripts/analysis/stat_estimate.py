import os
import argparse
import pandas as pd
import numpy as np
from scipy import interpolate


parser = argparse.ArgumentParser('statistical uncertainty estimate')
parser.add_argument('--prad-results', dest='prad_xs', help='path to PRad xs results', default='PRad_xs.txt')
parser.add_argument('--ep-xs', dest='ep_xs', help='path to ep xs estimate', default='ep_xs.dat')
parser.add_argument('-o', dest='outdir', help='output directory', default='.')

args = parser.parse_args()

target_thickness = 1.8e18  # atoms/cm^2, estimate from PRad
det_eff = 0.8   # detection efficiency


# from PRad results
prad_res = pd.read_csv(args.prad_xs, header=None, sep=r'\s+', comment='#')
prad_res.columns = ['E', 'theta', 'theta_bw', 'Q2', 'xs', 'xs_stat', 'xs_syst', 'f', 'f_stat', 'f_syst']

# projection on 0.7/1.4
mask1 = prad_res['E'] == 1101
mask2 = prad_res['E'] == 2143

xs_table = pd.read_csv(args.ep_xs, comment='#', header=None)
xs_table.columns = ['energy', 'angle', 'non_rad', 'rad', 'born']

# energy (MeV), current (nA), beam time (days), prad bins
beam_time_table = [
    (700, 20, 4, prad_res[mask1]),
    (1400, 70, 5, prad_res[mask1]),
    (2100, 70, 15, prad_res[mask2]),
]


stat_table = pd.DataFrame()
for e0, curr, days, ref in beam_time_table:
    mask = np.isclose(xs_table['energy'], e0)
    xsf = interpolate.interp1d(xs_table[mask]['angle'].values, xs_table[mask]['non_rad'].values/1e6)
    res = ref[['theta', 'theta_bw']].copy()
    # extra small angle bin
    extra = pd.DataFrame(columns=['theta', 'theta_bw'], data=[(ang, 0.025) for ang in [0.525, 0.575, 0.625, 0.675]])
    res = pd.concat([extra, res], axis=0)
    res.loc[:, 'energy (MeV)'] = e0
    res.loc[:, 'beam_current (nA)'] = curr
    res.loc[:, 'beam_time (days)'] = days
    res.loc[:, 'efficiency (%)'] = det_eff
    res.loc[:, 'target_thickness (1e18/cm^2)'] = target_thickness/1e18
    solid = (np.cos((res['theta'] - res['theta_bw'])/180.*np.pi) - np.cos((res['theta'] + res['theta_bw'])/180.*np.pi))*2.*np.pi
    # xs = np.reciprocal(np.square(e0 / ref['E'])) * ref['xs']
    xs = xsf(res['theta'].values)
    rates = xs*solid*(curr*1e-9/1.602e-19*target_thickness/1e27)
    yields = rates*days*24*3600*det_eff
    res.loc[:, 'xs (mb/sr)'] = xs
    res.loc[:, 'rates (Hz)'] = rates
    res.loc[:, 'yields'] = yields
    res.loc[:, 'acceptance'] = 1.0
    # extra bins only have 160/360 acceptance on azimuthal angle
    res.loc[res['theta'].isin(extra['theta']), 'acceptance'] = 160/360
    res.loc[:, 'stat (relative)'] = np.reciprocal(np.sqrt(res['yields']*res['acceptance']))

    print(res)
    stat_table = pd.concat([stat_table, res], axis=0)

stat_table.rename(columns={'theta': 'angular_bin_center', 'theta_bw': 'angular_bin_width'}).\
    to_csv(os.path.join(args.outdir, 'stat_estimate.csv'), index=False)
