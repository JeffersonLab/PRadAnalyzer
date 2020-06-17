import os
import errno
import pandas as pd
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
import argparse


parser = argparse.ArgumentParser('resolution estimate')
parser.add_argument('sim_res', help='simulation output file')
parser.add_argument('sim_gun', help='simulation output gun file')
parser.add_argument('-o', dest='outdir', help='output directory', default='.')
parser.add_argument('-e', '--energy', dest='energy', type=float, help='simulation beam energy', default=700)
parser.add_argument('-z', '--z-hycal', dest='hycal_z', type=float, help='HyCal z position (rel. to target)', default=2950 + 3000 - 88.9)
parser.add_argument('--prad-results', dest='prad_xs', help='path to PRad xs results', default='PRad_xs.txt')
parser.add_argument('--bootstrap-portion', dest='bport', help='portion of events for each resampling', type=int, default=5)
parser.add_argument('--bootstrap-times', dest='btimes', help='number of resampling times for bootstrapping', type=int, default=1000)

args = parser.parse_args()

energy = args.energy
hycal_z = args.hycal_z
base_res = 0.026


def calc_theta(x, y, z=hycal_z):
    r = np.sqrt(np.square(x) + np.square(y))
    return np.arctan(r/z)


def data_process(path, min_energy=20.):
    df = pd.read_csv(path)
    df = df[df['E'] > min_energy]
    df.loc[:, 'nhits'] = df.groupby('Event Number')['E'].transform('count')
    df.loc[:, 'matched'] = df['match'] > 0
    df.loc[:, 'nmatches'] = df.groupby('Event Number')['matched'].transform('sum')
    df.loc[:, 'totalE'] = df.groupby('Event Number')['E'].transform('sum')
    df.loc[:, 'scat_angle'] = calc_theta(df['x'], df['y'])/np.pi*180.
    df.loc[:, 'phi_angle'] = np.arctan2(df['x'], df['y'])/np.pi*180.
    return df


def elastic_cut(df, e0, base_res, sigma=3.0):
    energies = df['E'].values
    thetas = df['scat_angle'].values / 180 * np.pi
    expect = e0 / (1 + 2 * e0 / 938.272 * np.square(np.sin(thetas / 2.)))
    res = base_res/np.sqrt(energies/1000.)
    return (df['nmatches'] == 1) & (df['match'] > 0) &\
           (energies < expect*(1. + sigma*res)) & (energies > expect*(1. - sigma*res))


def range_cut(values, arange):
    return (values <= arange[1]) & (values >= arange[0])


def azimuthal_cut(angles, cuts):
    mask = np.array([False]*len(angles))
    for cut in cuts:
        mask |= (angles <= cut[1]) & (angles >= cut[0])
    return mask


def gaus(x, p0, p1, p2):
    return p0*np.exp(-(x-p1)**2/(2*p2**2))


# from PRad results
prad_res = pd.read_csv(args.prad_xs, header=None, sep=r'\s+', comment='#')
prad_res.columns = ['E', 'theta', 'theta_bw', 'Q2', 'xs', 'xs_stat', 'xs_syst', 'f', 'f_stat', 'f_syst']

# data processing
dfep = data_process(args.sim_res)
dfg = pd.read_csv(args.sim_gun)
dfg['scat_angle'] = dfg['theta']/np.pi*180.

# for plots
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
bprops = dict(boxstyle='round', facecolor='wheat', alpha=1.0)

# select elastic ep events, resolution depends on region (lead-glass or PbWO4)
res = np.where(dfep['cid'] < 1000, 1.0, 2.5)*base_res
elas = dfep[elastic_cut(dfep, energy, res, 3.5)]

# projections
mask1 = prad_res['E'] == 1101
mask2 = prad_res['E'] == 2143

beam_table = [
    (700, prad_res[mask1]),
    (1400, prad_res[mask1]),
    (2100, prad_res[mask2]),
]

dbins = np.arange(-0.2, 0.2, step=0.005)
didx = (dbins[1:] + dbins[:-1])/2.
bwidth = np.min(np.diff(didx))

angular_res = pd.DataFrame(columns=['energy', 'theta', 'theta_bw', 'angular_res', 'angular_bc_res'])
for eproj, pres in beam_table:
    extra_bins = np.array([(a - 0.025, a + 0.025) for a in [0.525, 0.575, 0.625, 0.675]])
    abins = np.vstack([pres['theta'] - pres['theta_bw'], pres['theta'] + pres['theta_bw']]).T
    abins = np.concatenate((extra_bins, abins))
    escale = energy/eproj
    for abin in abins:
        elbin = elas[range_cut(elas['scat_angle'], abin)]
        expval = elbin['scat_angle'].values
        # find the generated events
        evgen = dfg[dfg['Event Number'].isin(elbin['Event Number'].unique()) & (dfg['pid'] == 11)].copy()
        # only keep the electron that has the larger energy (there should be no duplicates for ep)
        realval = evgen.sort_values('E', ascending=False).drop_duplicates('Event Number')['scat_angle'].values

        # # angular resolution
        vals, _ = np.histogram(expval - realval, bins=dbins)
        res, _ = curve_fit(gaus, didx, vals, p0=[1000, 0., 0.05])
        # print('{:.4f}, {:.4f}, {:.5f}'.format(*abin, np.abs(res[2])))

        # bin center resolution, check by bootstrapping
        # randome pick 1/5 samples for 1000 times and check the statistical bin center difference
        diff = []
        for i in np.arange(args.btimes):
            # experimental bin center
            nsample = int(len(expval)/args.bport)
            vals, abins = np.histogram(np.random.choice(expval, nsample), np.linspace(*abin, 20))
            exp_bc = np.average((abins[-1:] + abins[1:])/2., weights=vals)
            # real bin center
            # vals, abins = np.histogram(realval, np.linspace(*abin, 20))
            # real_bc = np.average((abins[-1:] + abins[1:])/2., weights=vals)
            # geomterical center
            real_bc = np.average(abin)
            diff.append(exp_bc - real_bc)
        print('E = {:.1f} MeV, angle range = {:.3f} - {:.3f}, angular res. = {:.4f}, bin center res. = {:.7f}'.format(
            eproj, abin[0], abin[1], res[2]*escale, np.std(diff)*escale)
        )
        angular_res.loc[len(angular_res)] = (eproj, np.average(abin), (abin[1] - abin[0])/2., res[2]*escale, np.std(diff)*escale)


angular_res.to_csv(os.path.join(args.outdir, 'angular_res_estimate.csv'), index=False)
