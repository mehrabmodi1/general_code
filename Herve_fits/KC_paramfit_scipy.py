#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Output fitting with the scipy curve_fit function """
import os
import pathlib
# import sys
import logging

import numpy as np
from scipy.optimize import curve_fit
import scipy.io as io
import matplotlib.pyplot as plt

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ["Helvetica Neue"]


# Global variables
FRAME_RATE = 1 / 0.091
FRAME_RATE = 1 / 0.0999


def gen_optifunc(stim_lat, stim_dur):
    """ Generate the function to be optimized """
    def f_kc(t, a0, a1, a2, a3, t0, t1, t2, t3):
        """outputs the defined function
        Parameters
        ----------
        params : TODO

        Returns
        -------
        TODO
        """
        domain0 = (t <= stim_lat)
        domain1 = (t > stim_lat) & (t <= stim_lat + stim_dur)
        domain2 = (t > stim_lat + stim_dur)
        
        out = domain0 * a0
        out += domain1 * (
            a1 * (1.0 - np.exp(-(t - stim_lat) / t0))
            + a2 * (1.0 - np.exp(-(t - stim_lat) / t1))
            + a0
        )
        end1 = (
            a1 * (1.0 - np.exp(-stim_dur / t0))
            + a2 * (1.0 - np.exp(-stim_dur / t1))
            + a0
        )
        out += domain2 * (
            a3 * np.exp(-(t - stim_lat - stim_dur) / t2)
            + (end1 - a3) * np.exp(-(t - stim_lat - stim_dur) / t3)
        )
        return out

    def frac_longt(a0, a1, a2, a3, t0, t1, t2, t3):
        """TODO: Docstring for frac_longt.

        Parameters
        ----------
        t : TODO
        a0 : TODO
        a1 : TODO
        a2 : TODO
        a3 : TODO
        t0 : TODO
        t1 : TODO
        t2 : TODO
        t3 : TODO

        Returns
        -------
        TODO

        """
        thr_rat = 1.5

        totint = np.abs(
            a1 * t0 * (1 - np.exp(-stim_dur / t0))
            + a2 * t1 * (1 - np.exp(-stim_dur / t1))
            + (a0 + a1 + a2) * stim_dur
        )

        if t0 / t1 < thr_rat and t1 / t0 < thr_rat:
            frac_2exp = np.abs(
                a1 * (t0 - (t0 + stim_dur) * np.exp(-stim_dur / t0))
                + a2 * (t1 - (t1 + stim_dur) * np.exp(-stim_dur / t1))
            )
            return frac_2exp / totint, np.max([t0, t1])

        if t0 > t1:
            ti = t0
            amp = a1
        else:
            ti = t1
            amp = a2
        # logging.info("ti and amp", ti, amp)
        # logging.info("stim_dur", stim_dur)
        # logging.info("exp", np.exp(-stim_dur / ti))
        frac_longexp = np.abs(
            amp * (ti - (ti + stim_dur) * np.exp(-stim_dur / ti))
        )
        return frac_longexp / totint, ti

    return f_kc, frac_longt


def train(activity, stim_lat, stim_dur, fold_fit):
    
    stim_lat = float(stim_lat);
    stim_dur = float(stim_dur);
    
    """ performing the fit
    Parameters
    ----------
    params :
    activity : the actual array
    stim_lat : stimulus latency
    stim_dur : stimulus duration
    fn : filename where to save the fitting parameters
    """

    ncell = activity.shape[1]
    logging.info("number of cells %d", ncell)

    # mushroom body network
    act_dur = activity.shape[2]
    times = np.arange(act_dur) / FRAME_RATE
    logging.info(times)
    domain1 = times < stim_lat
    domain2 = (stim_lat <= times) & (times < stim_lat + stim_dur)
    domain3 = times >= stim_lat + stim_dur

    logging.info(domain1)
    logging.info(domain2)
    logging.info(domain3)

    # fig, ax = plt.subplots()
    # ax.plot(times, domain1)
    # ax.plot(times, domain2)
    # ax.plot(times, domain3)

    f_kc, frac_longt = gen_optifunc(stim_lat, stim_dur)
    # test_out = f_kc(times, 1.0, 2.0, 5.0, 3.0, 0.5, 2.5, 1.0, 10.0)
    # ax.plot(times, test_out)
    # plt.show()
    # sys.exit()

    # print(activity[0, :-1])
    # print(activity[0, :].shape)
    # for j in np.arange(3):
    fit_params = np.zeros((activity.shape[0], activity.shape[1], 8))
    for j in np.arange(activity.shape[1]):
        #fig, ax = plt.subplots(1, activity.shape[0], figsize=(12, 4))
        for i in np.arange(activity.shape[0]):
            try:
                popt, pcov = curve_fit(
                    f_kc,
                    times[:-1],
                    activity[i, j, :-1],
                    bounds=([-1e3, -1e3, -1e3, -1e3, 0.3, 0.3, 0.3, 0.3],
                            [1e3, 1e3, 1e3, 1e3, 30.0, 30.0, 30.0, 30.0])
                )
            except RuntimeError:
                continue
            fit_params[i, j, :] = popt

            #ax[i].plot(times[:-1], activity[i, j, :-1])
            #ax[i].plot(times[:-1], f_kc(times[:-1], *popt))
            
        #fig.tight_layout()
        #fig.savefig(os.path.join(fold_fit, "cell{}.svg".format(j)))

    fn_fits = os.path.join(fold_fit, "fit_params.npy")
    np.save(fn_fits, fit_params)
    

    # coefs = tf.Variable(init_c, name="filt_coefs")
    # xran = tf.reshape(tf.range(sig_len, dtype=tf.float32), (1, sig_len))
    # filt_dec = tf.exp(-tf.multiply(times, xran))
    # coeftimes = coefs * times[..., 0]
    # print(filt_dec.shape)
    # filt = tf.einsum("abc,abcd->abd", coeftimes, filt_dec)
    # xmushfft = tf.spectral.rfft(xmush)
    # print(xmush.shape)
    # filtfft = tf.spectral.rfft(filt)
    # outputh1 = tf.nn.relu(tf.spectral.irfft(xmushfft * filtfft))
    # init_c2 = tf.truncated_normal((ncell, 2))
    # coefs2 = tf.Variable(init_c2, name="filt_coefs2")
    # output = tf.einsum("ab,abc->ac", coefs2, outputh1)

    logging.info("starting optimization")


def dataset_analysis(activity, odor, nonan, param, fold_fit):
    """TODO: Docstring for dataset_analysis.

    Parameters
    ----------
    activity : TODO
    params : TODO
    lobe : TODO
    fly : TODO

    Returns
    -------
    TODO

    """
    # importing the dataset
    stim_lat = param['stim_latency'][0, 0][0]
    stim_durs = param['odor_duration'][0, nonan]      #list of all odor durations
    odor_ns = param['odor_n'][0, nonan]               #list of all odor numbers delivered, not odor indices
    dur_list = np.unique(stim_durs);
    odor_list = np.unique(odor_ns);
    n_odors = odor_list.shape;                        #number of odors presented
    n_odors = n_odors[0];
    stim_long = np.max(dur_list);       #longest stimulus duration delivered
    i_stilong = (stim_durs == stim_long)
    i_stilong = np.squeeze(i_stilong);
    
    logging.info(
        "activity shape (stim index, KC index, trace): %s",
        activity.shape
    )
   
    activity = activity[i_stilong, ...]     #activity traces only for longest dur trials, throwing away all other trials - not needed for fitting.     
    odor_long = odor[i_stilong]             #list of odor indices (not odor numbers) delivered on each long duration trial.
    
    logging.info("odors for the long stimuli: %s", odor_long)
    acti_perodor = np.zeros((3, activity.shape[1], activity.shape[2]))
    for odor_ni in range(0, n_odors):
        acti_perodor[odor_ni, ...] = np.mean(activity[odor_long == odor_ni], axis=0)        #mean response traces for each odor, long dur
        
    logging.info("time latency: %s", stim_lat)
    
    # looking at the synchronization
    n_stim = activity.shape[0]
    act_dur = activity.shape[2]
    times = np.arange(act_dur) / FRAME_RATE
    lline = n_stim // 2 + 1
    """
    fig, ax = plt.subplots(2, lline, figsize=(16, 4))
    xlim = (times > 20.0) & (times < 30.0)
    for i in np.arange(n_stim):
        ax[i // lline, i % lline].plot(
            times[xlim], np.mean(activity[i, ...], axis=0)[xlim]
        )
    #fig.tight_layout()
    #fig.savefig(os.path.join(fold_fit, "start_stim.svg"))

    fig, ax = plt.subplots(2, lline, figsize=(16, 4))
    xlim = (times > 70.0) & (times < 100.0)
    for i in np.arange(n_stim):
        ax[i // lline, i % lline].plot(
            times[xlim], np.mean(activity[i, ...], axis=0)[xlim]
        )
    fig.tight_layout()
    #fig.savefig(os.path.join(fold_fit, "end_stim.svg"))
    plt.show()
    # sys.exit()
    """
    
    stim_lat = stim_lat + 0.5;
    # stim_long = 56.5
    train(acti_perodor, stim_lat, stim_long, fold_fit)

    # opening the dataset, the fits and plotting them
    """
    fn_fits = os.path.join(fold_fit, "fit_params.npy")
    fits_params = np.load(fn_fits)
    f_kc, frac_longt = gen_optifunc(stim_lat, stim_long)
    
    fig, ax = plt.subplots()
    taus = np.max(fits_params[..., (4, 5)], axis=2)
    # print(taus)
    ax.hist(taus.flatten(), bins=20)
    plt.show()

    fits_all = fits_params.reshape((-1, 8))
    outs = np.zeros((fits_all.shape[0], 2))
    valid_par = np.zeros(fits_all.shape[0], dtype=bool)
    for i in np.arange(fits_all.shape[0]):
        outs[i, 0], outs[i, 1] = frac_longt(*fits_all[i, :])
        if outs[i, 0] > 0.1 and outs[i, 0] < 1.0 and outs[i, 1] < 25:
            valid_par[i] = True

    fig, ax = plt.subplots(2, 1, figsize=(4, 8))
    ax[0].semilogx(outs[:, 0], outs[:, 1], '.')
    ax[1].loglog(outs[:, 0], fits_all[:, 4] / fits_all[:, 5], '.')

    fig, ax = plt.subplots(figsize=(4, 3))
    ax.hist(outs[valid_par, 1], bins=20)
    ax.set_xlabel("$\\tau_1^\\mathrm{off}$ (s)")
    ax.set_ylabel("# cells")
    fig.tight_layout()
    fig.savefig(os.path.join(fold_fit, "hist_tauoff.svg"))
    plt.show()
    """

def main():
    """ Main function (supervises the optimization)
    """
   
    # configures the logging (level of messages)
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
    from pathlib import Path, PureWindowsPath
    
    base_path = PureWindowsPath("C:\Data\Data\Analysed_data\data_sharing\KC_long_trace");
    base_path = Path(base_path);
    obj_list = os.listdir(base_path);
    n_objs = len(obj_list);
    n_dirs = 0;
    for obj_n in range(0, (n_objs - 1) ):
        curr_path = "{}/{}".format(base_path, obj_list[obj_n]);
        if os.path.isdir(curr_path) == 1:
            n_dirs = n_dirs + 1;
   
    n_flies = n_dirs;
    for fly_n in range(1, (n_flies)):
        act_fn = "{}/fly{}/dFF_data.mat".format(base_path, fly_n)
        para_fn = "{}/fly{}/stim_mat".format(base_path, fly_n)
        activity = io.loadmat(act_fn)['dff_data_mat_f'].T
        nonan = ~np.isnan(activity[:, :, 0:400]).any(axis=(1, 2))
        activity = activity[nonan, ...]
        params = io.loadmat(para_fn)['stim_mat']

        
        odor = params['odor_n'][0, nonan]
        odor[odor == 3] = 0
        odor[odor == 10] = 1
        odor[odor == 11] = 2

        folder_fit = "{}/fly{}".format(base_path, fly_n)
        
        pathlib.Path(folder_fit).mkdir(exist_ok=True)
        dataset_analysis(activity, odor, nonan, params, folder_fit)


if __name__ == "__main__":
    # execute only if run as a script
    main()
