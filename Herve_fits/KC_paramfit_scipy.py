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
import matplotlib as plt

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

    
    f_kc, frac_longt = gen_optifunc(stim_lat, stim_dur)
   
    fit_params = np.zeros((activity.shape[0], activity.shape[1], 8))
    for j in np.arange(activity.shape[1]):
        
        for i in np.arange(activity.shape[0]):
            try:
                popt, pcov = curve_fit(
                    f_kc,
                    times[:-1],
                    activity[i, j, :-1],
                    bounds=([-1e3, -1e3, -1e3, -1e3, 0.3, 0.3, 0.3, 0.3],
                            [1e3, 1e3, 1e3, 1e3, 70.0, 70.0, 70.0, 70.0])
                )
            except RuntimeError:
                continue
            fit_params[i, j, :] = popt

           

    fn_fits = os.path.join(fold_fit, "fit_params.npy")
    fn_fits_a = os.path.join(fold_fit, "fit_data.npy")
    np.save(fn_fits, fit_params)
    np.save(fn_fits_a, activity)
    logging.info("starting optimization")


def dataset_analysis(activity, odor, nonan, param, fold_fit):
    
    # importing the dataset
    stim_lat = param['stimLatency'][0, 0][0]
    stim_durs = param['duration'][0, nonan]      #list of all odor durations
    odor_ns = param['odours'][0, nonan]               #list of all odor numbers delivered, not odor indices
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
    acti_perodor = np.zeros((n_odors, activity.shape[1], activity.shape[2]))
    for odor_n in np.arange(n_odors):
       
        odor_ni = odor_list[odor_n];
        curr_trs = [odor_long == odor_ni];
        curr_trs = curr_trs[0];
        curr_trs = np.squeeze(np.transpose(curr_trs));
        acti_perodor[odor_n, ...] = np.mean(activity[curr_trs, :, :], axis=0)        #mean response traces for each odor, long dur
        
        logging.info("time latency: %s", stim_lat)
        
        # looking at the synchronization
        n_stim = activity.shape[0]
        act_dur = activity.shape[2]
        times = np.arange(act_dur) / FRAME_RATE
        lline = n_stim // 2 + 1
       
        
        stim_lat = stim_lat + 0.5;
        # stim_long = 56.5
        
        train(acti_perodor, stim_lat, stim_long, fold_fit)
       
   

def main():
    """ Main function (supervises the optimization)
    """
   
    # configures the logging (level of messages)
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
    from pathlib import Path, PureWindowsPath
    
    base_path = PureWindowsPath("C:\Data\Data\Analysed_data\data_sharing\param_fitting\ApBpKC_hover_rtrain_prehabituated");       #this is the base path specified by the user 
    base_path = Path(base_path);
    obj_list = os.listdir(base_path);
    n_objs = len(obj_list);
    n_dirs = 0;
    for obj_n in range(0, (n_objs) ):
        curr_path = "{}/{}".format(base_path, obj_list[obj_n]);
        if os.path.isdir(curr_path) == 1:
            n_dirs = n_dirs + 1;
   
   
    n_flies = n_dirs;
    for fly_n in range(1, (n_flies + 1)):
        mes = "extracting params for fly {} of {} flies".format(fly_n,  n_flies);
        print(mes)
        act_fn = "{}/fly{}/dFF_data.mat".format(base_path, fly_n);
        para_fn = "{}/fly{}/stim_mat".format(base_path, fly_n);
        activity = io.loadmat(act_fn)['dff_data_mat_f'].T;
        
        #downsizing activity in frame_n dim to get rid of nans at the end. Matching all trials to shortest real data trial.
        n_nans_max = 0;
        tr_list = np.arange(activity.shape[0]);
        bad_trs = [];
        for tr_n in np.arange(activity.shape[0]):
            curr_vec = activity[tr_n, 0, :];
            n_nans = np.sum(np.isnan(curr_vec));
            if n_nans > (activity.shape[2]*0.7):
                bad_trs = bad_trs + [tr_n];           #building list of bad trs, with more than 70% of the frames as nans
            elif n_nans <= (activity.shape[2]*0.7):
                n_nans_max = np.max([n_nans, n_nans_max]);  #max number of nans in any trial. Cropping all trials at the end by this frame count
           
        nonan  = np.delete(tr_list, bad_trs, 0);      #getting rid of bd trs from list of good trs
        activity_trimmed = np.delete(activity, np.arange((activity.shape[2] - n_nans_max), activity.shape[2]), 2);
                
        activity = activity_trimmed[nonan, ...];
        params = io.loadmat(para_fn)['stim_mat'];
        
                
        #checking if current fly has already been analysed and recovering if partially analysed
        an_file_path = "{}/fly{}/fit_params.npy".format(base_path, fly_n);
       
        if os.path.isfile(an_file_path) == 1: 
            n_cells = activity.shape;
            n_cells = n_cells[1];
            prev_fits = np.load(an_file_path);
            #identifying last analysed cell in case of an interruption
            prev_fits = np.sum(prev_fits, axis = 2);
            #continue
        
        
        odor = params['odours'][0, nonan]
        
        
        folder_fit = "{}/fly{}".format(base_path, fly_n)
        
        pathlib.Path(folder_fit).mkdir(exist_ok=True)
        dataset_analysis(activity, odor, nonan, params, folder_fit)
        mes = "done extracting params for fly {} of {} flies".format(fly_n,  n_flies);
        print(mes)
    
    breakpoint()
if __name__ == "__main__":
    # execute only if run as a script
    main()
