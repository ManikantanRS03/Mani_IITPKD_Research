import SR_lib as sr
import numpy as np
import os
import pandas as pd
from scipy.interpolate import interp1d
import re
import json

def run(params_path):
    ## Read the Parse file
    params = sr.read_simulation_parameters(params_path)
    Vc = params.get('Vc')
    Pr = params.get('Pr')
    alpha = params.get('alpha')
    beta = params.get('beta')
    rho = params.get('rho')
    tf = params.get('tf')
    Af = params.get('Af')
    Asig = params.get('A')
    fsig = params.get('f')
    sd = params.get('sd')
    ts = params.get('ts')
    n = params.get('n')
    ens = params.get('ens')
    dsf = params.get('dsf')  # Downsampling factor
    bw = params.get('bw') 
    Ec = Vc/tf
    delF = sr.del_U(0, alpha, beta, tf)[0]

    ## Create the simulation directory
    id = params.get('id')
    output_dir = f'sim_{id}'
    os.makedirs(output_dir, exist_ok=True)
    raw_data_file = os.path.join(output_dir, "raw_data")

    def SR_matrics(params, f, A, raw_data_file, delw=5, delw_noise = 70):
        Vc = params.get('Vc')
        Pr = params.get('Pr')
        alpha = params.get('alpha')
        beta = params.get('beta')
        rho = params.get('rho')
        tf = params.get('tf')
        Af = params.get('Af')
        T = params.get('T')
        Asig = params.get('A')
        fsig = params.get('f')
        sd = params.get('sd')
        ts = params.get('ts')
        n = params.get('n')
        ens = params.get('ens')
        dsf = params.get('dsf')  # Downsampling factor
        bw = params.get('bw') 

        PSD = []
        PSR = []
        NF = []
        SNR = []
        COV = []
        sd_crt = []
        bw_crt = []

        t, Edrive, P = sr.read_tEP_data(raw_data_file + f'/PE_{f}hz_{A}V_sd{sd[0]:.2f}_ens{0}.txt')
        E0 = Edrive
        t, Edrive, P = sr.read_tEP_data(raw_data_file + f'/PE_{f}hz_{Vc*1.5:.2f}V_sd{sd[0]:.2f}_ens{0}.txt')
        Pswitch = P
        for i in range(len(sd)):
            Cov_ens = 0
            psd_ens = 0
            psd_noise_ens = 0
            for e in range(ens):
                t, Edrive, P = sr.read_tEP_data(raw_data_file + f'/PE_{f}hz_{A}V_sd{sd[i]:.2f}_ens{e}.txt')
                Cov_ens += np.cov(P,Pswitch)[1,0]    
                freq, psd = sr.signal_psd(P, f, 1/(dsf*ts), delw)
                freq, psd_noise = sr.signal_psd((Edrive - E0)*tf, f, 1/(dsf*ts), delw)
                psd_ens += psd
                psd_noise_ens += psd_noise

            P_ens, NF_ens, SNR_ens = sr.power_snr_noisefloor_iisc(freq, psd_ens/ens , f,  1/(dsf*ts), delw, delw_noise)
            # P_ens, NF_ens, SNR_ens, _, _ = sr.power_snr_noisefloor(freq, psd_ens/ens , f,  1/(dsf*ts), delw, delw_noise)

            ffro = int(2*len(freq)/ (1/(dsf*ts))*(f - delw_noise)) + 1
            fto = int(2*len(freq)/ (1/(dsf*ts))*(f + delw_noise)) + 1
            
            PSD.append(psd_ens/ens)
            PSR.append(P_ens)
            NF.append(NF_ens)
            SNR.append(SNR_ens)
            COV.append(Cov_ens/ens)
            sd_crt.append(sr.sd_correction(Edrive*tf,E0*tf))
            bw_crt.append(sr.bandwith_correction(freq, psd_noise_ens/ens)[1])

        sd_crt = np.array(sd_crt)
        bw_crt = np.array(bw_crt)
        Dext_delF = np.array(sr.Dext_delF(sd, bw, alpha, beta, rho, tf, Af))
        Dext_delF_crt = np.array(sr.Dext_delF(sd_crt, bw_crt, alpha, beta, rho, tf, Af))
        return sd_crt, bw_crt, Dext_delF, Dext_delF_crt, PSR, NF, SNR, COV

    def save_sr_metrics(output_dir, f, A,
                        sd_crt, bw_crt, Dext_delF, Dext_delF_crt,
                        PSR, NF, SNR, COV,
                        filename=None):
        """
        Save SR metrics arrays as JSON under output_dir/analysis.
        Returns the saved file path.
        """

        analysis_dir = os.path.join(output_dir, 'analysis')
        os.makedirs(analysis_dir, exist_ok=True)

        outname = filename if filename is not None else f"sr_metrics_f{f}Hz_A{A}V.json"
        outpath = os.path.join(analysis_dir, outname)

        data = {
            "sd_crt": np.asarray(sd_crt).tolist(),
            "bw_crt": np.asarray(bw_crt).tolist(),
            "Dext_delF": np.asarray(Dext_delF).tolist(),
            "Dext_delF_crt": np.asarray(Dext_delF_crt).tolist(),
            "PSR": np.asarray(PSR).tolist(),
            "NF": np.asarray(NF).tolist(),
            "SNR": np.asarray(SNR).tolist(),
            "COV": np.asarray(COV).tolist()
        }

        with open(outpath, "w") as jf:
            json.dump(data, jf, indent=2)

    for j in range(len(Asig)):
        for i in range(len(fsig)):
            f = fsig[i]
            A = Asig[j]
            sd_crt, bw_crt, Dext_delF, Dext_delF_crt, PSR, NF, SNR, COV = SR_matrics(params, f, A, raw_data_file, delw=5, delw_noise = 70)
            save_sr_metrics(output_dir, f, A, sd_crt, bw_crt, Dext_delF, Dext_delF_crt, PSR, NF, SNR, COV, filename=None)