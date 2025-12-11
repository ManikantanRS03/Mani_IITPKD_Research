# #function to load data
# import os
# def load_data(file_path, start_i, end_i):
#     PV_data = pd.read_excel(file_path)
#     P = np.array(PV_data.iloc[start_i:end_i, 3].to_list())*1e-2 # 1e-6/1e-2/1e-2 (C/m2 conversion)
#     V = np.array(PV_data.iloc[start_i:end_i, 2].to_list())
#     t = np.array(PV_data.iloc[start_i:end_i, 1].to_list())*1e-3
#     return P,V,t

# def save_sr_metrics(output_dir, f, A,
#                         sd_crt, bw_crt, Dext_delF, Dext_delF_crt,
#                         PSR, NF, SNR, COV,
#                         filename=None):
#         """
#         Save SR metrics arrays as JSON under output_dir/analysis.
#         Returns the saved file path.
#         """

#         analysis_dir = os.path.join(output_dir, 'analysis')
#         os.makedirs(analysis_dir, exist_ok=True)

#         outname = filename if filename is not None else f"sr_metrics_f{f}Hz_A{A}V.json"
#         outpath = os.path.join(analysis_dir, outname)

#         data = {
#             "sd_crt": np.asarray(sd_crt).tolist(),
#             "bw_crt": np.asarray(bw_crt).tolist(),
#             "Dext_delF": np.asarray(Dext_delF).tolist(),
#             "Dext_delF_crt": np.asarray(Dext_delF_crt).tolist(),
#             "PSR": np.asarray(PSR).tolist(),
#             "NF": np.asarray(NF).tolist(),
#             "SNR": np.asarray(SNR).tolist(),
#             "COV": np.asarray(COV).tolist()
#         }

#         with open(outpath, "w") as jf:
#             json.dump(data, jf, indent=2)

# def SR_matrics_exp(sd, V, P, V0, Pswitch, bw, f, alpha, beta, rho, tf, Af, delw=5, delw_noise = 70):
#     # calculating the matrics
#     PSD = []
#     PSR = []
#     NF = []
#     SNR = []
#     COV = []
#     sd_crt = []
#     bw_crt = []

#     for i in range(len(sd)):
#         Cov = np.cov(P[i],Pswitch)[1,0]    
#         freq, psd = sr.signal_psd(P[i], f, bw, delw)
#         freq, psd_noise = sr.signal_psd((V[i] - V0), f, bw, delw)
#         Pwr, Nf, Snr = sr.power_snr_noisefloor_iisc(freq, psd , f,  bw, delw, delw_noise)
    
#         PSD.append(psd)
#         PSR.append(Pwr)
#         NF.append(Nf)
#         SNR.append(Snr)
#         COV.append(Cov)
#         sd_crt.append(sr.sd_correction(V[i],V0))
#         bw_crt.append(sr.bandwith_correction(freq, psd_noise)[1])

#     sd_crt = np.array(sd_crt)
#     bw_crt = np.array(bw_crt)
#     Dext_delF = np.array(sr.Dext_delF(sd, bw, alpha, beta, rho, tf, Af))
#     Dext_delF_crt = np.array(sr.Dext_delF(sd_crt, bw, alpha, beta, rho, tf, Af))
#     return sd_crt, bw_crt, Dext_delF, Dext_delF_crt, PSR, NF, SNR, COV
# 

# IISc_data_path = "C:\\Drive\\FE_SR\\IITPKD_data_analysis\\sim_for_final_figs\\IISc_data\\"

# #load all the data's 75 Hz, 100Hz, 150 Hz
# sd = np.linspace(0.0,5.1,35)
# P_data = []
# V_data = []
# t_data = []

# start_i, end_i = [45, 45, 45], [20067, 17592, 15042] #index of the data
# f_data = [75, 100, 150]

# #index after the poling pulse
# fro = int(bw + 48)

# for j in range(0,3):
#     P_f = []
#     V_f = []
#     t_f = []
#     for i in range(0,len(sd)):
#         file_path = IISc_data_path + f'new_{f_data[j]}Hz\\{sd[i]:.2f}.xls'    
#         Ptemp, Vtemp, ttemp = load_data(file_path, start_i[j], end_i[j])    
#         P_f.append(np.array(Ptemp[fro:]))
#         V_f.append(np.array(Vtemp[fro:]))
#         t_f.append(np.array(ttemp[fro:]))
#     P_data.append(P_f)
#     V_data.append(V_f)
#     t_data.append(t_f)

# # reading Hyst data (for the super threshold signal)
# PV_data = pd.read_excel(IISc_data_path + "hyst_1000hz.xlsx")
# Pswitch = np.array(PV_data.iloc[46:43+2002, 3].to_list())*1e-2 # 1e-6/1e-2/1e-2 (C/m2 conversion)
# tswitch = np.array(PV_data.iloc[46:43+2002, 1].to_list())*1e-3

# for i in range(len(f_data)):
#     V = V_data[i]
#     V0 = V_data[i][0]
#     P = P_data[i]
    
#     # code to make the Pswitch signal
#     P_switch = list(Pswitch*2/3)*75
#     time_switch = np.linspace(t_data[i][0][0], t_data[i][0][-1], len(P_switch))
#     time_switch_gen = t_data[i][0]
#     interp_func = interp1d(time_switch, P_switch, kind='linear', fill_value='extrapolate')
#     P_switch_gen = interp_func(time_switch_gen)

#     sd_crt_exp, bw_crt_exp, Dext_delF_exp, Dext_delF_crt_exp, PSR_exp, NF_exp, SNR_exp, COV_exp = SR_matrics_exp(sd, V, P, V0, P_switch_gen, bw, f_data[i], alpha, beta, rho, tf, Af, delw=4, delw_noise = 70)
    
#     #sorting
#     idx_sort = np.argsort(Dext_delF_crt_exp)
#     Dext_delF_crt_exp = Dext_delF_crt_exp[idx_sort]
#     PSR_exp = np.array(PSR_exp)[idx_sort]
#     NF_exp = np.array(NF_exp)[idx_sort]
#     SNR_exp = np.array(SNR_exp)[idx_sort]
#     COV_exp = np.array(COV_exp)[idx_sort]

#     save_sr_metrics(IISc_data_path, f_data[i], 1.575 , sd_crt_exp, bw_crt_exp, Dext_delF_exp, Dext_delF_crt_exp, PSR_exp, NF_exp, SNR_exp, COV_exp, filename=None)