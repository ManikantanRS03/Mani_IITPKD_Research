import os
import numpy as np
import matplotlib.pyplot as plt
import SR_lib as sr
import scipy as sc
from matplotlib.backends.backend_pdf import PdfPages

class Analysis:
    def __init__(self):
        pass
    def t_V_P_plots(self, params, sd, f, A, ens, raw_data_file, pdf):
        #Making time series plot
        Vc = params.get('Vc')
        Pr = params.get('Pr')
        tf = params.get('tf')
        Ec = Vc/tf
        t, Edrive, P = sr.read_tEP_data(raw_data_file + f'\\PE_{f}hz_{A}V_sd{sd:.2f}_ens{ens}.txt')
        # First 100 cycles
        fro = 0
        to = int(len(t)/5)
        # Plotting
        plt.figure(figsize=(10, 8))
        # Subplot 1: Edrive
        plt.subplot(2, 1, 1)
        plt.plot(t[fro:to], Edrive[fro:to], label='Edrive')
        plt.axhline(y=-Ec, linestyle='--', color='b', label='-Ec')
        plt.axhline(y=Ec, linestyle='--', color='b', label='Ec')
        plt.xlabel("Time (s)")
        plt.ylabel("E (V/m)")
        plt.title(f'Edrive Signal: {f} Hz, {A} V, {sd:.2f}')
        # plt.ylim([-2 * Ec, 2 * Ec])
        plt.legend()

        # Subplot 2: Polarization
        plt.subplot(2, 1, 2)
        plt.plot(t[fro:to], P[fro:to], label='P')
        plt.axhline(y=Pr, linestyle='--', color='r', label='+Pr')
        plt.axhline(y=-Pr, linestyle='--', color='r', label='-Pr')
        plt.xlabel("Time (s)")
        plt.ylabel("P (C/m²)")
        plt.title(f'Polarization Response')
        plt.ylim([-1.5 * Pr, 1.5 * Pr])
        plt.legend()

        plt.tight_layout()
        pdf.savefig()  # saves the current figure into the PDF
        plt.close()
    
    def SR_matrics(self, params, raw_data_file, pdf, delw =5 , delw_noise = 70):
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
        
        for f in fsig:
            for A in Asig:
                t, Edrive, P = sr.read_tEP_data(raw_data_file + f'\\PE_{f}hz_{A}V_sd{sd[0]:.2f}_ens{0}.txt')
                E0 = Edrive
                for i in range(len(sd)):
                    Cov_ens = 0
                    psd_ens = 0
                    psd_noise_ens = 0
                    for e in range(ens):
                        t, Edrive, P = sr.read_tEP_data(raw_data_file + f'\\PE_{f}hz_{A}V_sd{sd[i]:.2f}_ens{e}.txt')
                        Cov_ens += np.cov(P,E0*tf)[1,0]    
                        freq, psd = sr.signal_psd(P, f, 1/(dsf*ts), delw)
                        freq, psd_noise = sr.signal_psd((Edrive - E0)*tf, f, 1/(dsf*ts), delw)
                        psd_ens += psd
                        psd_noise_ens += psd_noise

                    P_ens, NF_ens, SNR_ens, m, c = sr.power_snr_noisefloor(freq, psd_ens/ens , f,  1/(dsf*ts), delw, delw_noise)
                    ffro = int(2*len(freq)/ (1/(dsf*ts))*(f - delw_noise)) + 1
                    fto = int(2*len(freq)/ (1/(dsf*ts))*(f + delw_noise)) + 1
                    
                    plt.figure(figsize=(8, 5))
                    plt.semilogy(freq[ffro:fto], psd_ens[ffro:fto]/ens)
                    plt.semilogy(freq[ffro:fto], m*freq[ffro:fto]+c)
                    plt.xlabel("Frequency (Hz)")
                    plt.ylabel("PSD")
                    plt.title(f"PSD_{f}hz_{A}V_{sd[i]:.2f}")
                    plt.tight_layout()
                    pdf.savefig()  # saves the current figure into the PDF
                    plt.close()

                    PSD.append(psd_ens/ens)
                    PSR.append(P_ens)
                    NF.append(NF_ens)
                    SNR.append(SNR_ens)
                    COV.append(Cov_ens/ens)
                    sd_crt.append(sr.sd_correction(Edrive*tf,E0*tf))
                    bw_crt.append(sr.bandwith_correction(freq, psd_noise_ens/ens)[1])

            sd_crt = np.array(sd_crt)
            bw_crt = np.array(bw_crt)

            Dext_delF = np.array(sr.Dext_delF(sd_crt, bw_crt, alpha, beta, rho, tf, Af))
        return sd_crt, bw_crt, Dext_delF, PSR, NF, SNR, COV
    
    def SR_matric_exp(self, params, raw_data_file, pdf, delw =5 , delw_noise = 70):
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
        # ts = params.get('ts')
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

        fro = 250

        for k in range(len(fsig)): #incase multiple frequency needs to be simulated
            f = fsig[k]
            for j in range(len(Asig)): #incase multiple amplitude levels needs to be simulated
                A = Asig[j]
                tdrive, Vdrive, _ = sr.extract_tEP_data(raw_data_file + f"\\data_sd_{sd[0]}_itr{0}.txt")
                V0 =  Vdrive[fro:]
                ts = (tdrive[1] - tdrive[0])*1e-3
                for i in range(len(sd)):
                    Cov_ens = 0
                    psd_ens = 0
                    psd_noise_ens = 0
                    for e in range(ens):
                        if(i == 0): # this is an adjustment to account for only 1 ensemble for the 0 sd
                            e = 0
                        if((i == 13) & (e == 9)):
                            e = 8  # this is an adjustment 

                        tdrive, Vdrive, P = sr.extract_tEP_data(raw_data_file + f"\\data_sd_{sd[i]}_itr{e}.txt")
                        Vdrive = Vdrive[fro:]
                        P = P[fro:]
                        Cov_ens += np.cov(P,V0)[1,0]    
                        freq, psd = sr.signal_psd(P, f, 1/ts, delw)
                        freq, psd_noise = sr.signal_psd(Vdrive - V0, f, 1/ts, delw)
                        psd_ens += psd
                        psd_noise_ens += psd_noise

                    P_ens, NF_ens, SNR_ens, m, c = sr.power_snr_noisefloor(freq, psd_ens/ens , f,  1/ts, delw, delw_noise)
                    
                    ffro = int(2*len(freq)/ (1/ts)*(f - delw_noise)) + 1
                    fto = int(2*len(freq)/ (1/ts)*(f + delw_noise)) + 1
                    
                    plt.figure(figsize=(8, 5))
                    plt.semilogy(freq[ffro:fto], psd_ens[ffro:fto]/ens)
                    plt.semilogy(freq[ffro:fto], m*freq[ffro:fto]+c)
                    plt.xlabel("Frequency (Hz)")
                    plt.ylabel("PSD")
                    plt.title(f"PSD_{f}hz_{A}V_{sd[i]:.2f}")
                    plt.tight_layout()
                    pdf.savefig()  # saves the current figure into the PDF
                    plt.close()

                    PSD.append(psd_ens/ens)
                    PSR.append(P_ens)
                    NF.append(NF_ens)
                    SNR.append(SNR_ens)
                    COV.append(Cov_ens/ens)
                    sd_crt.append(sr.sd_correction(Vdrive,V0))
                    bw_crt.append(sr.bandwith_correction(freq, psd_noise_ens/ens)[1])

            sd_crt = np.array(sd_crt)
            bw_crt = np.array(bw_crt)

            Dext_delF = sr.Dext_delF(sd_crt, bw_crt, alpha, beta, rho, tf, Af)

        return sd_crt, bw_crt, Dext_delF, PSR, NF, SNR, COV

    
    def SR_matric_plots(self, pdf, sd_crt, Dext_delF, PSR, NF, SNR, COV, sd_opt, Dext_delF_opt):
        # Plotting sd**2
        plt.figure(figsize=(10, 5))
        plt.plot(sd_crt**2, 10*np.log10(np.array(PSR)), "--", color = "black"); plt.title("Power")
        plt.axvline(x=sd_opt**2, color='grey', linestyle='--', label=f'SD = {sd_opt:.3f}')
        plt.xlabel("sd^2")
        plt.ylabel("Power (dB)")
        plt.grid(alpha = 0.5)
        plt.tight_layout()
        pdf.savefig()  # saves the current figure into the PDF
        plt.close()

        plt.figure(figsize=(10, 5))
        plt.plot(sd_crt**2, 10*np.log10(SNR), "s-", color = "black"); plt.title("SNR")
        plt.axvline(x=sd_opt**2, color='grey', linestyle='--', label=f'SD = {sd_opt:.3f}')
        plt.xlabel("sd^2")
        plt.ylabel("SNR (dB/rad/s)")
        plt.grid(alpha = 0.5)
        plt.tight_layout()
        pdf.savefig()  # saves the current figure into the PDF
        plt.close()

        # plt.figure(figsize=(10, 5))
        # plt.plot(sd_crt**2, NF, "o-", color = "black"); plt.title("NF")
        # plt.xlabel("sd^2")
        # plt.ylabel("NF")
        # plt.grid(alpha = 0.5)
        # plt.tight_layout()
        # pdf.savefig()  # saves the current figure into the PDF
        # plt.close()
        
        plt.figure(figsize=(10, 5))
        plt.plot(sd_crt**2, COV, "-*", color = "black"); plt.title("Cross covariance")
        plt.axvline(x=sd_opt**2, color='grey', linestyle='--', label=f'SD = {sd_opt:.3f}')
        plt.xlabel("sd^2")
        plt.ylabel("C.C")
        plt.grid(alpha = 0.5)        
        plt.tight_layout()
        pdf.savefig()  # saves the current figure into the PDF
        plt.close()
        
        #############################################

        # Plotting Dext
        plt.figure(figsize=(10, 5))
        plt.plot(Dext_delF, 10*np.log10(np.array(PSR)), "--", color = "black"); plt.title("Power")
        plt.axvline(x=Dext_delF_opt, color='grey', linestyle='--', label=f'Dext_delF = {Dext_delF_opt:.3f}')
        plt.xlabel("Dext_delF")
        plt.ylabel("Power (dB)")
        plt.grid(alpha = 0.5)
        plt.tight_layout()
        pdf.savefig()  # saves the current figure into the PDF
        plt.close()

        plt.figure(figsize=(10, 5))
        plt.plot(Dext_delF, 10*np.log10(SNR), "s-", color = "black"); plt.title("SNR")
        plt.axvline(x=Dext_delF_opt, color='grey', linestyle='--', label=f'Dext_delF = {Dext_delF_opt:.3f}')
        plt.xlabel("Dext_delF")
        plt.ylabel("SNR (dB/rad/s)")
        plt.grid(alpha = 0.5)
        plt.tight_layout()
        pdf.savefig()  # saves the current figure into the PDF
        plt.close()
        
        plt.figure(figsize=(10, 5))
        plt.plot(Dext_delF, COV, "-*", color = "black"); plt.title("Cross covariance")
        plt.axvline(x=Dext_delF_opt, color='grey', linestyle='--', label=f'Dext_delF = {Dext_delF_opt:.3f}')
        plt.xlabel("Dext_delF")
        plt.ylabel("C.C")
        plt.grid(alpha = 0.5)        
        plt.tight_layout()
        pdf.savefig()  # saves the current figure into the PDF
        plt.close()

    def SR_matric_plots_overlay_exp(self, pdf, sd_crt, Dext_delF, PSR, SNR, COV, sd_crt_exp, Dext_delF_exp, PSR_exp, SNR_exp, COV_exp, sd_opt, Dext_delF_opt):
            # Plotting
            plt.figure(figsize=(12, 5))
            plt.subplot(1,2,1); plt.plot(sd_crt**2, 10*np.log10(np.array(PSR)), "--", color = "black", label = "num"); plt.title("Power")
            plt.axvline(x=sd_opt**2, color='grey', linestyle='--', label=f'SD = {sd_opt:.3f}')
            plt.xlabel("sd^2")
            plt.ylabel("Power (dB)")
            plt.grid(alpha = 0.5)
            plt.tight_layout()
            plt.legend()
            plt.subplot(1,2,2); plt.plot(sd_crt_exp**2, 10*np.log10(np.array(PSR_exp)), "s-", color = "blue", label = "exp")
            plt.axvline(x=sd_opt**2, color='grey', linestyle='--', label=f'SD = {sd_opt:.3f}')
            plt.xlabel("sd^2")
            plt.ylabel("Power (dB)")
            plt.grid(alpha = 0.5)
            plt.tight_layout()
            plt.legend()
            pdf.savefig()  # saves the current figure into the PDF
            plt.close()

            plt.figure(figsize=(12, 5))
            plt.subplot(1,2,1); plt.plot(sd_crt**2, 10*np.log10(SNR), "s-", color = "black", label = "num"); plt.title("SNR")
            plt.axvline(x=sd_opt**2, color='grey', linestyle='--', label=f'SD = {sd_opt:.3f}')
            plt.xlabel("sd^2")
            plt.ylabel("SNR (dB/rad/s)")
            plt.grid(alpha = 0.5)
            plt.tight_layout()
            plt.legend()
            plt.subplot(1,2,2); plt.plot(sd_crt_exp**2, 10*np.log10(SNR_exp), "s", color = "blue", label = "exp")
            plt.axvline(x=sd_opt**2, color='grey', linestyle='--', label=f'SD = {sd_opt:.3f}')
            plt.xlabel("sd^2")
            plt.ylabel("SNR (dB/rad/s)")
            plt.grid(alpha = 0.5)
            plt.tight_layout()
            plt.legend()
            pdf.savefig()  # saves the current figure into the PDF
            plt.close()

            plt.figure(figsize=(12, 5))
            plt.subplot(1,2,1); plt.plot(sd_crt**2, COV, "-*", color = "black", label = "num"); plt.title("Cross covariance")
            plt.axvline(x=sd_opt**2, color='grey', linestyle='--', label=f'SD = {sd_opt:.3f}')
            plt.xlabel("sd^2")
            plt.ylabel("C.C")
            plt.grid(alpha = 0.5)        
            plt.tight_layout()
            plt.legend()
            plt.subplot(1,2,2); plt.plot(sd_crt_exp**2, COV_exp, "-*", color = "blue", label = "exp"); plt.title("Cross covariance")
            plt.axvline(x=sd_opt**2, color='grey', linestyle='--', label=f'SD = {sd_opt:.3f}')
            plt.xlabel("sd^2")
            plt.ylabel("C.C")
            plt.grid(alpha = 0.5)        
            plt.tight_layout()
            plt.legend()
            pdf.savefig()  # saves the current figure into the PDF
            plt.close()

            #########################

            # Plotting with Dext_delF
            plt.figure(figsize=(12, 5))
            plt.subplot(1,2,1); plt.plot(Dext_delF, 10*np.log10(np.array(PSR)), "--", color = "black", label = "num"); plt.title("Power")
            plt.axvline(x=Dext_delF_opt, color='grey', linestyle='--', label=f'Dext_delF = {Dext_delF_opt:.3f}')
            plt.xlabel("Dext_delF")
            plt.ylabel("Power (dB)")
            plt.grid(alpha = 0.5)
            plt.tight_layout()
            plt.legend()
            plt.subplot(1,2,2); plt.plot(Dext_delF_exp, 10*np.log10(np.array(PSR_exp)), "s-", color = "blue", label = "exp")
            plt.axvline(x=Dext_delF_opt, color='grey', linestyle='--', label=f'Dext_delF = {Dext_delF_opt:.3f}')
            plt.xlabel("Dext_delF")
            plt.ylabel("Power (dB)")
            plt.grid(alpha = 0.5)
            plt.tight_layout()
            plt.legend()
            pdf.savefig()  # saves the current figure into the PDF
            plt.close()

            plt.figure(figsize=(12, 5))
            plt.subplot(1,2,1); plt.plot(Dext_delF, 10*np.log10(SNR), "s-", color = "black", label = "num"); plt.title("SNR")
            plt.axvline(x=Dext_delF_opt, color='grey', linestyle='--', label=f'Dext_delF = {Dext_delF_opt:.3f}')
            plt.xlabel("Dext_delF")
            plt.ylabel("SNR (dB/rad/s)")
            plt.grid(alpha = 0.5)
            plt.tight_layout()
            plt.legend()
            plt.subplot(1,2,2); plt.plot(Dext_delF_exp, 10*np.log10(SNR_exp), "s", color = "blue", label = "exp")
            plt.axvline(x=Dext_delF_opt, color='grey', linestyle='--', label=f'Dext_delF = {Dext_delF_opt:.3f}')
            plt.xlabel("Dext_delF")
            plt.ylabel("SNR (dB/rad/s)")
            plt.grid(alpha = 0.5)
            plt.tight_layout()
            plt.legend()
            pdf.savefig()  # saves the current figure into the PDF
            plt.close()

            plt.figure(figsize=(12, 5))
            plt.subplot(1,2,1); plt.plot(Dext_delF, COV, "-*", color = "black", label = "num"); plt.title("Cross covariance")
            plt.axvline(x=Dext_delF_opt, color='grey', linestyle='--', label=f'Dext_delF = {Dext_delF_opt:.3f}')
            plt.xlabel("Dext_delF")
            plt.ylabel("C.C")
            plt.grid(alpha = 0.5)        
            plt.tight_layout()
            plt.legend()
            plt.subplot(1,2,2); plt.plot(Dext_delF_exp, COV_exp, "-*", color = "blue", label = "exp"); plt.title("Cross covariance")
            plt.axvline(x=Dext_delF_opt, color='grey', linestyle='--', label=f'Dext_delF = {Dext_delF_opt:.3f}')
            plt.xlabel("Dext_delF")
            plt.ylabel("C.C")
            plt.grid(alpha = 0.5)        
            plt.tight_layout()
            plt.legend()
            pdf.savefig()  # saves the current figure into the PDF
            plt.close()

    def Ktime_opt_sd(self, params, sd_range, fsig, pdf):
        alpha = params.get('alpha')
        beta = params.get('beta')
        rho = params.get('rho')
        tf = params.get('tf')
        Af = params.get('Af')
        bw = params.get('bw')

        Ktime = sr.calculate_kramer_curve(sd_range, 0, alpha, beta, rho, bw, tf, Af)
        # Target Ktime in ms
        target_time = 1/(2*fsig)

        interp_fn = sc.interpolate.interp1d(Ktime, sd_range, bounds_error=False, fill_value="extrapolate")
        sd_opt = interp_fn(target_time)
        Dext_delF_opt = sr.Dext_delF(sd_opt, bw, alpha, beta, rho, tf, Af)[0]

        # Plot the curve
        plt.plot(sd_range, Ktime*1e3, linewidth = 2.5, color = "blue", label='Kramers Curve')
        # Draw horizontal line at target time
        plt.axhline(y=target_time*1e3, color='gray', linestyle='--', label=f'Ktime = {target_time*1e3:.2f} ms')
        # Draw vertical line at corresponding SD
        plt.axvline(x=sd_opt, color='black', linestyle='--', label=f'SD = {sd_opt:.3f}')
        # Axis labels and title
        plt.xlabel('Standard Deviation (SD)')
        plt.ylabel('Ktime (ms)')
        plt.title('Ktime vs SD at 0V')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        pdf.savefig()  # saves the current figure into the PDF
        plt.close()

        return target_time, sd_opt, Dext_delF_opt