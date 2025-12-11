import numpy as np
import scipy as sc
import re

def read_simulation_parameters(file_path):
    variables = {}
    with open(file_path, 'r') as file:
        for line in file:
            # Ignore comments and empty lines
            line = line.strip()
            if line.startswith('#') or not line:
                continue

            # Split the line into variable and value
            if '=' in line:
                var, value = line.split('=', 1)
                var = var.strip()
                value = value.strip()

                # Evaluate numeric and list/tuple values
                try:
                    evaluated_value = eval(value)
                    
                    # If tuple, convert to list
                    if isinstance(evaluated_value, tuple):
                        variables[var] = list(evaluated_value)
                    
                    # If list with 3 elements, use linspace
                    elif isinstance(evaluated_value, list) and len(evaluated_value) == 3:
                        start, end, points = evaluated_value
                        variables[var] = np.linspace(start, end, int(points))
                    
                    # For other types, keep as is
                    else:
                        variables[var] = evaluated_value
                
                except (NameError, SyntaxError):
                    variables[var] = value
    return variables

def read_tEP_data(file_path): # Reading the saved text file
    data = np.loadtxt(file_path)
    t = data[:, 0]
    Edrive = data[:, 1]
    P = data[:, 2]
    return t, Edrive, P

def extract_tEP_data(file_path): #read from text file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    data_started = False
    time_ms = []
    drive_voltage = []
    polarization = []

    for line in lines:
        # Detect start of data by finding the column header
        if "Point" in line and "Time" in line and "DRIVE Voltage" in line:
            data_started = True
            continue

        if data_started:
            # Split on tabs or multiple spaces
            values = re.split(r'\s+|\t+', line.strip())
            if len(values) >= 3 and values[0].isdigit():
                try:
                    
                    time_ms.append(float(values[1]))
                    drive_voltage.append(float(values[2]))
                    polarization.append(float(values[3])*1e-2)
                    
                except ValueError:
                    continue  # Skip malformed lines

    return np.array(time_ms), np.array(drive_voltage), np.array(polarization)

def log_print(message, log_file): # Define a function to print and log simultaneously
    print(message)
    with open(log_file, 'a') as log:
        log.write(message + '\n')

def TDGL(t, Pini, E, Dint, alpha, beta, rho): #solve TDGL equation
    P_fit = [Pini]
    dt = t[1] - t[0]
    for i in range(1,len(t)):
        df_dp = 2*alpha*P_fit[i-1] + 4*beta*P_fit[i-1]**3 - E[i]
        # noise_intr = (np.sqrt(2*Dint*rho)) * np.random.normal(0, np.sqrt(dt))/dt
        P_fit.append(P_fit[i-1] - dt/rho*df_dp) # + dt/rho*noise_intr)
    return P_fit

def TDGL_WienerProcess(t, Pini, E, Dint, Dext, alpha, beta, rho, bw): #solve TDGL equation
    P_fit = [Pini]
    E_drive = [E[0]]
    dt = t[1] - t[0]
    for i in range(1,len(t)):
        df_dp = 2*alpha*P_fit[i-1] + 4*beta*P_fit[i-1]**3 - E[i]
        noise_disc = (np.sqrt(2*Dext*rho) + np.sqrt(2*Dint*rho)) * np.random.normal(0, np.sqrt(dt))/dt
        P_fit.append(P_fit[i-1] - dt/rho*df_dp + dt/rho*noise_disc)
        E_drive.append(E[i]+noise_disc)
    return P_fit, E_drive

def sin_wave(A, frequency, n, dt): #generating sine wave (Amplitude, frequency, no cycles, timestep)
    T = 1 / frequency  # Period of the wave
    t = np.arange(0, n*T, dt)  # Time values
    y =  A*np.sin(2*np.pi*t*frequency)
    return t, y

def gen_noise(t, sd, coarse_step): #noise generation code, generate bandlimited noise
    t_start, t_end = t[0], t[-1]
    t_coarse = np.arange(t_start, t_end, coarse_step)
    # Generate WGN at coarse sampling rate
    coarse_noise = np.random.normal(0, sd, size=len(t_coarse))
    # Sample-and-hold interpolation (nearest neighbor)
    wgn = np.interp(t, t_coarse, coarse_noise, left=coarse_noise[0], right=coarse_noise[-1])
    return wgn

def upsample(t_drive, Vdrive, ts):
    # New time vector for oversampling
    t = np.arange(0, t_drive[-1] + ts, ts)
    # Interpolate
    interp_func = sc.interpolate.interp1d(t_drive, Vdrive, kind='linear', fill_value="extrapolate")
    Vsig = interp_func(t)
    return t, Vsig


def hysteresis(Vc, f, ts, alpha, beta, rho, tf): #to plot the hysteresis of the given sample for the given frequency
    t,V = sin_wave(2*Vc, f, 2, ts)
    E = V/tf
    P = TDGL(t, 0, E, alpha, beta, rho)
    return E,P

def mylsq(x, y):
    e1 = np.sum(x**2)
    e2 = np.sum(x)
    e3 = e2
    e4 = len(x)
    A = np.array([[e1, e2], [e3, e4]])

    e5 = np.sum(x * y)
    e6 = np.sum(y)
    B = np.array([e5, e6])
    # Solve the linear equation system
    X = np.linalg.solve(A, B)
    return X

def signal_psd(P, f_sig, fs, delw):
    freq, psd = sc.signal.periodogram(P, fs, window = "boxcar")
    #index to find power
    ffro = int(2*len(freq)/fs*(f_sig - delw)) + 1
    fto = int(2*len(freq)/fs*(f_sig + delw)) + 1
    return freq, psd

def power_snr_noisefloor(freq, psd , f_sig, fs, delw, delw_noise):    
    #index to find power
    ffro = int(2*len(freq)/fs*(f_sig - delw)) + 1
    fto = int(2*len(freq)/fs*(f_sig + delw)) + 1
    #index to find noise
    nfro1 = int(2*len(freq)/fs*(f_sig - delw_noise))+1
    nto1 = int(2*len(freq)/fs*(f_sig + delw_noise))+1
    m, c = mylsq(np.concatenate((freq[nfro1:ffro],freq[fto:nto1])), np.concatenate((psd[nfro1:ffro], psd[fto:nto1])))
    noise_floor = m*f_sig + c
    P_sig = np.trapz(psd[ffro:fto], freq[ffro:fto])
    Snr = P_sig/(noise_floor*2*delw)
    return P_sig, noise_floor, Snr,  m, c
    

def power_snr_noisefloor_iisc(freq, psd, f_sig, fs, delw, delw_noise):    
    #index to find power
    ffro = int(2*len(freq)/fs*(f_sig - delw)) + 1
    fto = int(2*len(freq)/fs*(f_sig + delw)) + 1
    #index to find noise
    nfro1 = int(2*len(freq)/fs*(f_sig - delw_noise))+1
    nto1 = int(2*len(freq)/fs*(f_sig + delw_noise))+1
    P_sig = np.trapz(psd[ffro:fto], freq[ffro:fto])
    P_noise = np.trapz(psd[nfro1:nto1], freq[nfro1:nto1]) - P_sig
    Snr = P_sig/P_noise
    return P_sig, P_noise, Snr
    
##################### correctrion #########################

def sd_correction(Vdrive, V0):
    return np.std(Vdrive - V0)

def bandwith_correction(freq, psd_noise, mov_avg_len = 100):
    psd_median = []
    for j in range(len(freq) - mov_avg_len):
        psd_median.append(np.median(psd_noise[j:j+mov_avg_len]))
    # code to find the 3dB point
    noise_floor = psd_median[0]
    indx_3db = np.where(psd_median < (noise_floor/2)) #find where the power is less than half
    if(len(indx_3db[0]) == 0):
         bw_crt= 2*freq[-1]
    else:
        bw_crt = 2*freq[indx_3db[0][0]]
    return [noise_floor, bw_crt]

def Dext_delF(sd, bw, alpha, beta, rho, tf, Af):
    return sd**2/(2*tf**2*rho*bw)/del_U(0, alpha, beta, tf)
    
##################### Ktime Calac #########################

def calculate_kramer_curve(sigma, bias, alpha, beta, rho, del_f, tf, Af):   
    R_F = rho*tf/Af # Resistance of the sample
    Dext = sigma**2/(2*R_F*del_f*tf*Af)

    P_A = sc.optimize.fsolve(dU_dx,-1,args=(bias, alpha, beta, tf))
    P_C = sc.optimize.fsolve(dU_dx,0,args=(bias, alpha, beta, tf))
    ddFA = 2*alpha + 12*beta*P_A**2
    ddFC = 2*alpha + 12*beta*P_C**2
    delta_U = abs(U(P_A, bias, alpha, beta, tf) - U(P_C, bias, alpha, beta, tf))

    tau_k = 2*np.pi*rho/(np.sqrt(np.abs(ddFA*ddFC)))*np.exp(delta_U/Dext)
    return tau_k

def del_U(bias, alpha, beta, tf):
    P_A = sc.optimize.fsolve(dU_dx,-1,args=(bias, alpha, beta, tf))
    P_C = sc.optimize.fsolve(dU_dx,0,args=(bias, alpha, beta, tf))
    ddFA = 2*alpha + 12*beta*P_A**2
    ddFC = 2*alpha + 12*beta*P_C**2
    return abs(U(P_A, bias, alpha, beta, tf) - U(P_C, bias, alpha, beta, tf))

def U(P,bias, alpha, beta, tf):
    return (alpha*P**2 +beta*P**4-P*bias/tf)

def dU_dx(P,bias, alpha, beta, tf):
    return 2*alpha*P +4*beta*P**3 - bias/tf

##################### landau coeff #########################

def find_alpha_beta(Vc, Pr, tf):
    Ec= Vc/tf
    alpha = -3*np.sqrt(3)*Ec/4/Pr
    beta = 3*np.sqrt(3)*Ec/8/(Pr**3)
    print(f"VC = {Vc} V")
    print(f"Pr = {Pr*10} uC/cm^2")
    print(f"Landau coeff, alpha = {alpha:.2e} mF^-1, beta = {beta:.2e} m^5F^-1C^-1")
    return alpha, beta