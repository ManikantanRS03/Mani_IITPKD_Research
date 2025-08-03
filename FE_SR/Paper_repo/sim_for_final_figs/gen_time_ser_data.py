import os
import time
import shutil
import numpy as np
import SR_lib as sr 

class Simulation():
    def __init__(self):
        pass

    def Run_one_sim(self, sd, e, A, f, ts, bw, n, dsf, Dint, alpha, beta, rho, tf, Pr, output_dir, log_file):
        sr.log_print(f"\nsd = {sd:.2f} running ...",log_file)
        start_time = time.time()
        
        # Data generation
        ncycle = n #this method is to generate data without exceeding the size limit
        
        t_down = []
        Edrive_down = []
        P_down = []
        
        tlast = 0
        Elast = 0
        Plast = -Pr
        offset = 0

        while(ncycle>0):
            sr.log_print(f"{ncycle} cycles remaining...", log_file)
            if(ncycle>50):
                t, Vdrive = sr.sin_wave(A, f, 50, ts) 
                t = np.concatenate(([tlast], t[1:]+tlast))
                noise = sr.gen_noise(t, sd, 1/bw)
                Edrive = (Vdrive + noise + offset) / tf
                Edrive = np.concatenate(([Elast], Edrive[1:]))
                P = sr.TDGL(t, Plast, Edrive, Dint, alpha, beta, rho)
                #updating the last elementens
                tlast = t[-1]
                Elast = Edrive[-1]
                Plast = P[-1]
                # Downsampling the data
                t_down_s = t[::dsf]
                Edrive_down_s = Edrive[::dsf]
                P_down_s = P[::dsf]
                #updating the ncycle remaining
                ncycle = ncycle - 50
            else:
                t, Vdrive = sr.sin_wave(A, f, ncycle, ts)
                t = np.concatenate(([tlast], t[1:]+tlast))
                noise = sr.gen_noise(t, sd, 1/bw)
                Edrive = (Vdrive + noise) / tf
                Edrive = np.concatenate(([Elast], Edrive[1:]))
                P = sr.TDGL(t, Plast, Edrive, Dint, alpha, beta, rho)
                #updating the last elementens
                tlast = t[-1]
                Elast = Edrive[-1]
                Plast = P[-1]
                # Downsampling the data
                t_down_s = t[::dsf]
                Edrive_down_s = Edrive[::dsf]
                P_down_s = P[::dsf]
                #updating the ncycle remaining
                ncycle = 0

            t_down = np.concatenate((t_down, t_down_s[1:]))
            Edrive_down = np.concatenate((Edrive_down, Edrive_down_s[1:]))
            P_down = np.concatenate((P_down, P_down_s[1:]))

        # Save data to text file
        data_dir = f'raw_data'
        os.makedirs(os.path.join(output_dir, data_dir), exist_ok=True)
        file_name = f'PE_{f}hz_{A}V_sd{sd:.2f}_ens{e}.txt'
        write_file = os.path.join(os.path.join(output_dir, data_dir), file_name)
        np.savetxt(write_file, np.column_stack((t_down, Edrive_down, P_down)), header='t Edrive P')
        
        # Time taken for this iteration
        end_time = time.time()
        elapsed_time = end_time - start_time
        sr.log_print(f"Ensemble {e}, sd {sd} completed in {elapsed_time:.2f} seconds.", log_file)


    def SR_sim_with_noise(self, params, output_dir):
        log_file = os.path.join(output_dir, 'log.txt')
        try:
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
        except:
            print("parameters missing simulation is killed, correct the input file")

        sr.log_print(f"sim_{id} folder created",log_file)
        sr.log_print("The simulation parameters are:",log_file)
        sr.log_print(str(params),log_file)
        
        kb = 1.380649e-23
        Dint = (kb*T)/(tf*Af)

        for k in range(len(fsig)): #incase multiple frequency needs to be simulated
            f = fsig[k]
            for j in range(len(Asig)): #incase multiple amplitude levels needs to be simulated
                A = Asig[j]
                #generate one supper threshold signal for cross cov calculation
                Vspr = Vc*1.5
                self.Run_one_sim(0, 0, Vspr, f, ts, bw, n, dsf, Dint, alpha, beta, rho, tf, Pr, output_dir, log_file)
                # Generating data
                sr.log_print(f"\nA = {A:.2f}, f = {f:.2f} running ...",log_file)
                for i in range(len(sd)):
                    for e in range(ens):
                        self.Run_one_sim(sd[i], e, A, f, ts, bw, n, dsf, Dint, alpha, beta, rho, tf, Pr, output_dir, log_file)
        sr.log_print("\nData generation done!", log_file)

    def SR_sim_weiner_process(self, params, output_dir):
        log_file = os.path.join(output_dir, 'log.txt')
        try:
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
        except:
            print("parameters missing simulation is killed, correct the input file")

        sr.log_print(f"sim_{id} folder created",log_file)
        sr.log_print("The simulation parameters are:",log_file)
        sr.log_print(str(params),log_file)

        kb = 1.380649e-23

        for k in range(len(fsig)): #incase multiple frequency needs to be simulated
            f = fsig[k]
            for j in range(len(Asig)): #incase multiple amplitude levels needs to be simulated
                A = Asig[j]
                # Generating data
                sr.log_print(f"\nA = {A:.2f}, f = {f:.2f} running ...",log_file)
                for i in range(len(sd)):
                    sr.log_print(f"\nsd = {sd[i]:.2f} running ...",log_file)
                    Dext = (i**2)/(2*bw*rho*tf**2)
                    Dint = (kb*T)/(tf*Af)

                    for e in range(ens):
                        start_time = time.time()
                        sr.log_print(f"\nens = {e:.2f} running ...",log_file)
                        # Data generation
                        ncycle = n #this method is to generate data without exceeding the size limit
                        
                        t_down = []
                        Edrive_down = []
                        P_down = []
                        
                        tlast = 0
                        Elast = 0
                        Plast = -Pr
                        offset = 0

                        while(ncycle>0):
                            sr.log_print(f"{ncycle} cycles remaining...", log_file)
                            if(ncycle>50):
                                t, Vdrive = sr.sin_wave(A, f, 50, ts) 
                                t = np.concatenate(([tlast], t[1:]+tlast))
                                Esig = (Vdrive + offset) / tf
                                P, Edrive = sr.TDGL_WienerProcess(t, Plast, Esig, Dint, Dext, alpha, beta, rho, bw)
                                Edrive = np.concatenate(([Elast], Edrive[1:]))
                                #updating the last elementens
                                tlast = t[-1]
                                Elast = Edrive[-1]
                                Plast = P[-1]
                                # Downsampling the data
                                t_down_s = t[::dsf]
                                Edrive_down_s = Edrive[::dsf]
                                P_down_s = P[::dsf]
                                #updating the ncycle remaining
                                ncycle = ncycle - 50
                            else:
                                t, Vdrive = sr.sin_wave(A, f, ncycle, ts) 
                                t = np.concatenate(([tlast], t[1:]+tlast))
                                Esig = (Vdrive + offset) / tf
                                P, Edrive = sr.TDGL_WienerProcess(t, Plast, Esig, Dint, Dext, alpha, beta, rho, bw)
                                Edrive = np.concatenate(([Elast], Edrive[1:]))
                                #updating the last elementens
                                tlast = t[-1]
                                Elast = Edrive[-1]
                                Plast = P[-1]
                                # Downsampling the data
                                t_down_s = t[::dsf]
                                Edrive_down_s = Edrive[::dsf]
                                P_down_s = P[::dsf]
                                #updating the ncycle remaining
                                ncycle = 0

                            t_down = np.concatenate((t_down, t_down_s[1:]))
                            Edrive_down = np.concatenate((Edrive_down, Edrive_down_s[1:]))
                            P_down = np.concatenate((P_down, P_down_s[1:]))

                        # Save data to text file
                        data_dir = f'raw_data'
                        os.makedirs(os.path.join(output_dir, data_dir), exist_ok=True)
                        file_name = f'PE_{f}hz_{A}V_sd{sd[i]:.2f}_ens{e}.txt'
                        write_file = os.path.join(os.path.join(output_dir, data_dir), file_name)
                        np.savetxt(write_file, np.column_stack((t_down, Edrive_down, P_down)), header='t Edrive P')
                        
                        # Time taken for this iteration
                        end_time = time.time()
                        elapsed_time = end_time - start_time
                        sr.log_print(f"Ensemble {e}, Iteration {i + 1} completed in {elapsed_time:.2f} seconds.", log_file)
        
        sr.log_print("\nData generation done!", log_file)

    def SR_sim_with_FEdrive(self, params, input_dir, output_dir):
        log_file = os.path.join(output_dir, 'log.txt')
        try:
            Vc = params.get('Vc')
            Pr = params.get('Pr')
            alpha = params.get('alpha')
            beta = params.get('beta')
            rho = params.get('rho')
            tf = params.get('tf')
            Af = params.get('Af')
            T = params.get('T')
            ts = params.get('ts')
            sd = params.get('sd')
            ens = params.get('ens')
            dsf = params.get('dsf')
            fsig = params.get('f')
            Asig = params.get('A')

        except:
            print("parameters missing simulation is killed, correct the input file")

        sr.log_print(f"sim_{id} folder created",log_file)
        sr.log_print("The simulation parameters are:",log_file)
        sr.log_print(str(params),log_file)
        Dint = 0

        for k in range(len(fsig)): #incase multiple frequency needs to be simulated
            f = fsig[k]
            for j in range(len(Asig)): #incase multiple amplitude levels needs to be simulated
                A = Asig[j]
                # Generating data
                sr.log_print(f"\nA = {A:.2f}, f = {f:.2f} running ...",log_file)
                for i in range(len(sd)):
                    for e in range(ens):
                        start_time = time.time()
                        if(i == 0): # this is an adjustment to account for only 1 ensemble for the 0 sd
                            e = 0
                        tdrive, Vdrive, Pexp = sr.extract_tEP_data(input_dir + f"\\data_sd_{sd[i]}_itr{e}.txt")
                        tdrive = tdrive*1e-3 # converting ms to s
                        t, Vsig = sr.upsample(tdrive, Vdrive, ts)
                        Esig = (Vsig) / tf
                        P = sr.TDGL(t, Pr, Esig, Dint, alpha, beta, rho)
                        # Downsampling the data
                        t_down = t[::dsf]
                        Esig_down = Esig[::dsf]
                        P_down = P[::dsf]
                        # Save data to text file
                        data_dir = f'raw_data'
                        os.makedirs(os.path.join(output_dir, data_dir), exist_ok=True)
                        file_name = f'PE_{f}hz_{A}V_sd{sd[i]:.2f}_ens{e}.txt'
                        write_file = os.path.join(os.path.join(output_dir, data_dir), file_name)

                        np.savetxt(write_file, np.column_stack((t_down, Esig_down, P_down)), header='t Edrive P')
                        # Time taken for this iteration
                        end_time = time.time()
                        elapsed_time = end_time - start_time
                        sr.log_print(f"Ensemble {e}, Iteration {i + 1} completed in {elapsed_time:.2f} seconds.", log_file)

            sr.log_print("\nData generation done!", log_file)