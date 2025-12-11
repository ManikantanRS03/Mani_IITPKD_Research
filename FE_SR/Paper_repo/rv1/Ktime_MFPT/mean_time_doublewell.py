from multiprocessing import Pool
import functools
import numpy as np
import scipy as sc
import json

def gen_noise(t, sd, bw, dt):
    t_coarse = np.arange(0, t, 1/bw)
    coarse_noise = np.random.normal(0, sd, size=len(t_coarse))
    wgn = np.interp(np.arange(0,t,dt), t_coarse, coarse_noise, left=coarse_noise[0], right=coarse_noise[-1])
    return wgn

def sim_ktime(t, noise, P_init, alpha, beta, rho, tf, dt, vbias):
    P = [P_init]
    for i in range(1,int(t/dt)):
        df_dp = 2*alpha*P[i-1] + 4*beta*P[i-1]**3 - (vbias + noise[i-1])/tf
        P.append(P[i-1] - dt/rho*df_dp)
        if(P[i]> 0.1):
            return True, i*dt, P[i], P
    return False, i*dt, P[i], P

def Ktime_anl(bias, alpha, beta, rho, Dext, tf):
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

def run_ensemble(e, D, sd, t, bw, dt, alpha, beta, rho, tf, vbias, ens):
    """Helper function to run a single ensemble simulation"""
    stat = False
    kt = 0
    P_in = -0.2
    while not stat:
        noise = gen_noise(t, sd, bw, dt)
        stat, kt_ens, P_in, _ = sim_ktime(t, noise, P_in, alpha, beta, rho, tf, dt, vbias)
        kt += kt_ens
    
    if e % 100 == 0:
        print(f"Ensemble {e}/{ens} done")
    return kt


if __name__ == '__main__':
    ######################################################
    alpha = -5.38e7
    beta = 1.04e9
    rho = 370
    tf = 255e-9
    bw = 1e4
    vrms = np.linspace(0.4, 0.65, 10)
    Dext = vrms**2/(2*rho*bw*tf**2)

    # Dext_DelF = np.linspace(0.7, 1.2, 20)
    # Dext = np.arange (4e4, 9e4, 0.5e4)
    
    # simulation parameters
    dt = 10e-9
    t = 5e-4
    vbias = 1.2
    ens = 500
    delF = del_U(0, alpha, beta, tf)[0]
    # Dext = Dext_DelF * delF
    ######################################################

    # Parallel execution
    kt_D = {}
    for D in Dext:
        print(f"Simulating for D/F = {D/delF:.2f} ...")
        sd = np.sqrt(2*D*bw*rho)*tf
        
        # Create a partial function with fixed parameters
        run_func = functools.partial(
            run_ensemble, 
            D=D, sd=sd, t=t, bw=bw, dt=dt, 
            alpha=alpha, beta=beta, rho=rho, tf=tf, vbias=vbias, ens=ens
        )
        
        # Run ensembles in parallel
        with Pool() as pool:
            kt_list = pool.map(run_func, range(ens+1))
        
        kt_D[D] = kt_list
        print(f"D/F = {D/delF:.2f}: kt = {np.mean(kt_list)*1e3:.2f} ms")

    ######################################################
    # Collect all parameters
    parameters = {
        "alpha": alpha,
        "beta": beta,
        "rho": rho,
        "tf": tf,
        "dt": dt,
        "t": t,
        "vbias": vbias,
        "ens": ens,
        "bw": bw,
        "Dext": [float(D) for D in Dext]
    }

    # Convert numpy arrays/floats to Python native types for JSON compatibility
    kt_D_serializable = {str(D): [float(kt) for kt in kt_list] for D, kt_list in kt_D.items()}

    # Combine parameters and results
    output_data = {
        "parameters": parameters,
        "results": kt_D_serializable
    }

    # Save to JSON file
    with open(f'kt_D_results_{vbias}_B_exp.json', 'w') as f:
        json.dump(output_data, f, indent=2)

    print("Data and parameters saved to kt_D_results.json")