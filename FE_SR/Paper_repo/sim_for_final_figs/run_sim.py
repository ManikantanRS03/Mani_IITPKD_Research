import SR_lib as sr 
from gen_time_ser_data import *
from analysis_data import *
import os

# ENABLE FLAGS
GEN_TDGL = 1
ANL = 0

## Read the Parse file
params_path = 'parameter_parse.txt'  # Parse input file
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

## Create the simulation directory
id = params.get('id')
output_dir = f'sim_{id}'
os.makedirs(output_dir, exist_ok=True)

## Run the data generation simulation (simulating TDGL with generated signal)
if(GEN_TDGL == 1):
    run = Simulation() #creating the object of the class simulation
    run.SR_sim_with_noise(params, output_dir) #run simulation type
    
    # copy the parameters to the destination
    destination = os.path.join(output_dir, 'parameters.txt')
    shutil.copy2(params_path, destination)


## Run the data analyais
if(ANL == 1):
    anl = Analysis() #creating the object of the class Analysis
    
    # raww data file
    raw_data_file = os.path.join(output_dir, "raw_data")
  
    # #generating timeseries plot
    # output_pdf_path = os.path.join(output_dir, "analysis", "all_PE_plots.pdf")
    # os.makedirs(os.path.dirname(output_pdf_path), exist_ok=True)
    
    # ens = 0
    # with PdfPages(output_pdf_path) as pdf:
    #     for std in sd:
    #         for f in fsig:
    #             for A in Asig:
    #                     anl.t_V_P_plots(params, std, f, A, ens, raw_data_file, pdf)
    # print("Time series plot generation done!")

    #generating PSD plot
    output_pdf_path = os.path.join(output_dir, "analysis", "all_PSD_plots.pdf")
    os.makedirs(os.path.dirname(output_pdf_path), exist_ok=True)
    with PdfPages(output_pdf_path) as pdf:
        sd_crt, bw_crt, Dext_delF, PSR, NF, SNR, COV = anl.SR_matrics(params, raw_data_file , pdf, delw =2 , delw_noise = 70)
    print("SR matric calculation done! PSD plot generation done!")

    #generate Ktime plot
    output_pdf_path = os.path.join(output_dir, "analysis", "Ktime_opt_sd.pdf")
    os.makedirs(os.path.dirname(output_pdf_path), exist_ok=True)
    with PdfPages(output_pdf_path) as pdf:
         sd_range = np.linspace(0.7, 0.9, 100)
         fsig = 100
         target_time, sd_opt, Dext_delF_opt = anl.Ktime_opt_sd(params, sd_range, fsig, pdf)
    print("Ktime calculation and plot generation done!")

    #generating SR matric plots
    output_pdf_path = os.path.join(output_dir, "analysis", "SR_matric.pdf")
    os.makedirs(os.path.dirname(output_pdf_path), exist_ok=True)
    with PdfPages(output_pdf_path) as pdf:
         anl.SR_matric_plots(pdf, sd_crt, Dext_delF, PSR, NF, SNR, COV, sd_opt, Dext_delF_opt)
    print("SR matric plot generation done!")
    

