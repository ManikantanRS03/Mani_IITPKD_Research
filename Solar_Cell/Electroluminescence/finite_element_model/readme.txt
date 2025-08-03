. parameter_gen_script ->
  is to generate the parameters such as front resistance network, back resistance network, rshunt, rseries etc

. solar_cell_fin_El ->
  simulation of EL imaging to generate EL images with the given parameters.

. solar_cell_fin_model ->
  take input from the file parameters.pkl, all the grid parameters, simulation parameters, probing points etc.
  put output in log files this will create a log file for each parameters.pkl file.. (can be later used for verification)

. solar_cell_grid_3X3 ->
  an example code to demonstrate the finite element modeling

. cell_img_to_gris ->
  code to find the identify the busbar and fingers from an actual solar cell image