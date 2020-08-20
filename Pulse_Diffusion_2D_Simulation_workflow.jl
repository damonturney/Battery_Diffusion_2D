# Work Flow
ss = Diffusion_2D.create_system_state();
ss,sim_data = Diffusion_2D.run_simulation(ss , 0:Int64(1E7), 1E4,  -0.01, 20, -0.01 , 10, 1E-6);  #0:Int64(1E7), 1E5, -0.01, 20, -0.01 , 10,  1E-6); # args: system_state, iterations, saved_iteration_spacing, electrode_voltage1, electrode_voltage1_ontime, electrode_voltage2, electrode_voltage2_ontime, dt 
#Diffusion_2D.plot_results_time_slice(ss,sim_data,5)
Diffusion_2D.make_movie(ss,sim_data,505)
#using FileIO
#cv_data = get(load("produced_data/20200803102000/20200803102016_dictionary_results.jld2"), "cyclicV_data",0);
#ss =      get(load("produced_data/20200803102000/20200803102016_dictionary_results.jld2"), "system_state"          ,0);
#ss =      get(load(data_file_name), "system_state",0);
#cv_data = get(load(data_file_name), "cyclicV_data",0);
#Diffusion_2D.plot_results_time_series(ss,cv_data)