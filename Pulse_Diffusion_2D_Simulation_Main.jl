"""
Started on Mon Aug 08 2020

2D Diffusion to a spike

To Do:


@author: damon
"""

module Diffusion_2D

    using Printf
    using FileIO
    using Dates

    if Sys.isapple()
        cd("/Users/damon/Desktop/BACKED_UP/WorkFiles/ReportsPublicationsPatents/20200308_Patent_GravityFlat_PulseCharged_Zinc_Electrodes/Science_Journal_Publication/20200820_Pulse_Diffusion_2D_Simulation//")
    else
        cd("/Users/damon/Desktop/BACKED_UP/WorkFiles/ReportsPublicationsPatents/20200308_Patent_GravityFlat_PulseCharged_Zinc_Electrodes/Science_Journal_Publication/20200820_Pulse_Diffusion_2D_Simulation//")
    end
    
    # Note: All the expert software engineers told me to not define anything in the 
    #       global scope of the modules. They said I should define all variables in
    #       the functions. Something about mutable global state being evil.

    #include("Battery_Simulation_1D_generally_useful_functions.jl")

    ############################################################
    #### Functions to create the initial system state and functions for V_anode and V_cathode
    include("Pulse_Diffusion_2D_create_system_state.jl")

    #function load_old_battery_state(path_to_old_output_file)
    #    println("Loading system state from the results of a previous simulation.")
    #    ss=get(load(path_to_old_output_file*"hi"), "system_state", 0)
    #end
    ############################################################


    ############################################################
    ###### A Function for the Cyclic Voltametry Operator 
    include("Pulse_Diffusion_2D_Simulation_run_simulation.jl")
    ############################################################


    ############################################################
    ###### A Function for Plotting Cyclic Voltametry Data
    include("Pulse_Diffusion_2D_Simulation_Plots.jl")
    ############################################################


end  

# Work Flow
ss = Diffusion_2D.create_system_state();
ss,sim_data = Diffusion_2D.run_simulation(ss , 0:Int64(1E7), 2E4, -0.01, 5E-7); # args: system_state, iterations, saved_iteration_spacing,  voltage_timeseries ,  dt 
plot_results_time_slice(ss,sim_data,104)
#using FileIO
#cv_data = get(load("produced_data/20200803102000/20200803102016_dictionary_results.jld2"), "cyclicV_data",0);
#ss =      get(load("produced_data/20200803102000/20200803102016_dictionary_results.jld2"), "system_state"          ,0);
#ss =      get(load(data_file_name), "system_state",0);
#cv_data = get(load(data_file_name), "cyclicV_data",0);
#Diffusion_2D.plot_results_time_series(ss,cv_data)

