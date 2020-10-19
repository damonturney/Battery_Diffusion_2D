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
   import Statistics.mean

   if Sys.isapple()
      cd("/Users/damon/Desktop/BACKED_UP/WorkFiles/ReportsPublicationsPatents/20200308_Patent_GravityFlat_PulseCharged_Zinc_Electrodes/Science_Journal_Publication/20200820_Pulse_Diffusion_2D_Simulation/")
   else
      cd("/home/damon/Simulations/20200715_Binary_Salt_Electroneutral/")
   end
   
   # Note: All the expert software engineers told me to not define anything in the 
   #       global scope of the modules. They said I should define all variables in
   #       the functions. Something about mutable global state being evil.

   #include("Battery_Simulation_1D_generally_useful_functions.jl") 

   ############################################################
   #### Functions to create the initial system state, load a previous system state, and functions for V_anode and V_cathode
   include("Pulse_Diffusion_2D_create_system_state.jl")
   ############################################################


   ############################################################
   ###### A Function for the Pulse Voltage Operator 
   include("Pulse_Diffusion_2D_Simulation_Pulse_Voltage.jl")
   ############################################################


   ############################################################
   ###### A Function for the Pulse Current Operator 
   include("Pulse_Diffusion_2D_Simulation_Pulse_Current.jl")
   ############################################################


   ############################################################
   ###### A Function for the HYBRID Pulse Current Operator (this allows you to do a voltage pulse then an OCV pulse)
   include("Pulse_Diffusion_2D_Simulation_Pulse_Current_AND_Voltage.jl")
   ############################################################


   ############################################################
   ###### A Function for pulsing current from a capacitor
   include("Pulse_Diffusion_2D_Simulation_Pulse_from_Capacitor.jl")
   ############################################################


   ############################################################
   ###### A Function for Plotting Data
   include("Pulse_Diffusion_2D_Simulation_Plots.jl")
   ############################################################


end  



