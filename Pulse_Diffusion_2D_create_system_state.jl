"""
Create the data structure for the system state


"""


#### Equilibrium voltage of the electrode with respect to a reference electrode in the bulk (assuming zero electrostatic fields). The reference electrode is in the "bulk'
function V_eq(conc_A, conc_B, conc_A_re, conc_B_re )  
   return( 8.3*300/96500*log( conc_A/conc_A_re * conc_B_re/conc_B ) )  #The reaction is:      A + e- <-> B        mu_A - F*V_eq = mu_B          R*T*ln(conc_A) - F*V_eq = R*T*ln(conc_B) + C         V_eq = RT/F ln(conc_A/conc_B) +C        The reference electrode is in the "bulk'
end

#### Equilibrium concentration at the interface of the electrode.  It returns the conc_A that is in equilibrium with the electrode voltage
function conc_A_eq(electrode_voltage, total_conc, conc_A_re, conc_B_re )
   return(1/(1/(exp(electrode_voltage/8.3/300*96500)*conc_A_re/conc_B_re)+ 1)*total_conc)
end 


struct system_state_structure
   parent_dictionary  ::String
   accumulated_simulation_time ::Array{Float64,1}
   Diffusivity        ::Float64
   dx                 ::Float64
   dy                 ::Float64
   locations_x        ::Array{Float64,2}
   locations_y        ::Array{Float64,2}
   num_x_mps          ::Int
   num_y_mps          ::Int
   spike_num_x_mps    ::Int  #number of meshpoints
   spike_num_y_mps    ::Int  #number of meshpoints
   conc_A             ::Array{Float64,2}
   total_conc         ::Float64
   electrode_voltage  ::Array{Float64,1}  #this is the working electrode's voltage wrt the reference electrode 
   reaction_k         ::Float64
   Beta               ::Float64
end


function create_system_state()
   ############## Define system CONSTANTS that don't mutate or change at all during the ############################
   ############## calculations                                                          ############################
   println("Creating fresh system state.")

   dx = 200E-9
   dy = 200E-9
   num_x_mps= 64  #number of meshpoints
   num_y_mps= 500  #number of meshpoints

   locations_x = ones(num_y_mps)*transpose(collect(range(0.0, num_x_mps*dx, length=num_x_mps)))
   locations_y = reverse(collect(range(0.0, num_y_mps*dy, length=num_y_mps))*transpose(ones(num_x_mps)),dims=1)

   ## Describing the electrolyte domain
   Diffusivity=1.E-10 #in units of m^2/s

   ## Describing the spike 
   spike_num_x_mps = 30
   spike_num_y_mps = 150

   ############## Define system arrays that WILL mutate during calculations but won't change in size and don't  ##############
   ############## have a dimension size that is dependent on the simulation duration or cycling procedure !     ##############
   conc_A   =    ones(num_y_mps,num_x_mps)*1000; #(mol / m3) 
   conc_A[num_y_mps-spike_num_y_mps+2:end,1:spike_num_x_mps-1] .= -1234.0
   total_conc  = conc_A[1,50]*3.0
   reaction_k  = 0.001
   Beta    =  0.5

   #ss stands for system state
   system_state=system_state_structure(
      "virginnewstate_dictionary_results.jld2",   #parent_dictionary_name
      [0.0],                #accumulated_simulation_time
      Diffusivity,         #Diffusivity       
      dx,                  #dx
      dy,                  #dy
      locations_x,         #locations_x     
      locations_y,         #locations_y 
      num_x_mps,           #num_x_mps 
      num_y_mps,           #num_y_mps
      spike_num_x_mps,     #spike_num_x_mps
      spike_num_y_mps,     #spike_num_y_mps
      conc_A,              #concentration of A  
      total_conc,          #concentration of A  
      [0.0],               #this is the working electrode's voltage wrt the reference electrode 
      reaction_k,          #reaction_k
      Beta                 #Beta
   )

   return(system_state)
end




function load_previous_state(file_datenumber,iteration_number=0)
   ss       = get(FileIO.load("produced_data/"*file_datenumber*"_dictionary_results.jld2"), "system_state",0);
   sim_data = get(FileIO.load("produced_data/"*file_datenumber*"_dictionary_results.jld2"), "sim_data",0);
   if iteration_number != 0
      ss.conc_A[:,:] = 1.0*sim_data.conc_A_saved[iteration_number,:,:]
   end
   return(ss,sim_data)
end
