####################################################################################################
###### Pulse Voltage From Capacitor
####################################################################################################

###### Create the data type that will hold all information for the pulse voltage operator (including time series)
struct pulse_voltage_from_capacitor_simulation_data_structure
   start_time                      ::Array{String,1}
   stop_time                       ::Array{String,1}
   input_ss                        ::system_state_structure
   data_dictionary_name            ::String
   simulation_duration             ::Float64
   dt_biggest                      ::Float64
   saved_dt_spacing                ::Float64
   time_saved                      ::Array{Float64,1}
   electrode_voltage_saved         ::Array{Float64,1}
   overvoltage_saved               ::Array{Float64,2}
   current_density_saved           ::Array{Float64,2}
   current_density_time_average    ::Array{Float64,1}
   superficial_cd_saved            ::Array{Float64,1}
   superficial_cd_time_average     ::Array{Float64,1}
   main_loop_iteration_saved       ::Array{Float64,1}
   Charge_Passed_saved             ::Array{Float64,2}
   conc_A_saved                    ::Array{Float64,3}
   capacitor_voltage_saved         ::Array{Float64,1}
   capacitor_voltage               ::Array{Float64,1}
end


######## A function to run the simulation                            
function pulse_voltage_from_capacitor(ss, simulation_duration, dt_biggest, saved_dt_spacing, external_current_per_m2, capacitor_Farads_per_m2 , capacitor_initial_voltage, capacitor_voltage_ontime, ocv_ontime)   # units are always mks  e.g. current = A/m2,  Flux = moles/m2/s,  concentration = moles/m3
   start_time=Dates.format(Dates.now(),"yyyymmddHHMMSS")   #start_time

   if saved_dt_spacing <= dt_biggest; saved_dt_spacing = dt_biggest; end;
   save_data_time_thresholds = [ 0.0; 1.0*dt_biggest; 2.0*dt_biggest; 3.0*dt_biggest;  collect(range(0.0 , step=saved_dt_spacing , stop=simulation_duration - 0.5*saved_dt_spacing))[2:end] ]  #the "- 0.5*saved_dt_spacing" allows me to have just one timestep saved at the very end of the simulation (previously it would save the 2nd to final timestep and the final timestep)
   save_data_time_thresholds = cat(save_data_time_thresholds,1E9,dims=1);  #this final 1E9 second threshold will never be triggered.  It exists just so that the "time[1] + 1E-10 >= save_data_time_thresholds[simdata_i]" logic can be executed
   save_data_time_thresholds = save_data_time_thresholds .+ ss.accumulated_simulation_time[1]
   
   simdata=pulse_voltage_from_capacitor_simulation_data_structure(
      [start_time]                                                                            #start_time
      ,["running"]                                                                            #stop_time
      ,ss                                                                                     #input system state
      ,start_time * "_dictionary_results.jld2"                                                #data_dictionary_name
      ,simulation_duration                                                                    #simulation_duration
      ,dt_biggest                                                                             #dt_biggest
      ,saved_dt_spacing                                                                       #saved_dt_spacing
      ,zeros(length(save_data_time_thresholds))                                               #time_saved 
      ,zeros(length(save_data_time_thresholds))                                               #electrode_voltage_saved             
      ,zeros(length(save_data_time_thresholds),ss.num_x_mps + ss.spike_num_y_mps + 1 )        #overvoltage_saved
      ,zeros(length(save_data_time_thresholds),ss.num_x_mps + ss.spike_num_y_mps + 1 )        #current_density_saved    WE COUNT THE CORNER MESHPOINTS TWICE BECAUSE THEY HAVE INTERFACE POINTING IN THE DY DIRECTION AND THE DX DIRECTION
      ,zeros(ss.num_x_mps + ss.spike_num_y_mps + 1)                                           #current_density_time_average
      ,zeros(length(save_data_time_thresholds))                                               #superficial_cd_saved
      ,[0.0]                                                                                  #superficial_cd_time_average
      ,zeros(length(save_data_time_thresholds))                                               #main_loop_iteration_saved  
      ,zeros(length(save_data_time_thresholds),ss.num_x_mps + ss.spike_num_y_mps + 1)         #Charge_Passed_saved       WE COUNT THE CORNER MESHPOINT TWICE BECAUSE THEY HAVE INTERFACE POINTING IN THE DY DIRECTION AND THE DX DIRECTION
      ,zeros(length(save_data_time_thresholds),ss.num_y_mps,ss.num_x_mps)                     #conc_A_saved
      ,zeros(length(save_data_time_thresholds))                                               #capacitor_voltage_saved
      ,[capacitor_initial_voltage]                                                            #capacitor_voltage
   )

   println(" ");println(" ");println(" ");println(" ")
   println("Starting Pulse Voltage Mode.")
   println("data saved to produced_data/"*simdata.data_dictionary_name) 
   println("D dt/dx^2 is ",ss.Diffusivity*simdata.dt_biggest/ss.dx/ss.dx, " and must be less than 0.5")

   # Set some values that we need for the simulation loops 
   short_y_num = ss.num_y_mps - ss.spike_num_y_mps
   short_x_num = ss.num_x_mps - ss.spike_num_x_mps 
   conc_A_along_surface = zeros(ss.num_x_mps + ss.spike_num_y_mps + 1 )
   conc_A_along_surface[1:ss.spike_num_x_mps]                                       = ss.conc_A[short_y_num+1,1:ss.spike_num_x_mps]
   conc_A_along_surface[ss.spike_num_x_mps+1:ss.spike_num_x_mps+ss.spike_num_y_mps] = ss.conc_A[short_y_num+1:ss.num_y_mps,ss.spike_num_x_mps]
   conc_A_along_surface[ss.spike_num_x_mps+ss.spike_num_y_mps+1:end]                = ss.conc_A[ss.num_y_mps,ss.spike_num_x_mps:ss.num_x_mps]
   conc_B = 1.0*ss.conc_A
   conc_B_along_surface                                                             = 1.0*conc_A_along_surface[:]
   conc_A_eq_along_surface                                                          = [ conc_A_eq(ss.electrode_voltage[1], conc_B[end,end], ss.conc_A[1,50], conc_B[1,50] ) ]  #This calculation is pointless in the Binary Salt situation.
   r_x                                                   = ss.Diffusivity * dt_biggest / ss.dx / ss.dx
   r_y                                                   = ss.Diffusivity * dt_biggest / ss.dy / ss.dy
   voltage_eq_along_surface                              = V_eq.(conc_A_along_surface, conc_B_along_surface, ss.conc_A[1,50], conc_B[1,50])  #The reference electrode is located at [1,50]
   overvoltage                                           = ss.electrode_voltage[1] .- voltage_eq_along_surface[:]
   current_density                                       = Current_Density.(ss.reaction_k, ss.Beta, conc_A_along_surface[:], conc_B_along_surface[:], overvoltage[:]) #(A/m2)
   Charge_Passed                                         = 0.0*conc_A_along_surface
   superficial_current_density                           = [0.0]
   molar_flux                                            = 0.0*conc_A_along_surface
   time                                                  = [ss.accumulated_simulation_time[1]]
   simdata_i                                             = 1
   dt_reduction_factor                                   = 1
   conc_next_timestep                                    = copy(ss.conc_A)
   conc_A_reduced_dt                                     = copy(ss.conc_A)
   conc_A_along_surface_trial                            = copy(conc_A_along_surface)
   conc_increment_dx                                     = 0.0*ss.conc_A
   conc_increment_dy                                     = 0.0*ss.conc_A
   superficial_current_density_target                    = [0.0]
   superficial_current_density_target_previous           = [-1E6]
   superficial_current_density_target_previous_previous  = [-1E6]
   electrode_voltage_previous                            = [-1E6] 
   electrode_voltage_previous_previous                   = [-1E6]
   simdata.capacitor_voltage_saved[1]                    = capacitor_initial_voltage
   electrode_voltage_limit                               = 0.3
   
   ########### This for loop increments time.   It's an EXPLICIT simulation. It uses a forward Euler time marching scheme. 
   main_loop_iteration = 0
   while time[1] <=  ss.accumulated_simulation_time[1] + simulation_duration

      # Save data from the previous timestep 
      if time[1] + 1E-13 >= save_data_time_thresholds[simdata_i]  #the + 1E-13 is because the computer can't store perfect numbers, e.g. 3E-5 can only be stored as 3.0000000000000004e-5
         @printf(":%-9i   real_time:%+0.3e   dt_red_fac:%-3i    conc_eq:%+0.7e   conc_corner:%+0.7e    elctrd_volt:%+0.7e   V_eq_corner:%+0.7e     sup_cd:%+0.4e\n", main_loop_iteration , time[1], dt_reduction_factor, conc_A_eq_along_surface[1], ss.conc_A[end,30],  ss.electrode_voltage[1], voltage_eq_along_surface[180], superficial_current_density[1])
         record_output(ss, simdata, simdata_i, time[1], main_loop_iteration, ss.electrode_voltage[1],  current_density, superficial_current_density[1], Charge_Passed, overvoltage, conc_A_along_surface, simdata.capacitor_voltage[1])
         simdata_i = simdata_i + 1
      end

      # Increment the loop and time identifiers for the upcoming timestep
      main_loop_iteration = main_loop_iteration + 1
      time[1] = time[1] + simdata.dt_biggest

      # Update the voltage or current density target
      superficial_current_density_target_previous_previous[1] = superficial_current_density_target_previous[1]
      superficial_current_density_target_previous[1]          = superficial_current_density_target[1]
      electrode_voltage_previous_previous[1]                  = electrode_voltage_previous[1]
      electrode_voltage_previous[1]                           = ss.electrode_voltage[1]
      if mod(time[1], capacitor_voltage_ontime + ocv_ontime) <= capacitor_voltage_ontime
         ss.electrode_voltage[1] = simdata.capacitor_voltage[1]
         simdata.capacitor_voltage[1] = simdata.capacitor_voltage[1] + simdata.dt_biggest * ( external_current_per_m2 - superficial_current_density[1] ) / capacitor_Farads_per_m2
      else
         simdata.capacitor_voltage[1] = simdata.capacitor_voltage[1] + simdata.dt_biggest * ( external_current_per_m2 -         0.0                    ) / capacitor_Farads_per_m2
         # Put the battery in open circuit
         # Every few microseconds the "Arbin" circuitry adjusts ss.electrode_voltage[1] so that superficial_current_density equals superficial_current_density_target[1], then afterwards the current will stray while ss.electrode_voltage[1] is held constant until the next time the Arbin enforces the correct current
         # Predict (guess) the ss.electrode_voltage[1] value that might create the target current density (superficial_current_density_target[1]) based on linear regression of the previous two values of ss.electrode_voltage[1]
         if superficial_current_density_target[1] == superficial_current_density_target_previous_previous[1]
            ss.electrode_voltage[1] = electrode_voltage_previous[1] + (electrode_voltage_previous[1] - electrode_voltage_previous_previous[1])
            if ss.electrode_voltage[1] > +abs(electrode_voltage_limit); ss.electrode_voltage[1] = +abs(electrode_voltage_limit); end;
            if ss.electrode_voltage[1] < -abs(electrode_voltage_limit); ss.electrode_voltage[1] = -abs(electrode_voltage_limit); end;
         end 
         # Correct (aka adjust) the value of ss.electrode_voltage[1] until it creates exactly the correct current, then afterwards hold ss.electrode_voltage[1] steady whil the current strays away from the target value until the Arbin once again enforced the correct current
         voltage_eq_along_surface[:] = V_eq.(conc_A_along_surface, conc_B_along_surface, ss.conc_A[1,50], conc_B[1,50])  #The reference electrode is located at [1,50]
         overvoltage[:] = ss.electrode_voltage[1] .- voltage_eq_along_surface[:]
         current_density[:] = Current_Density.(ss.reaction_k, ss.Beta, conc_A_along_surface[:], conc_B_along_surface[:], overvoltage[:]) #(A/m2)
         superficial_current_density[1] = mean(current_density[:])*length(current_density[:])/ss.num_x_mps
         superficial_current_density_error = superficial_current_density - superficial_current_density_target
         while (abs(superficial_current_density_error[1]) > 0.1) & (abs(ss.electrode_voltage[1]) < abs(electrode_voltage_limit)) # the 0.3 is a voltage limit, to prevent voltage from blowing up.  But since this simulation uses OCV, voltage should never blow up.
            ss.electrode_voltage[1] = ss.electrode_voltage[1] - superficial_current_density_error[1]*1E-7
            overvoltage[:] = ss.electrode_voltage[1] .- voltage_eq_along_surface
            current_density[:] = Current_Density.(ss.reaction_k, ss.Beta, conc_A_along_surface[:], conc_B_along_surface[:], overvoltage[:])  #(A/m2)
            superficial_current_density[1]       = mean(current_density[:])*length(current_density[:])/ss.num_x_mps
            superficial_current_density_error = superficial_current_density - superficial_current_density_target
            #@printf("loop:%5.0i   superficial_current_density:%+0.7e   electrode_voltage:%+0.7e   target_cd:%+0.7e   superficial_cd_error:%+0.7e\n", main_loop_iteration, superficial_current_density , ss.electrode_voltage[1] , superficial_current_density_target , superficial_current_density_error)
         end
      end
      
      conc_A_eq_along_surface[1] = conc_A_eq(ss.electrode_voltage[1], conc_B[end,end], ss.conc_A[1,50], conc_B[1,50] ) 
      molar_flux[:] = current_density[:]/96500.0

      ####### Calculate the change in concentration at the interfacial locations.  "Reduced" means the reduction in dt so that overshooting instability doesn't kill the simulation.
      conc_A_along_surface_trial[:]    = conc_A_along_surface[:]
      conc_next_timestep[:,:] = ss.conc_A[:,:]
      ## Trial calculation of the next timestep's interfacial concentrations
      for i_x in 1:ss.spike_num_x_mps - 1
         conc_next_timestep[short_y_num+1, i_x]                = 2*r_y*ss.conc_A[short_y_num      , i_x]          + (1 - 2*r_y)*ss.conc_A[short_y_num+1, i_x]                                                            +   2*r_y*ss.dy*molar_flux[i_x]/ss.Diffusivity
      end
      i_x = ss.spike_num_x_mps
      conc_next_timestep[short_y_num+1, i_x]                   =4/3*r_y*ss.conc_A[short_y_num, i_x]               + (1 - 2*r_y)*ss.conc_A[short_y_num+1, i_x]          +2/3*r_y*ss.conc_A[short_y_num+2, i_x]         + 2/3*r_y*ss.dy*molar_flux[ss.spike_num_x_mps]/ss.Diffusivity  #This line is for the corner mesh point.  The 0.5 factor is because this meshpoint only has electrode interface on half it's area.
      conc_next_timestep[ss.num_y_mps,  i_x]                   =  2*r_y*ss.conc_A[ss.num_y_mps - 1, i_x]          + (1 - 2*r_y)*ss.conc_A[ss.num_y_mps, i_x]                                                             +   2*r_y*ss.dy*molar_flux[ss.spike_num_x_mps + ss.spike_num_y_mps+1]/ss.Diffusivity  #This line is for the corner mesh point.  The 0.5 factor is because this meshpoint only has electrode interface on half it's area.
      for i_x in ss.spike_num_x_mps+1:ss.num_x_mps                                           
         conc_next_timestep[ss.num_y_mps      , i_x]           =  2*r_y*ss.conc_A[ss.num_y_mps - 1  , i_x]        + (1 - 2*r_y)*ss.conc_A[ss.num_y_mps, i_x]                                                             +   2*r_y*ss.dy*molar_flux[i_x + ss.spike_num_y_mps + 1]/ss.Diffusivity
      end
      conc_increment_dy = conc_next_timestep[:,:] - ss.conc_A[:,:]       #This line is for all mesh points 
      conc_next_timestep[:,:] = ss.conc_A[:,:]
      i_y = short_y_num+1
      conc_next_timestep[i_y, ss.spike_num_x_mps]              =2/3*r_x*ss.conc_A[i_y , ss.spike_num_x_mps-1]     + (1 - 2*r_x)*ss.conc_A[i_y , ss.spike_num_x_mps]    +4/3*r_x*ss.conc_A[i_y , ss.spike_num_x_mps+1] - 2/3*r_x*ss.dx*(-molar_flux[ss.spike_num_x_mps+1]/ss.Diffusivity)
      for i_y in short_y_num+2:ss.num_y_mps-1
         conc_next_timestep[i_y , ss.spike_num_x_mps]          =  2*r_x*ss.conc_A[i_y , ss.spike_num_x_mps + 1]   + (1 - 2*r_x)*ss.conc_A[i_y , ss.spike_num_x_mps]                                                      -   2*r_x*ss.dx*(-molar_flux[ss.spike_num_x_mps + i_y - short_y_num]/ss.Diffusivity)
      end
      i_y = ss.num_y_mps
      conc_next_timestep[i_y , ss.spike_num_x_mps]             =  2*r_x*ss.conc_A[i_y , ss.spike_num_x_mps+1]     + (1 - 2*r_x)*ss.conc_A[i_y , ss.spike_num_x_mps]                                                      -   2*r_x*ss.dx*(-molar_flux[ss.spike_num_x_mps + ss.spike_num_y_mps + 1]/ss.Diffusivity)
      conc_increment_dx = conc_next_timestep[:,:] - ss.conc_A[:,:]
      conc_A_along_surface_trial[1:ss.spike_num_x_mps]                                       = conc_A_along_surface_trial[1:ss.spike_num_x_mps]                                       + conc_increment_dy[short_y_num+1,1:ss.spike_num_x_mps]            + conc_increment_dx[short_y_num+1,1:ss.spike_num_x_mps] 
      conc_A_along_surface_trial[ss.spike_num_x_mps+1:ss.spike_num_x_mps+ss.spike_num_y_mps] = conc_A_along_surface_trial[ss.spike_num_x_mps+1:ss.spike_num_x_mps+ss.spike_num_y_mps] + conc_increment_dy[short_y_num+1:ss.num_y_mps,ss.spike_num_x_mps] + conc_increment_dx[short_y_num+1:ss.num_y_mps,ss.spike_num_x_mps] 
      conc_A_along_surface_trial[ss.spike_num_x_mps+ss.spike_num_y_mps+1:end]                = conc_A_along_surface_trial[ss.spike_num_x_mps+ss.spike_num_y_mps+1:end]                + conc_increment_dy[ss.num_y_mps,ss.spike_num_x_mps:ss.num_x_mps]  + conc_increment_dx[ss.num_y_mps,ss.spike_num_x_mps:ss.num_x_mps]
      ## Next, calculate the reduction of dt_biggest that is necessary to avoid the interfacial concentration overshooting the equilibrium value
      dt_reduction_factor = Int( maximum( cat( 1,round.( abs.( 3*(conc_A_along_surface .- conc_A_along_surface_trial) ./ conc_A_along_surface ) ),dims=1)))  # The 3 is to say we don't want the concentration changing by more than 30% of the distance to its equilibrium.   The ( (conc_A_eq_along_surface .- conc_A_along_surface) / conc_A_along_surface ) is to avoid the situation where conc_A_eq_along_surface hasn't been changing in time thus (conc_A_eq_along_surface .- conc_A_along_surface) is super small and tiny noise in current_density causes overshoot.
      ## Next, calculate the time-evolution of interfacial concentration with "reduced" timesteps
      dt_reduced = simdata.dt_biggest / dt_reduction_factor
      r_x_reduced = ss.Diffusivity * dt_reduced / ss.dx / ss.dx
      r_y_reduced = ss.Diffusivity * dt_reduced / ss.dy / ss.dy   
      conc_next_timestep[:,:] = ss.conc_A[:,:]
      conc_A_reduced_dt[:,:]  = ss.conc_A[:,:]
      overvoltage_reduced_dt_average                 = 0.0*overvoltage[:]                
      current_density_reduced_dt_average             = 0.0*current_density[:]            
      conc_A_along_surface_reduced_dt_average        = 0.0*conc_A_along_surface[:]
      for t = 1:dt_reduction_factor
         #Fickian change in interfacial concentrations due to y-direction gradients
         conc_next_timestep[:,:] = conc_A_reduced_dt[:,:]
         for i_x in 1:ss.spike_num_x_mps - 1
            conc_next_timestep[short_y_num+1, i_x]                = 2*r_y_reduced*conc_A_reduced_dt[short_y_num      , i_x]          + (1 - 2*r_y_reduced)*conc_A_reduced_dt[short_y_num+1, i_x]                                                                       +   2*r_y_reduced*ss.dy*molar_flux[i_x]/ss.Diffusivity
         end
         i_x = ss.spike_num_x_mps
         conc_next_timestep[short_y_num+1, i_x]                   =4/3*r_y_reduced*conc_A_reduced_dt[short_y_num, i_x]               + (1 - 2*r_y_reduced)*conc_A_reduced_dt[short_y_num+1, i_x]        +2/3*r_y_reduced*conc_A_reduced_dt[short_y_num+2, i_x]         + 2/3*r_y_reduced*ss.dy*molar_flux[ss.spike_num_x_mps]/ss.Diffusivity  #This line is for the corner mesh point.  The 0.5 factor is because this meshpoint only has electrode interface on half it's area.
         conc_next_timestep[ss.num_y_mps,  i_x]                   =  2*r_y_reduced*conc_A_reduced_dt[ss.num_y_mps - 1, i_x]          + (1 - 2*r_y_reduced)*conc_A_reduced_dt[ss.num_y_mps, i_x]                                                                        +   2*r_y_reduced*ss.dy*molar_flux[ss.spike_num_x_mps + ss.spike_num_y_mps+1]/ss.Diffusivity  #This line is for the corner mesh point.  The 0.5 factor is because this meshpoint only has electrode interface on half it's area.
         for i_x in ss.spike_num_x_mps+1:ss.num_x_mps                                           
            conc_next_timestep[ss.num_y_mps      , i_x]           =  2*r_y_reduced*conc_A_reduced_dt[ss.num_y_mps - 1  , i_x]        + (1 - 2*r_y_reduced)*conc_A_reduced_dt[ss.num_y_mps, i_x]                                                                        +   2*r_y_reduced*ss.dy*molar_flux[i_x + ss.spike_num_y_mps + 1]/ss.Diffusivity
         end
         conc_increment_dy = conc_next_timestep[:,:] - conc_A_reduced_dt[:,:]       #This line is for all mesh points 
         #Fickian change in interfaical concentrations due to x-direction gradients
         conc_next_timestep[:,:] = conc_A_reduced_dt[:,:]
         i_y = short_y_num+1
         conc_next_timestep[i_y, ss.spike_num_x_mps]              =2/3*r_x_reduced*conc_A_reduced_dt[i_y , ss.spike_num_x_mps-1]     + (1 - 2*r_x_reduced)*conc_A_reduced_dt[i_y , ss.spike_num_x_mps]  +4/3*r_x_reduced*conc_A_reduced_dt[i_y , ss.spike_num_x_mps+1] - 2/3*r_x_reduced*ss.dx*(-molar_flux[ss.spike_num_x_mps+1]/ss.Diffusivity)
         for i_y in short_y_num+2:ss.num_y_mps-1
            conc_next_timestep[i_y , ss.spike_num_x_mps]          =  2*r_x_reduced*conc_A_reduced_dt[i_y , ss.spike_num_x_mps + 1]   + (1 - 2*r_x_reduced)*conc_A_reduced_dt[i_y , ss.spike_num_x_mps]                                                                 -   2*r_x_reduced*ss.dx*(-molar_flux[ss.spike_num_x_mps + i_y - short_y_num]/ss.Diffusivity)
         end
         i_y = ss.num_y_mps
         conc_next_timestep[i_y , ss.spike_num_x_mps]             =  2*r_x_reduced*conc_A_reduced_dt[i_y , ss.spike_num_x_mps+1]     + (1 - 2*r_x_reduced)*conc_A_reduced_dt[i_y , ss.spike_num_x_mps]                                                                 -   2*r_x_reduced*ss.dx*(-molar_flux[ss.spike_num_x_mps + ss.spike_num_y_mps + 1]/ss.Diffusivity)
         conc_increment_dx = conc_next_timestep[:,:] - conc_A_reduced_dt[:,:]
         # Add the x-direction and y-direction effects together!
         conc_A_reduced_dt[:,:] = conc_A_reduced_dt[:,:] + conc_increment_dy[:,:] + conc_increment_dx[:,:]

         #Calculate for the next timestep (dt_reduced) the values of conc_A_along_surface, voltage_eq_along_surface, overvoltage, current_density, conc_next_timestep,  for each "reduced" time step for the interfacial mesh points
         conc_A_along_surface[1:ss.spike_num_x_mps]                                       = conc_A_reduced_dt[short_y_num+1,1:ss.spike_num_x_mps]
         conc_A_along_surface[ss.spike_num_x_mps+1:ss.spike_num_x_mps+ss.spike_num_y_mps] = conc_A_reduced_dt[short_y_num+1:ss.num_y_mps,ss.spike_num_x_mps]
         conc_A_along_surface[ss.spike_num_x_mps+ss.spike_num_y_mps+1:end]                = conc_A_reduced_dt[ss.num_y_mps,ss.spike_num_x_mps:ss.num_x_mps]
         conc_B_along_surface[:] = conc_A_along_surface[:]
         voltage_eq_along_surface[:] = V_eq.(conc_A_along_surface, conc_B_along_surface, conc_A_reduced_dt[1,50], conc_A_reduced_dt[1,50])  #The reference electrode is located at [1,50]
         overvoltage[:] = ss.electrode_voltage[1] .- voltage_eq_along_surface
         current_density[:] = Current_Density.(ss.reaction_k, ss.Beta, conc_A_along_surface[:], conc_B_along_surface[:], overvoltage[:])  #(A/m2)
         molar_flux[:] = current_density[:]/96500.0

         #Build the time-averaged values during these "reduced" time steps
         conc_A_along_surface_reduced_dt_average[:]     = conc_A_along_surface_reduced_dt_average[:]     + conc_A_along_surface[:]        / dt_reduction_factor       
         overvoltage_reduced_dt_average[:]              = overvoltage_reduced_dt_average[:]              + overvoltage[:]                 / dt_reduction_factor                
         current_density_reduced_dt_average[:]          = current_density_reduced_dt_average[:]          + current_density[:]             / dt_reduction_factor           

      end   #end of for t = 1:dt_reduction_factor


      #Revert to the time-averaged values of conc_A_along_surface, overvoltage, current_density, molar_flux, superficial_current_density along the interface during this timestep
      conc_A_along_surface[:]                = conc_A_along_surface_reduced_dt_average[:]
      overvoltage[:]                         = overvoltage_reduced_dt_average[:]
      current_density[:]                     = current_density_reduced_dt_average[:]
      simdata.current_density_time_average[:]= simdata.current_density_time_average[:] + current_density[:] / (simulation_duration/dt_biggest)
      superficial_current_density[1]         = mean(current_density[:])*length(current_density[:])/ss.num_x_mps
      simdata.superficial_cd_time_average[1] = simdata.superficial_cd_time_average[1] + superficial_current_density[1]/(simulation_duration/dt_biggest)
      molar_flux[:]                          = current_density[:]/96500.0
      Charge_Passed[:]                       = Charge_Passed[:] + current_density*simdata.dt_biggest

      ##### Fickian change in concentration due to y-direction gradients.     FOR DISCRETIZATION HELP SEE 2014 DISSERTATION BY ZANGANA !! IF there’s a flux boundary condition at x = 0 then: du/dt = (2*D*u_i+1,j - 2*D*u_i,j )/deltax^2 - 2*D*g/deltax,  u_i,j+1 = u_i,j + 2r*u_i+1,j - 2r*u_i,j - 2r*g*deltax   where g is du/dx on the boundary at x =0 and the 2 comes from the fact that we divide by 0.5deltax, not full deltax (because of boundary at x=0).  If there's a LHS "half boundary" due to a corner on the spike then:  du/dt = ( D*( u_i+1,j - u_i,j )/deltax - (0.5*D*g + 0.5*D*( u_i,j - u_i-1,j )/deltax) )/0.75/deltax,  du/dt = 4/3*D*( u_i+1,j - u_i,j )/deltax^2 - 2/3*D*( u_i,j - u_i-1,j )/deltax^2 - 2/3*D*g/deltax ,   u_i,j+1 = u_i,j + 4/3*r*( u_i+1,j - u_i,j ) - 2/3*r*( u_i,j - u_i-1,j ) - 2/3*r*g*deltax = 2/3*r*u_i-1,j + (1 - 2*r)*u_i,j + 4/3*r*u_i+1,j - 2/3*r*g*deltax       and if the boundary condition is a RHS "half boundary" then this is: u_i,j+1 = 4/3*r*u_i-1,j + (1 - 2*r)*u_i,j + 2/3*r*u_i+1,j + 2/3*r*g*deltax
      conc_next_timestep[:,:] = ss.conc_A[:,:]
      for i_x in 1:ss.spike_num_x_mps - 1
         conc_next_timestep[2:short_y_num, i_x]                        =   r_y*ss.conc_A[1:short_y_num - 1, i_x]                      + (1 - 2*r_y)*ss.conc_A[2:short_y_num, i_x]                          +    r_y*ss.conc_A[3:short_y_num+1 , i_x]
         #conc_next_timestep[short_y_num+1, i_x]                        = 2*r_y*ss.conc_A[short_y_num      , i_x]                     + (1 - 2*r_y)*ss.conc_A[short_y_num+1, i_x]                                                                                         +   2*r_y*ss.dy*molar_flux[i_x]/ss.Diffusivity
      end
      i_x = ss.spike_num_x_mps
      conc_next_timestep[2:short_y_num, i_x]                           =   r_y*ss.conc_A[1:short_y_num - 1, i_x]                      + (1 - 2*r_y)*ss.conc_A[2:short_y_num, i_x]                          +    r_y*ss.conc_A[3:short_y_num+1, i_x]
      #conc_next_timestep[short_y_num+1, i_x]                           =4/3*r_y*ss.conc_A[short_y_num, i_x]                          + (1 - 2*r_y)*ss.conc_A[short_y_num+1, i_x]                          +2/3*r_y*ss.conc_A[short_y_num+2, i_x]                         + 2/3*r_y*ss.dy*molar_flux[ss.spike_num_x_mps]/ss.Diffusivity  #This line is for the corner mesh point.  The 0.5 factor is because this meshpoint only has electrode interface on half it's area.
      conc_next_timestep[short_y_num+2:ss.num_y_mps-1,  i_x]           =   r_y*ss.conc_A[short_y_num+1:ss.num_y_mps-2, i_x]           + (1 - 2*r_y)*ss.conc_A[short_y_num+2:ss.num_y_mps-1, i_x]           +    r_y*ss.conc_A[short_y_num+3:ss.num_y_mps , i_x]
      #conc_next_timestep[ss.num_y_mps,  i_x]                           = 2*r_y*ss.conc_A[ss.num_y_mps - 1, i_x]                      + (1 - 2*r_y)*ss.conc_A[ss.num_y_mps, i_x]                                                                                          +   2*r_y*ss.dy*molar_flux[ss.spike_num_x_mps + ss.spike_num_y_mps+1]/ss.Diffusivity  #This line is for the corner mesh point.  The 0.5 factor is because this meshpoint only has electrode interface on half it's area.
      for i_x in ss.spike_num_x_mps+1:ss.num_x_mps                                           
         conc_next_timestep[2:ss.num_y_mps - 1, i_x]                   =   r_y*ss.conc_A[1:ss.num_y_mps - 2, i_x]                     + (1 - 2*r_y)*ss.conc_A[2:ss.num_y_mps - 1, i_x]                     +    r_y*ss.conc_A[3:ss.num_y_mps , i_x]
         #conc_next_timestep[ss.num_y_mps      , i_x]                   =2*r_y*ss.conc_A[ss.num_y_mps - 1  , i_x]                     + (1 - 2*r_y)*ss.conc_A[ss.num_y_mps, i_x]                                                                                          +   2*r_y*ss.dy*molar_flux[i_x + ss.spike_num_y_mps + 1]/ss.Diffusivity
      end
      conc_increment_dy = conc_next_timestep[:,:] - ss.conc_A[:,:]       #This line is for all mesh points 
      ##### Fickian change in concentration due to x-direction gradients
      conc_next_timestep[:,:] = ss.conc_A[:,:]
      for i_y in 1:short_y_num
         conc_next_timestep[i_y, 1]                                    =   r_x*ss.conc_A[i_y , 2]                                     + (1 - 2*r_x)*ss.conc_A[i_y , 1]                                     +    r_x*ss.conc_A[i_y , 2]                                    #this is due to the "periodic" boundary condition
         conc_next_timestep[i_y, end]                                  =   r_x*ss.conc_A[i_y , end-1]                                 + (1 - 2*r_x)*ss.conc_A[i_y , end]                                   +    r_x*ss.conc_A[i_y , end-1]                                #this is due to the "periodic" boundary condition
         conc_next_timestep[i_y, 2:ss.num_x_mps - 1]                   =   r_x*ss.conc_A[i_y , 1:ss.num_x_mps - 2]                    + (1 - 2*r_x)*ss.conc_A[i_y , 2:ss.num_x_mps - 1]                    +    r_x*ss.conc_A[i_y , 3:ss.num_x_mps ]
      end
      i_y = short_y_num+1
      conc_next_timestep[i_y, 1]                                       =   r_x*ss.conc_A[i_y , 2]                                     + (1 - 2*r_x)*ss.conc_A[i_y , 1]                                     +    r_x*ss.conc_A[i_y , 2]                                    #this is due to the "periodic" boundary condition
      conc_next_timestep[i_y, 2:ss.spike_num_x_mps-1]                  =   r_x*ss.conc_A[i_y , 1:ss.spike_num_x_mps-2]                + (1 - 2*r_x)*ss.conc_A[i_y , 2:ss.spike_num_x_mps - 1]              +    r_x*ss.conc_A[i_y , 3:ss.spike_num_x_mps ]  
      #conc_next_timestep[i_y, ss.spike_num_x_mps]                      =2/3*r_x*ss.conc_A[i_y , ss.spike_num_x_mps-1]                + (1 - 2*r_x)*ss.conc_A[i_y , ss.spike_num_x_mps]                    +4/3*r_x*ss.conc_A[i_y , ss.spike_num_x_mps+1]                 - 2/3*r_x*ss.dx*(-molar_flux[ss.spike_num_x_mps+1]/ss.Diffusivity)
      conc_next_timestep[i_y, ss.spike_num_x_mps+1:ss.num_x_mps-1]     =   r_x*ss.conc_A[i_y , ss.spike_num_x_mps:ss.num_x_mps-2]     + (1 - 2*r_x)*ss.conc_A[i_y , ss.spike_num_x_mps+1:ss.num_x_mps-1]   +    r_x*ss.conc_A[i_y , ss.spike_num_x_mps+2:ss.num_x_mps ]  
      conc_next_timestep[i_y, end]                                     =   r_x*ss.conc_A[i_y , end-1]                                 + (1 - 2*r_x)*ss.conc_A[i_y , end]                                   +    r_x*ss.conc_A[i_y , end-1]                                #this is due to the "periodic" boundary condition
      for i_y in short_y_num+2:ss.num_y_mps-1
         #conc_next_timestep[i_y , ss.spike_num_x_mps]                  = 2*r_x*ss.conc_A[i_y , ss.spike_num_x_mps + 1]               + (1 - 2*r_x)*ss.conc_A[i_y , ss.spike_num_x_mps]                                                                                   -  2*r_x*ss.dx*(-molar_flux[ss.spike_num_x_mps + i_y - short_y_num]/ss.Diffusivity)
         conc_next_timestep[i_y , ss.spike_num_x_mps+1:ss.num_x_mps-1] =   r_x*ss.conc_A[i_y , ss.spike_num_x_mps:ss.num_x_mps - 2]   + (1 - 2*r_x)*ss.conc_A[i_y , ss.spike_num_x_mps+1:ss.num_x_mps-1]   +    r_x*ss.conc_A[i_y , ss.spike_num_x_mps+2:ss.num_x_mps] 
         conc_next_timestep[i_y , ss.num_x_mps]                        =   r_x*ss.conc_A[i_y , ss.num_x_mps - 1]                      + (1 - 2*r_x)*ss.conc_A[i_y , ss.num_x_mps]                          +    r_x*ss.conc_A[i_y , ss.num_x_mps - 1]                     #this is due to the "periodic" boundary condition
      end
      i_y = ss.num_y_mps
      #conc_next_timestep[i_y , ss.spike_num_x_mps]                     = 2*r_x*ss.conc_A[i_y , ss.spike_num_x_mps+1]                 + (1 - 2*r_x)*ss.conc_A[i_y , ss.spike_num_x_mps]                                                                                   -  2*r_x*ss.dx*(-molar_flux[ss.spike_num_x_mps + ss.spike_num_y_mps + 1]/ss.Diffusivity)
      conc_next_timestep[i_y , ss.spike_num_x_mps+1:ss.num_x_mps-1]    =   r_x*ss.conc_A[i_y , ss.spike_num_x_mps:ss.num_x_mps - 2]   + (1 - 2*r_x)*ss.conc_A[i_y , ss.spike_num_x_mps+1:ss.num_x_mps-1]   +    r_x*ss.conc_A[i_y , ss.spike_num_x_mps+2:ss.num_x_mps]   
      conc_next_timestep[i_y , ss.num_x_mps]                           =   r_x*ss.conc_A[i_y , ss.num_x_mps - 1]                      + (1 - 2*r_x)*ss.conc_A[i_y , ss.num_x_mps]                          +    r_x*ss.conc_A[i_y , ss.num_x_mps - 1]                     #this is due to the "periodic" boundary condition
      conc_increment_dx = conc_next_timestep[:,:] - ss.conc_A[:,:]

      ##### Add the Fickian change in concentration due to x-direction and y-direction gradients together!
      ss.conc_A[:,:] = ss.conc_A[:,:] + conc_increment_dy[:,:] + conc_increment_dx[:,:]
      conc_B[:,:] = ss.conc_A[:,:]

      ##### Pull the values for interfacial concentration from the calculations that used finer time steps
      ss.conc_A[short_y_num+1,1:ss.spike_num_x_mps]              = conc_A_along_surface[1:ss.spike_num_x_mps]
      ss.conc_A[short_y_num+1:ss.num_y_mps,ss.spike_num_x_mps]   = conc_A_along_surface[ss.spike_num_x_mps+1:ss.spike_num_x_mps+ss.spike_num_y_mps]
      ss.conc_A[ss.num_y_mps,ss.spike_num_x_mps:ss.num_x_mps]    = conc_A_along_surface[ss.spike_num_x_mps+ss.spike_num_y_mps+1:end]
      conc_B_along_surface[:]                                    = conc_A_along_surface[:]

   end   #end of time stepping: while time[1] <=  ss.accumulated_simulation_time[1] + simulation_duration

   #Print and store the final data      
   @printf(":%-9i   real_time:%+0.3e   conc_eq:%+0.7e   conc_corner:%+0.7e    elctrd_volt:%+0.7e\n", main_loop_iteration , time[1], conc_A_eq_along_surface[1], ss.conc_A[end,30],  ss.electrode_voltage[1] )
   ss.accumulated_simulation_time[1] = time[1]
   ss.parent_operation_dictionary[1] = simdata.data_dictionary_name
   record_output(ss, simdata, simdata_i, time[1], main_loop_iteration, ss.electrode_voltage[1],  current_density, superficial_current_density[1], Charge_Passed, overvoltage, conc_A_along_surface, simdata.capacitor_voltage[1])
   ### Save the results to hard disk
   simdata.stop_time[1]=Dates.format(Dates.now(),"yyyymmddHHMMSS")
   save("produced_data/"*simdata.data_dictionary_name, Dict("system_state"=>ss,"simdata"=>simdata))

   ### Print some basic output
   print("simulated duration is ",simulation_duration, " seconds\n")
   println("using FileIO")
   println("simdata = get(load(\"produced_data/"*simdata.data_dictionary_name*"\"), \"simdata\",0);")
   println("ss =       get(load(\"produced_data/"*simdata.data_dictionary_name*"\"), \"system_state\",0);")
   println("ss,simdata = Diffusion_2D.load_previous_state(\""*simdata.data_dictionary_name[1:14]*"\");")

   ### Give the caller the resulting data
   return(ss,simdata)

end  ## end of function pulse_voltage(...) 




####### A function to record a time series of the battery state during the CV operation
function record_output(ss, simdata, simdata_i , time, main_loop_iteration, electrode_voltage, current_density, superficial_current_density, Charge_Passed, overvoltage, conc_A_along_surface, capacitor_voltage)
      simdata.time_saved[simdata_i]               = time
      simdata.main_loop_iteration_saved[simdata_i]= main_loop_iteration
      simdata.electrode_voltage_saved[simdata_i]  = electrode_voltage
      simdata.current_density_saved[simdata_i,:]  = current_density[:]
      simdata.superficial_cd_saved[simdata_i]     = superficial_current_density  #superficial area to real area
      simdata.Charge_Passed_saved[simdata_i,:]    = Charge_Passed[:]   # coulombs per m3
      simdata.overvoltage_saved[simdata_i,:]      = overvoltage[:]
      simdata.conc_A_saved[simdata_i,:,:]         = ss.conc_A[:,:]
      simdata.capacitor_voltage_saved[simdata_i]          = capacitor_voltage
end
