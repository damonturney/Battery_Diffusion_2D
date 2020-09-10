####################################################################################################
###### Pulse Voltage Operator
####################################################################################################


###### Create the data type that will hold all information for the pulse voltage operator (including time series)
struct pulse_voltage_simulation_data_structure
   start_time                      ::Array{String,1}
   stop_time                       ::Array{String,1}
   data_dictionary_name            ::String
   dt                              ::Float64
   iterations                      ::UnitRange{Int64}
   saved_iteration_spacing         ::Float64
   electrode_voltage_saved         ::Array{Float64,1}
   overvoltage_saved               ::Array{Float64,2}
   current_density_saved           ::Array{Float64,2}
   superficial_cd_saved            ::Array{Float64,1}
   iterations_saved                ::Array{Float64,1}
   time_real_saved                 ::Array{Float64,1}
   main_loop_iteration_saved       ::Array{Float64,1}
   Charge_Passed_saved             ::Array{Float64,2}
   conc_A_saved                      ::Array{Float64,3}
   conc_A_along_surface            ::Array{Float64,2}
end


######## A function to create and initialize a single instance of the data structure 
function pulse_voltage(ss, iterations, saved_iteration_spacing, electrode_voltage1, electrode_voltage1_ontime, electrode_voltage2, electrode_voltage2_ontime, dt)   # iterations=0:Int64(1E2)  , saved_iteration_spacing=1E0,  current_density=? ,  dt=5e-4
   println(" ");println(" ");println(" ");println(" ")
   println("Starting Pulse Voltage Mode.")

   start_time=Dates.format(Dates.now(),"yyyymmddHHMMSS")   #start_time

   iterations_saved = Array(saved_iteration_spacing:saved_iteration_spacing:length(iterations))
   if any(iterations_saved.==4)==false;  iterations_saved=cat(4,iterations_saved,dims=1);  end;
   if any(iterations_saved.==3)==false;  iterations_saved=cat(3,iterations_saved,dims=1);  end;
   if any(iterations_saved.==2)==false;  iterations_saved=cat(2,iterations_saved,dims=1);  end;
   if any(iterations_saved.==1)==false;  iterations_saved=cat(1,iterations_saved,dims=1);  end;
   if length(iterations) != iterations_saved[end]; iterations_saved=cat(iterations_saved,length(iterations),dims=1);  end;

   simdata=pulse_voltage_simulation_data_structure(
      [start_time]                                                                #start_time
      ,["running"]                                                                #stop_time
      ,start_time * "_dictionary_results.jld2"                                    #data_dictionary_name
      ,dt                                                                         #dt
      ,iterations                                                                 #iterations
      ,saved_iteration_spacing                                                    #saved_iteration_spacing
      ,zeros(length(iterations_saved))                                            #electrode_voltage_saved             
      ,zeros(length(iterations_saved),ss.num_x_mps + ss.spike_num_y_mps + 1 )     #overvoltage_saved
      ,zeros(length(iterations_saved),ss.num_x_mps + ss.spike_num_y_mps + 1 )     #current_density_saved    WE COUNT THE CORNER MESHPOINTS TWICE BECAUSE THEY HAVE INTERFACE POINTING IN THE DY DIRECTION AND THE DX DIRECTION
      ,zeros(length(iterations_saved))                                            #superficial_cd_saved
      ,iterations_saved                                                           #iterations_saved
      ,iterations_saved*dt                                                        #time_real_saved     
      ,zeros(length(iterations_saved))                                            #main_loop_iteration_saved  
      ,zeros(length(iterations_saved),ss.num_x_mps + ss.spike_num_y_mps + 1)      #Charge_Passed_saved       WE COUNT THE CORNER MESHPOINT TWICE BECAUSE THEY HAVE INTERFACE POINTING IN THE DY DIRECTION AND THE DX DIRECTION
      ,zeros(length(iterations_saved),ss.num_y_mps,ss.num_x_mps)                  #conc_A_saved
      ,zeros(length(iterations_saved),ss.num_x_mps + ss.spike_num_y_mps + 1 )     #conc_A_along_surface
   )

   short_y_num = ss.num_y_mps - ss.spike_num_y_mps
   short_x_num = ss.num_x_mps - ss.spike_num_x_mps

   conc_A_along_surface = zeros(ss.num_x_mps + ss.spike_num_y_mps + 1)  #WE COUNT THE CORNER MESHPOINTS TWICE BECAUSE THEY HAVE INTERFACE POINTING IN THE DY DIRECTION AND THE DX DIRECTION
   conc_A_along_surface[1:ss.spike_num_x_mps] = ss.conc_A[short_y_num+1,1:ss.spike_num_x_mps]
   conc_A_along_surface[ss.spike_num_x_mps+1:ss.spike_num_x_mps+ss.spike_num_y_mps] = ss.conc_A[short_y_num+1:ss.num_y_mps,ss.spike_num_x_mps]
   conc_A_along_surface[ss.spike_num_x_mps+ss.spike_num_y_mps+1:end] = ss.conc_A[ss.num_y_mps,ss.spike_num_x_mps:ss.num_x_mps]
   conc_B_along_surface = ss.total_conc .- conc_A_along_surface
   voltage_eq_along_surface = V_eq.(conc_A_along_surface, conc_B_along_surface, ss.conc_A[1,50], ss.total_conc - ss.conc_A[1,50])
   ss.electrode_voltage[1] = electrode_voltage1
   overvoltage = ss.electrode_voltage[1] .- voltage_eq_along_surface
   current_density = -96500*ss.reaction_k .* sqrt.(conc_A_along_surface[:].*conc_B_along_surface[:]).* ( exp.(-(1.0 .- ss.Beta)*96500/8.3/300 .*overvoltage ) .-  exp.(ss.Beta*300/8.3/300 .*overvoltage) )  #(A/m2)
   current_density_from_conc_A = 1.0*current_density
   Charge_Passed = 0.0*current_density
   molar_flux = current_density/96500.0

   record_pulse_voltage_output(ss, simdata, 1, 1, current_density, Charge_Passed, overvoltage, conc_A_along_surface, ss.electrode_voltage[1])

   print("simulated duration is ",simdata.dt[1]*simdata.iterations[end], " seconds\n") 
   println("data saved to produced_data/"*simdata.data_dictionary_name) 
   println("D dt/dx^2 is ",ss.Diffusivity*simdata.dt[1]/ss.dx/ss.dx, " and must be less than 0.5")

   conc_next_timestep = 1.0*ss.conc_A

   #### This for loop increments time
   ## It's an EXPLICIT simulation. It uses a forward Euler time marching scheme. 
   for main_loop_iteration in 2:length(simdata.iterations)

      r_x = ss.Diffusivity * simdata.dt / ss.dx / ss.dx
      r_y = ss.Diffusivity * simdata.dt / ss.dy / ss.dy

      voltage_eq_along_surface[:] = V_eq.(conc_A_along_surface, conc_B_along_surface, ss.conc_A[1,50], ss.total_conc - ss.conc_A[1,50])
      if mod(main_loop_iteration*simdata.dt, electrode_voltage1_ontime+electrode_voltage2_ontime) <= electrode_voltage1_ontime
         ss.electrode_voltage[1] = electrode_voltage1
      else
         ss.electrode_voltage[1] = electrode_voltage2
      end
      overvoltage[:] = ss.electrode_voltage[1] .- voltage_eq_along_surface
      current_density[:] = -96500*ss.reaction_k .* sqrt.(conc_A_along_surface[:].*conc_B_along_surface[:]) .* ( exp.(-(1.0 .- ss.Beta)*96500/8.3/300 .*overvoltage ) .-  exp.(ss.Beta*300/8.3/300 .*overvoltage) )   #(A/m2)
      #println(current_density)
      molar_flux[:] = current_density/96500.0   
      
      # FOR DISCRETIZATION HELP SEE 2014 DISSERTATION BY ZANGANA !!
      # First I calculate the increase in concentration due to y-direction gradients
      for i_x in 1:ss.spike_num_x_mps - 1
         conc_next_timestep[2:short_y_num, i_x]                        =  r_y*ss.conc_A[1:short_y_num - 1, i_x]                      + (1 - 2*r_y)*ss.conc_A[2:short_y_num, i_x]                          +   r_y*ss.conc_A[3:short_y_num+1 , i_x]
         conc_next_timestep[short_y_num+1, i_x]                        =2*r_y*ss.conc_A[short_y_num      , i_x]                      + (1 - 2*r_y)*ss.conc_A[short_y_num+1, i_x]                          + 2*r_y*ss.dy*molar_flux[i_x]/ss.Diffusivity
      end                                           
      for i_x in ss.spike_num_x_mps:ss.num_x_mps                                           
         conc_next_timestep[2:ss.num_y_mps - 1, i_x]                   =  r_y*ss.conc_A[1:ss.num_y_mps - 2, i_x]                     + (1 - 2*r_y)*ss.conc_A[2:ss.num_y_mps - 1, i_x]                     +   r_y*ss.conc_A[3:ss.num_y_mps , i_x]
         conc_next_timestep[ss.num_y_mps      , i_x]                   =2*r_y*ss.conc_A[ss.num_y_mps - 1  , i_x]                     + (1 - 2*r_y)*ss.conc_A[ss.num_y_mps      , i_x]                     + 2*r_y*ss.dy*molar_flux[i_x+ss.spike_num_y_mps+1]/ss.Diffusivity
      end
      #Calculate how much the concentration increases due to y-direction gradients
      conc_increment_dy = conc_next_timestep[:,:] - ss.conc_A[:,:]
      conc_increment_dy[short_y_num+1,ss.spike_num_x_mps] = conc_increment_dy[short_y_num+1,ss.spike_num_x_mps] + 2*r_y*ss.dy*molar_flux[ss.spike_num_x_mps]/ss.Diffusivity

      # Next I calculate the increase in concentration due to x-direction gradients
      for i_y in 1:short_y_num+1
         conc_next_timestep[i_y, 1]                                    =  r_x*ss.conc_A[i_y , 2]                                     + (1 - 2*r_x)*ss.conc_A[i_y , 1]                                     +   r_x*ss.conc_A[i_y , 2]     #the first term is due to the "periodic" boundary condition
         conc_next_timestep[i_y, end]                                  =  r_x*ss.conc_A[i_y , end-1]                                 + (1 - 2*r_x)*ss.conc_A[i_y , end]                                   +   r_x*ss.conc_A[i_y , end-1] #the first term is due to the "periodic" boundary condition
         conc_next_timestep[i_y, 2:ss.num_x_mps - 1]                   =  r_x*ss.conc_A[i_y , 1:ss.num_x_mps - 2]                    + (1 - 2*r_x)*ss.conc_A[i_y , 2:ss.num_x_mps - 1]                    +   r_x*ss.conc_A[i_y , 3:ss.num_x_mps ]
      end
      for i_y in short_y_num+2:ss.num_y_mps
         conc_next_timestep[i_y , ss.spike_num_x_mps]                  =2*r_x*ss.conc_A[i_y , ss.spike_num_x_mps + 1]                + (1 - 2*r_x)*ss.conc_A[i_y , ss.spike_num_x_mps      ]              + 2*r_x*ss.dx*molar_flux[ss.spike_num_x_mps+i_y-short_y_num]/ss.Diffusivity
         conc_next_timestep[i_y , ss.spike_num_x_mps+1:ss.num_x_mps-1] =  r_x*ss.conc_A[i_y , ss.spike_num_x_mps:ss.num_x_mps - 2]   + (1 - 2*r_x)*ss.conc_A[i_y , ss.spike_num_x_mps+1:ss.num_x_mps-1]   +   r_x*ss.conc_A[i_y , ss.spike_num_x_mps+2:ss.num_x_mps] #this last term is due to the "periodic" boundary condition
         conc_next_timestep[i_y , ss.num_x_mps]                        =  r_x*ss.conc_A[i_y , ss.num_x_mps - 1]                      + (1 - 2*r_x)*ss.conc_A[i_y , ss.num_x_mps-1]                        +   r_x*ss.conc_A[i_y , ss.num_x_mps - 1]
      end
      #Calculate how much the concentration increases due to y-direction gradients
      conc_increment_dx = conc_next_timestep[:,:] - ss.conc_A[:,:]
      conc_increment_dx[short_y_num+1,ss.spike_num_x_mps] = conc_increment_dx[short_y_num+1,ss.spike_num_x_mps] + 2*r_y*ss.dy*molar_flux[ss.spike_num_x_mps+1]/ss.Diffusivity

      # Finally, I add the increase in concentration due to x-direction and y-direction gradients together!
      ss.conc_A[:,:] = ss.conc_A[:,:] + conc_increment_dy[:,:] + conc_increment_dx[:,:]

      #Now I check if concentrations dropped below zero and then calculate the real current density  (the real current density on the corner of the spike will be ~2x higher than the Butler-Volmer boundary condition, which was calculated 40 lines above, because the concentration of the corner meshpoint gets reduced by the x and y direction boundary conditions simultaneously)
      conc_A_along_surface[1:ss.spike_num_x_mps] = ss.conc_A[short_y_num+1,1:ss.spike_num_x_mps]
      conc_A_along_surface[ss.spike_num_x_mps+1:ss.spike_num_x_mps+ss.spike_num_y_mps] = ss.conc_A[short_y_num+1:ss.num_y_mps,ss.spike_num_x_mps]
      conc_A_along_surface[ss.spike_num_x_mps+ss.spike_num_y_mps+1:end] = ss.conc_A[ss.num_y_mps,ss.spike_num_x_mps:ss.num_x_mps]
      conc_B_along_surface[:] = ss.total_conc .- conc_A_along_surface
      indices_below_zero = -1000.0 .< conc_A_along_surface[:] .< 0.001
      ss.conc_A[ -1000.0 .< ss.conc_A[:] .< 0.001] .= 0.001
      if length(indices_below_zero) > 0
         current_density_from_conc_A[1:ss.spike_num_x_mps]                                       = -96500*ss.Diffusivity*(ss.conc_A[short_y_num+0,1:ss.spike_num_x_mps]               - ss.conc_A[short_y_num+1,1:ss.spike_num_x_mps]            )/ss.dy
         current_density_from_conc_A[ss.spike_num_x_mps+1:ss.spike_num_x_mps+ss.spike_num_y_mps] = -96500*ss.Diffusivity*(ss.conc_A[short_y_num+1:ss.num_y_mps,ss.spike_num_x_mps+1]  - ss.conc_A[short_y_num+1:ss.num_y_mps,ss.spike_num_x_mps] )/ss.dx
         current_density_from_conc_A[ss.spike_num_x_mps+ss.spike_num_y_mps+1:end]                = -96500*ss.Diffusivity*(ss.conc_A[ss.num_y_mps-1,ss.spike_num_x_mps:ss.num_x_mps]   - ss.conc_A[ss.num_y_mps,ss.spike_num_x_mps:ss.num_x_mps]  )/ss.dy
         current_density[indices_below_zero .== 1] = current_density_from_conc_A[indices_below_zero .== 1]  # Only fix the current density when it drops below zero.  My dx and dy values are 200nm, thus the first ~100 time steps (of 1E-5 seconds) do not allow the concentration boundary layer to penetrate far enough into the first grid space to create the correct value of current density based on F*D*delta_conc/dx, thus I only want to base the current density on F*D*delta_conc/dx when I absolutely have to (when conc drops below zero)
         conc_A_along_surface[1:ss.spike_num_x_mps] = ss.conc_A[short_y_num+1,1:ss.spike_num_x_mps]
         conc_A_along_surface[ss.spike_num_x_mps+1:ss.spike_num_x_mps+ss.spike_num_y_mps] = ss.conc_A[short_y_num+1:ss.num_y_mps,ss.spike_num_x_mps]
         conc_A_along_surface[ss.spike_num_x_mps+ss.spike_num_y_mps+1:end] = ss.conc_A[ss.num_y_mps,ss.spike_num_x_mps:ss.num_x_mps]
         conc_B_along_surface[:] = ss.total_conc .- conc_A_along_surface
      end


      Charge_Passed[:] = Charge_Passed[:] + current_density*simdata.dt

      ##### Save data for post-analysis
      for k in findall(simdata.iterations_saved.==main_loop_iteration)
         @printf(":%-4i   real_time:%+0.7e \n", main_loop_iteration , main_loop_iteration*simdata.dt[1] )
         record_pulse_voltage_output(ss, simdata, k, main_loop_iteration, current_density, Charge_Passed, overvoltage, conc_A_along_surface, ss.electrode_voltage[1])
      end

   end ## this is the end of the for loop of time steps

   ### Save the results
   simdata.stop_time[1]=Dates.format(Dates.now(),"yyyymmddHHMMSS")
   save("produced_data/"*simdata.data_dictionary_name, Dict("system_state"=>ss,"simdata"=>simdata))


   ### Print some basic output
   print("simulated duration is ",simulation_duration, " seconds\n")
   println("using FileIO")
   println("simdata = get(load(\"produced_data/"*simdata.data_dictionary_name*"\"), \"simdata\",0);")
   println("ss =       get(load(\"produced_data/"*simdata.data_dictionary_name*"\"), \"system_state\",0);")
   println("ss,simdata = Diffusion_2D.load_previous_state(\""*simdata.data_dictionary_name[1:14]*"\");")

   return(ss,simdata)

end  ## end of function run_simulation





####### A function to record a time series of the battery state during the CV operation
function record_pulse_voltage_output(ss, simdata, k, main_loop_iteration, current_density, Charge_Passed, overvoltage, conc_A_along_surface, electrode_voltage)
   simdata.time_real_saved[k]          = simdata.iterations[main_loop_iteration]*simdata.dt[1]
   simdata.current_density_saved[k,:]  = current_density[:]
   simdata.superficial_cd_saved[k]     = mean(current_density[:])*length(current_density[:])/ss.num_x_mps  #superficial area to real area
   simdata.Charge_Passed_saved[k,:]    = Charge_Passed[:]   # coulombs per m3
   simdata.electrode_voltage_saved[k]  = electrode_voltage
   simdata.overvoltage_saved[k,:]      = overvoltage[:]
   simdata.conc_A_saved[k,:,:]           = ss.conc_A[:,:]
   simdata.conc_A_along_surface[k,:]   = conc_A_along_surface[:]
end






