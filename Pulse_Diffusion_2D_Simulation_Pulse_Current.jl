####################################################################################################
###### Pulse Current Operator
####################################################################################################


###### Create the data type that will hold all information for the pulse current operator (including time series)
struct pulse_current_simulation_data_structure
   start_time                      ::Array{String,1}
   stop_time                       ::Array{String,1}
   data_dictionary_name            ::String
   dt                              ::Float64
   iterations                      ::UnitRange{Int64}
   saved_iteration_spacing         ::Float64
   electrode_voltage_saved         ::Array{Float64,1}
   overvoltage_saved               ::Array{Float64,2}
   current_density_saved           ::Array{Float64,2}
   iterations_saved                ::Array{Float64,1}
   time_real_saved                 ::Array{Float64,1}
   main_loop_iteration_saved       ::Array{Float64,1}
   Charge_Passed_saved             ::Array{Float64,2}
   conc_saved                      ::Array{Float64,3}
   conc_A_along_surface            ::Array{Float64,2}
end




######## A function to create and initialize a single instance of the data structure 
function pulse_current(ss, iterations, saved_iteration_spacing, current_density_target1, current_density1_ontime, current_density_target2, current_density2_ontime, dt)   # iterations=0:Int64(1E2)  , saved_iteration_spacing=1E0,  current_density=? ,  dt=5e-4
   println(" ");println(" ");println(" ");println(" ")
   println("Creating CyclicV_operator.")

   start_time=Dates.format(Dates.now(),"yyyymmddHHMMSS")   #start_time

   iterations_saved = Array(saved_iteration_spacing:saved_iteration_spacing:length(iterations))
   if any(iterations_saved.==4)==false;  iterations_saved=cat(4,iterations_saved,dims=1);  end;
   if any(iterations_saved.==3)==false;  iterations_saved=cat(3,iterations_saved,dims=1);  end;
   if any(iterations_saved.==2)==false;  iterations_saved=cat(2,iterations_saved,dims=1);  end;
   if any(iterations_saved.==1)==false;  iterations_saved=cat(1,iterations_saved,dims=1);  end;
   if length(iterations) != iterations_saved[end]; iterations_saved=cat(iterations_saved,length(iterations),dims=1);  end;

   sim_data=pulse_current_simulation_data_structure(
      [start_time]                                                                #start_time
      ,["running"]                                                                #stop_time
      ,start_time * "_dictionary_results.jld2"                                    #data_dictionary_name
      ,dt                                                                         #dt
      ,iterations                                                                 #iterations
      ,saved_iteration_spacing                                                    #saved_iteration_spacing
      ,zeros(length(iterations_saved))                                            #electrode_voltage_saved             
      ,zeros(length(iterations_saved),ss.num_x_mps + ss.spike_num_y_mps + 1 )     #overvoltage_saved
      ,zeros(length(iterations_saved),ss.num_x_mps + ss.spike_num_y_mps + 1 )     #current_density_saved    WE COUNT THE CORNER MESHPOINTS TWICE BECAUSE THEY HAVE INTERFACE POINTING IN THE DY DIRECTION AND THE DX DIRECTION
      ,iterations_saved                                                           #iterations_saved
      ,iterations_saved*dt                                                        #time_real_saved     
      ,zeros(length(iterations_saved))                                            #main_loop_iteration_saved  
      ,zeros(length(iterations_saved),ss.num_x_mps + ss.spike_num_y_mps + 1)      #Charge_Passed_saved       WE COUNT THE CORNER MESHPOINT TWICE BECAUSE THEY HAVE INTERFACE POINTING IN THE DY DIRECTION AND THE DX DIRECTION
      ,zeros(length(iterations_saved),ss.num_y_mps,ss.num_x_mps)                  #conc_saved
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
   current_density_target =current_density_target1
   overvoltage = ss.electrode_voltage[1] .- voltage_eq_along_surface
   current_density = -96500*ss.reaction_k .* sqrt.(conc_A_along_surface[:].*conc_B_along_surface[:]).* ( exp.(-(1.0 .- ss.Beta)*96500/8.3/300 .*overvoltage ) .-  exp.(ss.Beta*300/8.3/300 .*overvoltage) )  #(A/m2)
   current_density_error = mean(current_density) - current_density_target
   while abs(current_density_error) > 0.1
      overvoltage = ss.electrode_voltage[1] .- voltage_eq_along_surface
      current_density = -96500*ss.reaction_k .* sqrt.(conc_A_along_surface[:].*conc_B_along_surface[:]).* ( exp.(-(1.0 .- ss.Beta)*96500/8.3/300 .*overvoltage ) .-  exp.(ss.Beta*300/8.3/300 .*overvoltage) )  #(A/m2)
      current_density_error = mean(current_density) - current_density_target
      ss.electrode_voltage[1] = ss.electrode_voltage[1] - mean(current_density_error)*1E-7
      #@printf("mean_current_density:%+0.7e   electrode_voltage:%+0.7e     cd_error:%0.7e \n", mean(current_density) , ss.electrode_voltage[1],  current_density_error)
      #sleep(0.5)
   end
   Charge_Passed = 0.0*current_density
   molar_flux = current_density/96500.0

   record_pulse_current_output(ss, sim_data, 1, 1, current_density, Charge_Passed, overvoltage, conc_A_along_surface, ss.electrode_voltage[1])

   print("simulated duration is ",sim_data.dt[1]*sim_data.iterations[end], " seconds\n") 
   println("data saved to produced_data/"*sim_data.data_dictionary_name) 
   println("D dt/dx^2 is ",ss.Diffusivity*sim_data.dt[1]/ss.dx/ss.dx, " and must be less than 0.5")

   conc_next_timestep = 0.0*ss.conc_A

   #### This for loop increments time
   ## It's an EXPLICIT simulation. It uses a forward Euler time marching scheme. 
   for main_loop_iteration in 2:length(sim_data.iterations)

      r_x = ss.Diffusivity * sim_data.dt / ss.dx / ss.dx
      r_y = ss.Diffusivity * sim_data.dt / ss.dy / ss.dy
      conc_next_timestep[1,:] = 1.0*ss.conc_A[1,:] 

      conc_A_along_surface[1:ss.spike_num_x_mps] = ss.conc_A[short_y_num+1,1:ss.spike_num_x_mps]
      conc_A_along_surface[ss.spike_num_x_mps+1:ss.spike_num_x_mps+ss.spike_num_y_mps] = ss.conc_A[short_y_num+1:ss.num_y_mps,ss.spike_num_x_mps]
      conc_A_along_surface[ss.spike_num_x_mps+ss.spike_num_y_mps+1:end] = ss.conc_A[ss.num_y_mps,ss.spike_num_x_mps:ss.num_x_mps]
      conc_B_along_surface[:] = ss.total_conc .- conc_A_along_surface
      #println(conc_A_along_surface/ss.conc_A[1,50])
      voltage_eq_along_surface[:] = V_eq.(conc_A_along_surface, conc_B_along_surface, ss.conc_A[1,50], ss.total_conc - ss.conc_A[1,50])
      electrode_voltage1, electrode_voltage1_ontime, electrode_voltage2, electrode_voltage2_ontime
      if mod(main_loop_iteration*sim_data.dt, current_density1_ontime + current_density2_ontime) <= current_density1_ontime
         current_density_target = current_density_target1
      else
         current_density_target = current_density_target2
      end
      overvoltage[:] = ss.electrode_voltage[1] .- voltage_eq_along_surface
      current_density = -96500*ss.reaction_k .* sqrt.(conc_A_along_surface[:].*conc_B_along_surface[:]).* ( exp.(-(1.0 .- ss.Beta)*96500/8.3/300 .*overvoltage ) .-  exp.(ss.Beta*300/8.3/300 .*overvoltage) )  #(A/m2)
      current_density_error = mean(current_density) - current_density_target
      while abs(current_density_error) > 0.1
         overvoltage = ss.electrode_voltage[1] .- voltage_eq_along_surface
         current_density = -96500*ss.reaction_k .* sqrt.(conc_A_along_surface[:].*conc_B_along_surface[:]).* ( exp.(-(1.0 .- ss.Beta)*96500/8.3/300 .*overvoltage ) .-  exp.(ss.Beta*300/8.3/300 .*overvoltage) )  #(A/m2)
         current_density_error = mean(current_density) - current_density_target
         ss.electrode_voltage[1] = ss.electrode_voltage[1] - current_density_error*1E-7
         #@printf("loop:%5.0i   mean_current_density:%+0.7e   electrode_voltage:%+0.7e \n", main_loop_iteration, mean(current_density) , ss.electrode_voltage[1] )
         #sleep(0.5)
      end
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

      # Finally, I add the increase in concentration due to x-direction and y-direction gradients together!
      conc_increment_dx = conc_next_timestep[:,:] - ss.conc_A[:,:]
      conc_increment_dx[short_y_num+1,ss.spike_num_x_mps] = conc_increment_dx[short_y_num+1,ss.spike_num_x_mps] + 2*r_y*ss.dy*molar_flux[ss.spike_num_x_mps+1]/ss.Diffusivity

      ss.conc_A[:,:] = ss.conc_A[:,:] + conc_increment_dy[:,:] + conc_increment_dx[:,:]

      Charge_Passed[:] = Charge_Passed[:] + current_density*sim_data.dt

      ##### Save data for post-analysis
      for k in findall(sim_data.iterations_saved.==main_loop_iteration)
         @printf(":%-4i   real_time:%+0.7e \n", main_loop_iteration , main_loop_iteration*sim_data.dt[1] )
         record_pulse_current_output(ss, sim_data, k, main_loop_iteration, current_density, Charge_Passed, overvoltage, conc_A_along_surface, ss.electrode_voltage[1])
      end

   end ## this is the end of the for loop of time steps

   ### Save the results
   sim_data.stop_time[1]=Dates.format(Dates.now(),"yyyymmddHHMMSS")
   save("produced_data/"*sim_data.data_dictionary_name, Dict("system_state"=>ss,"sim_data"=>sim_data))


   ### Print some very basic output
   print("simulated duration is ",sim_data.iterations[end]*sim_data.dt[1], " seconds\n")
   println("using FileIO")
   println("sim_data = get(load(\"produced_data/"*sim_data.data_dictionary_name*"\"), \"sim_data\",0);")
   println("ss =       get(load(\"produced_data/"*sim_data.data_dictionary_name*"\"), \"system_state\",0);")
   println("Number of saved timesteps: "*string(length(sim_data.iterations_saved)))

   return(ss,sim_data)

end  ## end of function run_simulation




####### A function to record a time series of the battery state during the CV operation
function record_pulse_current_output(ss, sim_data, k, main_loop_iteration, current_density, Charge_Passed, overvoltage, conc_A_along_surface, electrode_voltage)
      sim_data.time_real_saved[k]          = sim_data.iterations[main_loop_iteration]*sim_data.dt[1]
      sim_data.current_density_saved[k,:]  = current_density[:]
      sim_data.Charge_Passed_saved[k,:]    = Charge_Passed[:]   # coulombs per m3
      sim_data.electrode_voltage_saved[k]  = electrode_voltage
      sim_data.overvoltage_saved[k,:]      = overvoltage[:]
      sim_data.conc_saved[k,:,:]           = ss.conc_A[:,:]
      sim_data.conc_A_along_surface[k,:]   = conc_A_along_surface[:]
end

