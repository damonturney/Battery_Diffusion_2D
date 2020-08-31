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
   superficial_cd_saved            ::Array{Float64,1}
   iterations_saved                ::Array{Float64,1}
   time_real_saved                 ::Array{Float64,1}
   main_loop_iteration_saved       ::Array{Float64,1}
   Charge_Passed_saved             ::Array{Float64,2}
   conc_A_saved                    ::Array{Float64,3}
   conc_A_along_surface            ::Array{Float64,2}
end




######## A function to create and initialize a single instance of the data structure 
function pulse_current(ss, iterations, saved_iteration_spacing, superficial_current_density_target1, current_density1_ontime, superficial_current_density_target2, current_density2_ontime, dt)   # iterations=0:Int64(1E2)  , saved_iteration_spacing=1E0,  current_density=? ,  dt=5e-4
   println(" ");println(" ");println(" ");println(" ")
   println("Starting Pulse Current Mode.")

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

   #Calculate the equilibrium voltage along the surface
   conc_A_along_surface = zeros(ss.num_x_mps + ss.spike_num_y_mps + 1)  #WE COUNT THE CORNER MESHPOINTS TWICE BECAUSE THEY HAVE INTERFACE POINTING IN THE DY DIRECTION AND THE DX DIRECTION
   conc_A_along_surface[1:ss.spike_num_x_mps]                                       = ss.conc_A[short_y_num+1,1:ss.spike_num_x_mps]
   conc_A_along_surface[ss.spike_num_x_mps+1:ss.spike_num_x_mps+ss.spike_num_y_mps] = ss.conc_A[short_y_num+1:ss.num_y_mps,ss.spike_num_x_mps]
   conc_A_along_surface[ss.spike_num_x_mps+ss.spike_num_y_mps+1:end]                = ss.conc_A[ss.num_y_mps,ss.spike_num_x_mps:ss.num_x_mps]
   conc_B_along_surface                                                             = ss.total_conc .- conc_A_along_surface
   voltage_eq_along_surface = V_eq.(conc_A_along_surface, conc_B_along_surface, ss.conc_A[1,50], ss.total_conc - ss.conc_A[1,50])
   
   superficial_current_density_target = 1.0*superficial_current_density_target1
   overvoltage = ss.electrode_voltage[1] .- voltage_eq_along_surface[:]
   current_density = -96500*ss.reaction_k .*abs.(overvoltage[:]).* sqrt.(conc_A_along_surface[:].*conc_B_along_surface[:]).* ( exp.(-(1.0 .- ss.Beta)*96500/8.3/300 .*overvoltage ) .-  exp.(ss.Beta*300/8.3/300 .*overvoltage) )  #(A/m2)
   current_density_from_conc_A = 1.0*current_density
   superficial_current_density = mean(current_density[:])*length(current_density[:])/ss.num_x_mps
   superficial_current_density_error = superficial_current_density - superficial_current_density_target
   while abs(superficial_current_density_error) > 0.1
      overvoltage = ss.electrode_voltage[1] .- voltage_eq_along_surface
      current_density = -96500*ss.reaction_k .*abs.(overvoltage[:]).* sqrt.(conc_A_along_surface[:].*conc_B_along_surface[:]).* ( exp.(-(1.0 .- ss.Beta)*96500/8.3/300 .*overvoltage ) .-  exp.(ss.Beta*300/8.3/300 .*overvoltage) )  #(A/m2)
      superficial_current_density = mean(current_density[:])*length(current_density[:])/ss.num_x_mps
      superficial_current_density_error = superficial_current_density - superficial_current_density_target
      ss.electrode_voltage[1] = ss.electrode_voltage[1] - superficial_current_density_error*1E-7
      #@printf("superficial_current_density:%+0.7e   electrode_voltage:%+0.7e     superficial_cd_error:%0.7e \n", superficial_current_density , ss.electrode_voltage[1],  superficial_current_density_error)
      sleep(0.5)
   end
   Charge_Passed = 0.0*current_density
   molar_flux = current_density/96500.0 
   
   record_pulse_current_output(ss, sim_data, 1, 1, current_density, superficial_current_density, Charge_Passed, overvoltage, conc_A_along_surface, ss.electrode_voltage[1])

   print("simulated duration is ",sim_data.dt*sim_data.iterations[end], " seconds\n") 
   println("data saved to produced_data/"*sim_data.data_dictionary_name) 
   println("D dt/dx^2 is ",ss.Diffusivity*sim_data.dt/ss.dx/ss.dx, " and must be less than 0.5")

   conc_next_timestep = 1.0*ss.conc_A
   superficial_current_density_target_previous           = -1E6
   superficial_current_density_target_previous_previous  = -1E6
   electrode_voltage_previous_previous                   = -1E6
   electrode_voltage_previous                            = -1E6
   conc_A_eq_along_surface_previous                      = 1.0*conc_A_eq_along_surface


   ########### This for loop increments time.   It's an EXPLICIT simulation. It uses a forward Euler time marching scheme. 
   for main_loop_iteration in 2:length(sim_data.iterations)

      r_x = ss.Diffusivity * sim_data.dt / ss.dx / ss.dx
      r_y = ss.Diffusivity * sim_data.dt / ss.dy / ss.dy

      voltage_eq_along_surface[:] = V_eq.(conc_A_along_surface, conc_B_along_surface, ss.conc_A[1,50], ss.total_conc - ss.conc_A[1,50])

      #Update the current density target
      superficial_current_density_target_previous_previous = 1.0*superficial_current_density_target_previous
      superficial_current_density_target_previous          = 1.0*superficial_current_density_target
      if mod(main_loop_iteration*sim_data.dt, current_density1_ontime + current_density2_ontime) <= current_density1_ontime
         superficial_current_density_target = 1.0*superficial_current_density_target1
      else
         superficial_current_density_target = 1.0*superficial_current_density_target2
      end
      
      # Calculate the electrode voltage by predictor-corrector
      electrode_voltage_previous_previous = 1.0*electrode_voltage_previous
      electrode_voltage_previous          = 1.0*ss.electrode_voltage[1]
      if superficial_current_density_target == superficial_current_density_target_previous_previous
         ss.electrode_voltage[1] = electrode_voltage_previous + (electrode_voltage_previous - electrode_voltage_previous_previous)
      end
      overvoltage[:] = ss.electrode_voltage[1] .- voltage_eq_along_surface
      current_density = -96500*ss.reaction_k .*abs.(overvoltage[:]).*sqrt.(conc_A_along_surface[:].*conc_B_along_surface[:]).* ( exp.(-(1.0 .- ss.Beta)*96500/8.3/300 .*overvoltage[:] ) .-  exp.(ss.Beta*300/8.3/300 .*overvoltage[:]) )  #(A/m2)
      superficial_current_density = mean(current_density[:])*length(current_density[:])/ss.num_x_mps
      superficial_current_density_error = superficial_current_density - superficial_current_density_target
      while abs(superficial_current_density_error) > 0.1
         overvoltage[:] = ss.electrode_voltage[1] .- voltage_eq_along_surface
         current_density[:] = -96500*ss.reaction_k .*abs.(overvoltage[:]).*sqrt.(conc_A_along_surface[:].*conc_B_along_surface[:]).* ( exp.(-(1.0 .- ss.Beta)*96500/8.3/300 .*overvoltage[:] ) .-  exp.(ss.Beta*300/8.3/300 .*overvoltage[:]) )  #(A/m2)
         superficial_current_density       = mean(current_density[:])*length(current_density[:])/ss.num_x_mps
         superficial_current_density_error = superficial_current_density - superficial_current_density_target
         ss.electrode_voltage[1] = ss.electrode_voltage[1] - superficial_current_density_error*1E-7
         #@printf("loop:%5.0i   superficial_current_density:%+0.7e   electrode_voltage:%+0.7e   target_cd:%+0.7e   superficial_cd_error:%+0.7e\n", main_loop_iteration, superficial_current_density , ss.electrode_voltage[1] , superficial_current_density_target , superficial_current_density_error)
         #sleep(0.5)
      end
      #@printf("loop.%5.0i   superficial_current_density:%+0.7e   electrode_voltage:%+0.7e   target_cd:%+0.7e   superficial_cd_error:%+0.7e\n", main_loop_iteration, superficial_current_density , ss.electrode_voltage[1] , superficial_current_density_target , superficial_current_density_error)
      
      # Calculate the equilibrium conc_A, and find out which mesh points have concentration above the equilibrium value
      d_V_eq_d_conc_A = ( V_eq.(conc_A_along_surface[50] .+ 0.1, conc_B_along_surface[50] .- 0.1, ss.conc_A[1,50], ss.total_conc - ss.conc_A[1,50]) .- voltage_eq_along_surface[50] ) / 0.1   
      d_conc_A_d_V_eq = 1 / d_V_eq_d_conc_A
      conc_A_eq_along_surface = conc_A_along_surface[50] .+ d_conc_A_d_V_eq .* overvoltage[50]
      conc_A_eq_along_surface_previous[:] = conc_A_eq_along_surface[:]
      indices_above_eq = conc_A_along_surface .> conc_A_eq_along_surface
      for k in findall(sim_data.iterations_saved.==main_loop_iteration)
         #println(ss.conc_A[end-2,29:37])
         #println(ss.conc_A[end-1,29:37])
         #println(ss.conc_A[end-0,29:37])
         @printf("loop:%5.0i   conc_A_corner:%+0.7e   conc_A_eq:%+0.7e\n", main_loop_iteration, ss.conc_A[end,30], conc_A_eq_along_surface)
      end


      molar_flux[:] = current_density[:]/96500.0   
      
      # FOR DISCRETIZATION HELP SEE 2014 DISSERTATION BY ZANGANA !! IF there’s a flux boundary condition at x = 0 then: du/dt = (2*D*u_i+1,j - 2*D*u_i,j )/deltax^2 - 2*D*g/deltax,  u_i,j+1 = u_i,j + 2r*u_i+1,j - 2r*u_i,j - 2r*g*deltax   where g is du/dx on the boundary at x =0 and the 2 comes from the fact that we divide by 0.5deltax, not full deltax (because of boundary at x=0).  If there's a LHS "half boundary" due to a corner on the spike then:  du/dt = ( D*( u_i+1,j - u_i,j )/deltax - (0.5*D*g + 0.5*D*( u_i,j - u_i-1,j )/deltax) )/0.75/deltax,  du/dt = 4/3*D*( u_i+1,j - u_i,j )/deltax^2 - 2/3*D*( u_i,j - u_i-1,j )/deltax^2 - 2/3*D*g/deltax ,   u_i,j+1 = u_i,j + 4/3*r*( u_i+1,j - u_i,j ) - 2/3*r*( u_i,j - u_i-1,j ) - 2/3*r*g*deltax = 2/3*r*u_i-1,j + (1 - 2*r)*u_i,j + 4/3*r*u_i+1,j - 2/3*r*g*deltax       and if the boundary condition is a RHS "half boundary" then this is: u_i,j+1 = 4/3*r*u_i-1,j + (1 - 2*r)*u_i,j + 2/3*r*u_i+1,j + 2/3*r*g*deltax
      # First I calculate the increase in concentration due to y-direction gradients
      for i_x in 1:ss.spike_num_x_mps - 1
         conc_next_timestep[2:short_y_num, i_x]                        =   r_y*ss.conc_A[1:short_y_num - 1, i_x]                      + (1 - 2*r_y)*ss.conc_A[2:short_y_num, i_x]                          +    r_y*ss.conc_A[3:short_y_num+1 , i_x]
         conc_next_timestep[short_y_num+1, i_x]                        = 2*r_y*ss.conc_A[short_y_num      , i_x]                      + (1 - 2*r_y)*ss.conc_A[short_y_num+1, i_x]                                                                                         +   2*r_y*ss.dy*molar_flux[i_x]/ss.Diffusivity
      end
      i_x = ss.spike_num_x_mps
      conc_next_timestep[2:short_y_num, i_x]                           =   r_y*ss.conc_A[1:short_y_num - 1, i_x]                      + (1 - 2*r_y)*ss.conc_A[2:short_y_num, i_x]                          +    r_y*ss.conc_A[3:short_y_num+1, i_x]
      conc_next_timestep[short_y_num+1, i_x]                           =4/3*r_y*ss.conc_A[short_y_num, i_x]                           + (1 - 2*r_y)*ss.conc_A[short_y_num+1, i_x]                          +2/3*r_y*ss.conc_A[short_y_num+2, i_x]                         + 2/3*r_y*ss.dy*molar_flux[ss.spike_num_x_mps]/ss.Diffusivity  #This line is for the corner mesh point.  The 0.5 factor is because this meshpoint only has electrode interface on half it's area.
      conc_next_timestep[short_y_num+2:ss.num_y_mps-1,  i_x]           =   r_y*ss.conc_A[short_y_num+1:ss.num_y_mps-2, i_x]           + (1 - 2*r_y)*ss.conc_A[short_y_num+2:ss.num_y_mps-1, i_x]           +    r_y*ss.conc_A[short_y_num+3:ss.num_y_mps , i_x]
      conc_next_timestep[ss.num_y_mps,  i_x]                           = 2*r_y*ss.conc_A[ss.num_y_mps - 1, i_x]                       + (1 - 2*r_y)*ss.conc_A[ss.num_y_mps, i_x]                                                                                          +   2*r_y*ss.dy*molar_flux[ss.spike_num_x_mps + ss.spike_num_y_mps+1]/ss.Diffusivity  #This line is for the corner mesh point.  The 0.5 factor is because this meshpoint only has electrode interface on half it's area.
      for i_x in ss.spike_num_x_mps+1:ss.num_x_mps                                           
         conc_next_timestep[2:ss.num_y_mps - 1, i_x]                   =   r_y*ss.conc_A[1:ss.num_y_mps - 2, i_x]                     + (1 - 2*r_y)*ss.conc_A[2:ss.num_y_mps - 1, i_x]                     +    r_y*ss.conc_A[3:ss.num_y_mps , i_x]
         conc_next_timestep[ss.num_y_mps      , i_x]                   =2*r_y*ss.conc_A[ss.num_y_mps - 1  , i_x]                      + (1 - 2*r_y)*ss.conc_A[ss.num_y_mps, i_x]                                                                                          +   2*r_y*ss.dy*molar_flux[i_x + ss.spike_num_y_mps + 1]/ss.Diffusivity
      end
      #Calculate how much the concentration increases due to y-direction gradients         
      conc_increment_dy = conc_next_timestep[:,:] - ss.conc_A[:,:]       #This line is for all mesh points 

      # Next I calculate the increase in concentration due to x-direction gradients
      for i_y in 1:short_y_num
         conc_next_timestep[i_y, 1]                                    =   r_x*ss.conc_A[i_y , 2]                                     + (1 - 2*r_x)*ss.conc_A[i_y , 1]                                     +    r_x*ss.conc_A[i_y , 2]                                    #this is due to the "periodic" boundary condition
         conc_next_timestep[i_y, end]                                  =   r_x*ss.conc_A[i_y , end-1]                                 + (1 - 2*r_x)*ss.conc_A[i_y , end]                                   +    r_x*ss.conc_A[i_y , end-1]                                #this is due to the "periodic" boundary condition
         conc_next_timestep[i_y, 2:ss.num_x_mps - 1]                   =   r_x*ss.conc_A[i_y , 1:ss.num_x_mps - 2]                    + (1 - 2*r_x)*ss.conc_A[i_y , 2:ss.num_x_mps - 1]                    +    r_x*ss.conc_A[i_y , 3:ss.num_x_mps ]
      end
      i_y = short_y_num+1
      conc_next_timestep[i_y, 1]                                       =   r_x*ss.conc_A[i_y , 2]                                     + (1 - 2*r_x)*ss.conc_A[i_y , 1]                                     +    r_x*ss.conc_A[i_y , 2]                                    #this is due to the "periodic" boundary condition
      conc_next_timestep[i_y, 2:ss.spike_num_x_mps-1]                  =   r_x*ss.conc_A[i_y , 1:ss.spike_num_x_mps-2]                + (1 - 2*r_x)*ss.conc_A[i_y , 2:ss.spike_num_x_mps - 1]              +    r_x*ss.conc_A[i_y , 3:ss.spike_num_x_mps ]  
      conc_next_timestep[i_y, ss.spike_num_x_mps]                      =2/3*r_x*ss.conc_A[i_y , ss.spike_num_x_mps-1]                 + (1 - 2*r_x)*ss.conc_A[i_y , ss.spike_num_x_mps]                    +4/3*r_x*ss.conc_A[i_y , ss.spike_num_x_mps+1]                 - 2/3*r_x*ss.dx*(-molar_flux[ss.spike_num_x_mps+1]/ss.Diffusivity)
      conc_next_timestep[i_y, ss.spike_num_x_mps+1:ss.num_x_mps-1]     =   r_x*ss.conc_A[i_y , ss.spike_num_x_mps:ss.num_x_mps-2]     + (1 - 2*r_x)*ss.conc_A[i_y , ss.spike_num_x_mps+1:ss.num_x_mps-1]   +    r_x*ss.conc_A[i_y , ss.spike_num_x_mps+2:ss.num_x_mps ]  
      conc_next_timestep[i_y, end]                                     =   r_x*ss.conc_A[i_y , end-1]                                 + (1 - 2*r_x)*ss.conc_A[i_y , end]                                   +    r_x*ss.conc_A[i_y , end-1]                                #this is due to the "periodic" boundary condition
      for i_y in short_y_num+2:ss.num_y_mps-1
         conc_next_timestep[i_y , ss.spike_num_x_mps]                  = 2*r_x*ss.conc_A[i_y , ss.spike_num_x_mps + 1]                + (1 - 2*r_x)*ss.conc_A[i_y , ss.spike_num_x_mps]                                                                                   -  2*r_x*ss.dx*(-molar_flux[ss.spike_num_x_mps + i_y - short_y_num]/ss.Diffusivity)
         conc_next_timestep[i_y , ss.spike_num_x_mps+1:ss.num_x_mps-1] =   r_x*ss.conc_A[i_y , ss.spike_num_x_mps:ss.num_x_mps - 2]   + (1 - 2*r_x)*ss.conc_A[i_y , ss.spike_num_x_mps+1:ss.num_x_mps-1]   +    r_x*ss.conc_A[i_y , ss.spike_num_x_mps+2:ss.num_x_mps] 
         conc_next_timestep[i_y , ss.num_x_mps]                        =   r_x*ss.conc_A[i_y , ss.num_x_mps - 1]                      + (1 - 2*r_x)*ss.conc_A[i_y , ss.num_x_mps]                          +    r_x*ss.conc_A[i_y , ss.num_x_mps - 1]                     #this is due to the "periodic" boundary condition
      end
      i_y = ss.num_y_mps
      conc_next_timestep[i_y , ss.spike_num_x_mps]                     = 2*r_x*ss.conc_A[i_y , ss.spike_num_x_mps+1]                  + (1 - 2*r_x)*ss.conc_A[i_y , ss.spike_num_x_mps]                                                                                   -  2*r_x*ss.dx*(-molar_flux[ss.spike_num_x_mps + ss.spike_num_y_mps + 1]/ss.Diffusivity)
      conc_next_timestep[i_y , ss.spike_num_x_mps+1:ss.num_x_mps-1]    =   r_x*ss.conc_A[i_y , ss.spike_num_x_mps:ss.num_x_mps - 2]   + (1 - 2*r_x)*ss.conc_A[i_y , ss.spike_num_x_mps+1:ss.num_x_mps-1]   +    r_x*ss.conc_A[i_y , ss.spike_num_x_mps+2:ss.num_x_mps]   
      conc_next_timestep[i_y , ss.num_x_mps]                           =   r_x*ss.conc_A[i_y , ss.num_x_mps - 1]                      + (1 - 2*r_x)*ss.conc_A[i_y , ss.num_x_mps]                          +    r_x*ss.conc_A[i_y , ss.num_x_mps - 1]                     #this is due to the "periodic" boundary condition
      #Calculate how much the concentration increases due to y-direction gradients
      conc_increment_dx = conc_next_timestep[:,:] - ss.conc_A[:,:]

      # Finally, I add the increase in concentration due to x-direction and y-direction gradients together!
      ss.conc_A[:,:] = ss.conc_A[:,:] + conc_increment_dy[:,:] + conc_increment_dx[:,:]

      #Now I check if concentrations crossed the equilibrium concentration due to the discretized Fickian diffusion. If so, set conc_A equal to the equilibrium, and calculate current density (the real current density on the corner of the spike will be ~2x higher than the Butler-Volmer boundary condition, which was calculated 40 lines above, because the concentration of the corner meshpoint gets reduced by the x and y direction boundary conditions simultaneously)
      conc_A_along_surface[1:ss.spike_num_x_mps]                                       = ss.conc_A[short_y_num+1,1:ss.spike_num_x_mps]
      conc_A_along_surface[ss.spike_num_x_mps+1:ss.spike_num_x_mps+ss.spike_num_y_mps] = ss.conc_A[short_y_num+1:ss.num_y_mps,ss.spike_num_x_mps]
      conc_A_along_surface[ss.spike_num_x_mps+ss.spike_num_y_mps+1:end]                = ss.conc_A[ss.num_y_mps,ss.spike_num_x_mps:ss.num_x_mps]
      conc_B_along_surface[:]                                                          = ss.total_conc .- conc_A_along_surface
      indices_below_eq = conc_A_along_surface .< conc_A_eq_along_surface
      indices_crossed_eq = (indices_below_eq .& indices_above_eq) .| (.~indices_below_eq .& .~indices_above_eq)
      if length(findall(indices_crossed_eq)) > 0
         ##### The line below is not physics-based.  You need to come up with a way to fix this problem based on physical equations.
         conc_A_along_surface[indices_crossed_eq] .= conc_A_eq_along_surface_previous[indices_crossed_eq] .+ 0.85*( conc_A_eq_along_surface .- conc_A_eq_along_surface_previous[indices_crossed_eq] ) 
         ss.conc_A[short_y_num+1,1:ss.spike_num_x_mps]            = conc_A_along_surface[1:ss.spike_num_x_mps]
         ss.conc_A[short_y_num+1:ss.num_y_mps,ss.spike_num_x_mps] = conc_A_along_surface[ss.spike_num_x_mps+1:ss.spike_num_x_mps+ss.spike_num_y_mps]
         ss.conc_A[ss.num_y_mps,ss.spike_num_x_mps:ss.num_x_mps]  = conc_A_along_surface[ss.spike_num_x_mps+ss.spike_num_y_mps+1:end]
         conc_B_along_surface[:]                                  = ss.total_conc .- conc_A_along_surface
      end

      Charge_Passed[:] = Charge_Passed[:] + current_density*sim_data.dt

      ##### Save data for post-analysis
      for k in findall(sim_data.iterations_saved.==main_loop_iteration)
         #println(ss.conc_A[end-2,29:37])
         #println(ss.conc_A[end-1,29:37])
         #println(ss.conc_A[end-0,29:37])
         #@printf("loop:%5.0i   conc_A_corner:%+0.7e   conc_A_eq:%+0.7e\n", main_loop_iteration, ss.conc_A[end,30], conc_A_eq_along_surface)
         #println(conc_A_along_surface[ss.spike_num_x_mps + ss.spike_num_y_mps - 1 : ss.spike_num_x_mps + ss.spike_num_y_mps + 2 ])
         @printf(":%-4i   real_time:%+0.7e \n", main_loop_iteration , main_loop_iteration*sim_data.dt )
         record_pulse_current_output(ss, sim_data, k, main_loop_iteration, current_density, superficial_current_density, Charge_Passed, overvoltage, conc_A_along_surface, ss.electrode_voltage[1])
         @printf("loop.%5.0i   superficial_current_density:%+0.7e   electrode_voltage:%+0.7e   target_cd:%+0.7e   superficial_cd_error:%+0.7e\n", main_loop_iteration, superficial_current_density , ss.electrode_voltage[1] , superficial_current_density_target , superficial_current_density_error)
      end

   end ## this is the end of the for loop of time steps

   ### Save the results
   sim_data.stop_time[1]=Dates.format(Dates.now(),"yyyymmddHHMMSS")
   save("produced_data/"*sim_data.data_dictionary_name, Dict("system_state"=>ss,"sim_data"=>sim_data))


   ### Print some very basic output
   print("simulated duration is ",sim_data.iterations[end]*sim_data.dt, " seconds\n")
   println("using FileIO")
   println("sim_data = get(load(\"produced_data/"*sim_data.data_dictionary_name*"\"), \"sim_data\",0);")
   println("ss =       get(load(\"produced_data/"*sim_data.data_dictionary_name*"\"), \"system_state\",0);")
   println("Number of saved timesteps: "*string(length(sim_data.iterations_saved)))

   return(ss,sim_data)

end  ## end of function run_simulation




####### A function to record a time series of the battery state during the CV operation
function record_pulse_current_output(ss, sim_data, k, main_loop_iteration, current_density, superficial_current_density, Charge_Passed, overvoltage, conc_A_along_surface, electrode_voltage)
      sim_data.time_real_saved[k]          = sim_data.iterations[main_loop_iteration]*sim_data.dt
      sim_data.current_density_saved[k,:]  = current_density[:]
      sim_data.superficial_cd_saved[k]     = superficial_current_density  #superficial area to real area
      sim_data.Charge_Passed_saved[k,:]    = Charge_Passed[:]   # coulombs per m3
      sim_data.electrode_voltage_saved[k]  = electrode_voltage
      sim_data.overvoltage_saved[k,:]      = overvoltage[:]
      sim_data.conc_A_saved[k,:,:]         = ss.conc_A[:,:]
      sim_data.conc_A_along_surface[k,:]   = conc_A_along_surface[:]
end


