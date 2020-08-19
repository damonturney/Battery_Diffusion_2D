


#This function will be used in many location
function create_staggered_linear_spline(input)
   input_diff=diff(input,dims=1)
   temp    =input[1:end-1] + input_diff/2
   temp0   =input[1]       - input_diff[1]/2
   temp_end=input[end]     + input_diff[end]/2
   output=[temp0  temp' temp_end]
   return(output[:])
 end

 
 function smooth_gradient(input)
   output=input.*1.0;
   output[1]=input[2]-input[1];
   output[2]=(input[3]-input[1])/2
   for i in 3:length(output)-2;
       output[i]=(-input[i+2]+8*input[i+1]-8*input[i-1]+input[i-2])/12;
   end
   output[end-1]=(input[end]-input[end-2])/2;
   output[end]=input[end]-input[end-1];
   return(output)
end



function update_conc_staggered2(dx, dt, velocity_staggered, conc_old, reaction_rate_source_sink, porosity_stagrd, porosity_nonstagrd)   #Porosity only enters the reaction_rate term.   It cancels out of the velocity convergence/divergence term.
    array_size=length(conc_old);
    conc_new::Array{BigFloat,1} =  zeros(BigFloat,array_size);
    for i in 2:array_size-1
        conc_new[i] = conc_old[i]   + dt*( porosity_stagrd[i]*velocity_staggered[i]*(conc_old[i-1]+conc_old[i])/2/dx               -  porosity_stagrd[i+1]*velocity_staggered[i+1]*(conc_old[i]+conc_old[i+1])/2/dx  + reaction_rate_source_sink[i]  ) / porosity_nonstagrd[i] #*dx/dx    reaction_rate has units of amps per m^3
    end
    conc_new[1] =     conc_old[1]   + dt*(                                    0                                                    -  porosity_stagrd[2]*velocity_staggered[2]*(conc_old[1]+conc_old[2])/2/dx        + reaction_rate_source_sink[1]  ) / porosity_nonstagrd[1] #*dx/dx       reaction_rate has units of amps per m^3
    conc_new[end] =   conc_old[end] + dt*( porosity_stagrd[end-1]*velocity_staggered[end-1]*(conc_old[end-1]+conc_old[end])/2/dx   -                                       0                                         + reaction_rate_source_sink[end]) / porosity_nonstagrd[end]  #*dx/dx    reaction_rate has units of amps per m^3
                                                           #end-1#
                                                           #yes the above "end-1" should REALLY be "end-1" because the "end" location of velocity_staggered is on the RHS!  you want the LHS!  the velocity_staggered array has more elements than the concentration array!
    return(conc_new)
end

