############# Functions to plot results from   ###################################################
############# Al:AlCl3:EMImCl:graphite battery ###################################################
using Plots, LaTeXStrings, Plots.PlotMeasures, Dates, Printf
import Statistics.mean
import PyPlot
import Printf
using LaTeXStrings

 
function plot_results_time_slice(ss,sim_data,k)
   xs = collect(range(0.0, 4*ss.num_x_mps*ss.dx, length=4*ss.num_x_mps))*1E6
   ys = reverse(ss.locations_y[:,1])*1E6
   conc_data_cat = cat(sim_data.conc_saved[k,:,:],reverse(sim_data.conc_saved[k,:,:],dims=2),sim_data.conc_saved[k,:,:],reverse(sim_data.conc_saved[k,:,:],dims=2);dims=2)/1000 #to convert from moles/m3 to moles/L
   plot_figure_han = PyPlot.figure(figsize=((4*ss.num_x_mps+100)/100,(ss.num_y_mps+100)/100))
   PyPlot.rc("font", size=9)
   plot_axis_han = plot_figure_han.add_axes([0.12,0.29,0.84,0.7])
   contour_cb_axes_han = plot_figure_han.add_axes([0.05,0.18,0.9,0.03])
   scatt_cb_axes_han   = plot_figure_han.add_axes([0.05,0.07,0.9,0.03])
   plot_axis_han.set_facecolor((0.5, 0.5, 0.5))
   fig_x, fig_y, fig_dx, fig_dy = PyPlot.get_current_fig_manager().window.geometry().getRect()  #left border location , top border location, width, height     
   PyPlot.get_current_fig_manager().window.setGeometry(0,fig_y,fig_dx,fig_dy) 
   levels = collect(range(0.6,stop=1,length=50))
   contourf_han = plot_axis_han.contourf(xs,ys,conc_data_cat,levels,cmap="plasma", vmin=0.6, vmax=1, zorder=1)
   contour_han = plot_axis_han.contour(xs,ys,conc_data_cat,levels,linewidths=0.5, colors="k", vmin=0.6, vmax=1, alpha=0.25, zorder=2)
   ylim=plot_axis_han.get_ylim()
   xlim=plot_axis_han.get_xlim()
   contour_cb_han = PyPlot.colorbar(contourf_han, cax=contour_cb_axes_han, orientation="horizontal", ticks=[0.6,0.7,0.8,0.9,1.0])
   contour_cb_han.set_label("concentration (moles/L)", labelpad=1)
   plot_axis_han.set_xlabel("microns", labelpad=1)
   plot_axis_han.set_ylabel("microns", labelpad=1)
   #plot_axis_han.text(0.1,0.85,                 " time         ave. current density      overvoltage")
   #plot_axis_han.text(0.1,0.15,Printf.@sprintf("%3.2fs, %8.0f mAh/cm2, %8i mV", sim_data.time_real_saved[k], abs(mean(sim_data.current_density_saved[k,:])/10) , abs(sim_data.electrode_voltage_saved[k]*1000)),family="monospace")
   plot_axis_han.annotate("electrode", (0.5,0.95), xycoords="axes fraction", horizontalalignment="center")
   plot_axis_han.annotate(" time        ave. current density      overvoltage",                               (4, 11), xycoords="axes points")
   plot_axis_han.annotate(Printf.@sprintf("%2.2f s",     sim_data.time_real_saved[k] ),                       (35, 1), horizontalalignment="right", xycoords="axes points")
   plot_axis_han.annotate(Printf.@sprintf("%4.0f mA/cm", mean(sim_data.current_density_saved[k,:])/10)*L"^2", (130,1), horizontalalignment="right", xycoords="axes points")
   plot_axis_han.annotate(Printf.@sprintf("%i mV" ,      sim_data.electrode_voltage_saved[k]*1000),           (200,1), horizontalalignment="right", xycoords="axes points")
   surface_xs1 = [collect(1:ss.spike_num_x_mps)*ss.dx; repeat([ss.spike_num_x_mps*ss.dx],inner=ss.spike_num_y_mps) ; collect(ss.spike_num_x_mps:ss.num_x_mps)*ss.dx ; ]*1E6
   surface_ys1 = ys[[repeat([ss.num_y_mps - ss.spike_num_y_mps],inner=ss.spike_num_x_mps) ; collect(ss.num_y_mps - ss.spike_num_y_mps:ss.num_y_mps) ; repeat([ss.num_y_mps],inner=ss.num_x_mps-ss.spike_num_x_mps)]]
   surface_xs4 = [surface_xs1 ; -reverse(surface_xs1) .+ 2*surface_xs1[end] ; surface_xs1 .+ 2*surface_xs1[end] ; -reverse(surface_xs1) .+ 4*surface_xs1[end] ]
   surface_ys4 = [surface_ys1 ; reverse(surface_ys1) ; surface_ys1 ; reverse(surface_ys1)]
   current_density4 = [sim_data.current_density_saved[k,:] ; reverse(sim_data.current_density_saved[k,:]) ; sim_data.current_density_saved[k,:] ; reverse(sim_data.current_density_saved[k,:])]/10 #divide by 10 to convert to mA/cm2
   scatt_han = plot_axis_han.scatter(surface_xs4, surface_ys4, c=abs.(current_density4), s=7, cmap="winter", edgecolors="none", vmin=0, vmax=250, zorder=3) 
   plot_axis_han.set_xlim(xlim)
   plot_axis_han.set_ylim([0,ys[end,1]*1.1])
   scatter_cb_han = PyPlot.colorbar(scatt_han, cax=scatt_cb_axes_han, orientation="horizontal", ticks=[0,50,100,150,200,250])
   scatter_cb_han.set_label("current density "*L"(mA/cm^2)", labelpad=1)
end




function make_movie(ss,sim_data)
   rm("produced_data/video_images/",recursive=true,force=true)
   rm("produced_data/"*sim_data.data_dictionary_name[1:14]*"_movie.mp4")
   mkdir("produced_data/video_images/")
   PyPlot.close("all")
   for frame_num in 1:30
      plot_results_time_slice(ss,sim_data,frame_num)
      PyPlot.savefig("produced_data/video_images/"*Printf.@sprintf("%04d", frame_num)*".png", dpi=300)
      PyPlot.close("all") 
   end
   bash_command ="ffmpeg -start_number 4 -i produced_data/video_images/%04d.png -vcodec libx265 -x265-params \"lossless=1\" -preset slow -vf format=yuv420p produced_data/"*sim_data.data_dictionary_name[1:14]*"_movie.mp4"
   run(`bash -c $bash_command`)
end


function plot_results_time_series(ss,sim_data)
   plot_folder = ss.system_state_folder_name*"/CV_time_series_"
   num_digits_cell_voltage_plot =  round(Int64,-log10((maximum(sim_data.cell_voltage_saved) - minimum(sim_data.cell_voltage_saved))/10 )) + 2 # The +2 is for the "1." characters in front of everything
   xformatter_cell_voltage=x->if typeof(x)==Float64 @sprintf("%0.10e",x)[[1:num_digits_cell_voltage_plot;]] end
   yformatter=y->if typeof(y)==Float64 @sprintf("%0.10e",y)[[1:4; end-3:end]] end
            plot( sim_data.cell_voltage_saved,sim_data.eta_cathode_saved[:,20],                        markershape=:circle, linealpha=0.0,markersize=5,legend=:topleft,label="eta_cathode")
            plot!(sim_data.cell_voltage_saved,sim_data.V_elctrlt_saved[:,20],                          markershape=:circle, linealpha=0.0,markersize=5,legend=:topleft,label="Voltage elctrlt")
   savefig(plot!(sim_data.cell_voltage_saved,sim_data.eta_anode_saved[:,20],                          markershape=:circle, linealpha=0.0,markersize=5,legend=:topleft,label="eta_Al",   xticks=range(minimum(sim_data.cell_voltage_saved), length=10, stop=maximum(sim_data.cell_voltage_saved)),xrotation=45,xformatter=xformatter_cell_voltage ), plot_folder*"overvoltages.png");
   savefig(plot( sim_data.cell_voltage_saved, sim_data.i_sol_sup_stagrd_saved[:,ss.num_cathode_nodes+10], markershape=:circle, linealpha=0.0,markersize=5,legend=:none,title="I vs cell V", xticks=range(minimum(sim_data.cell_voltage_saved), length=10, stop=maximum(sim_data.cell_voltage_saved)),xrotation=45,xformatter=xformatter_cell_voltage ), plot_folder*"cell_voltage.png");
end

function plot_results_image(ss,sim_data,k)
end


###########To make movies
#theta_animation(1,length(sim_data.cell_voltage_saved),0.0,1.0)
function theta_vs_cv_animation(ss,sim_data,first_time,last_time,thetamin,thetamax)
        anim = @animate for k=first_time:last_time
                layout1 = grid(2,1,heights=[0.7,0.3])
                # plot(  ss.pcl_x_nodes_fine*ss.radius_cathode_pcl*1e6,sim_data.theta_fine_saved[k,:,1],size=(600,400),line=(2,:blue),markeralpha=0,legend=:none,ylims=(thetamin,thetamax),layout=layout1,subplot=1,xaxis = ("Location in Graphite Particle (um)", font(8)),yaxis = ("\\theta (n.d.)", (0,1), font(8)),bottom_margin=0mm,left_margin=-2.5mm)
                # plot!(-ss.pcl_x_nodes_fine*ss.radius_cathode_pcl*1e6,sim_data.theta_fine_saved[k,:,1],size=(600,400),line=(2,:blue),markeralpha=0,legend=:none,ylims=(thetamin,thetamax),layout=layout1,subplot=1,xaxis = ("Location in Graphite Particle (um)", font(8)),yaxis = ("\\theta (n.d.)", (0,1), font(8)),bottom_margin=0mm,left_margin=-2.5mm)
                plot(    ss.pcl_x_nodes[:]  *ss.radius_cathode_pcl*1e6,  sim_data.theta_saved[k,:,10], linealpha=0.0, marker=(:circle,5,:blue, stroke(0)),legend=:none,ylims=(thetamin,thetamax),xaxis = ("Location in Graphite Particle (um)", font(20)),yaxis = ("\\theta (n.d.)", (0,1), font(20)),bottom_margin=0mm,left_margin=2mm,layout=layout1,subplot=1,size=(840,1000))
                plot!(  -ss.pcl_x_nodes[:]  *ss.radius_cathode_pcl*1e6,  sim_data.theta_saved[k,:,10], linealpha=0.0, marker=(:circle,5,:blue, stroke(0)),legend=:none,ylims=(thetamin,thetamax),layout=layout1,subplot=1)
                plot!( [ ss.pcl_x_nodes[end]*ss.radius_cathode_pcl*1e6],[sim_data.theta_eq_intfc_saved[k,10]],        linealpha=0.0, marker=(:circle,6,:green,stroke(0)),legend=:none,ylims=(thetamin,thetamax),layout=layout1,subplot=1)
                plot!( [-ss.pcl_x_nodes[end]*ss.radius_cathode_pcl*1e6],[sim_data.theta_eq_intfc_saved[k,10]],        linealpha=0.0, marker=(:circle,6,:green,stroke(0)),legend=:none,ylims=(thetamin,thetamax),layout=layout1,subplot=1)
                plot!(sim_data.cell_voltage_saved,sim_data.i_sol_sup_stagrd_saved[:,ss.num_cathode_nodes+10],linealpha=0.0,marker=(:circle, 6, stroke(0),:blue),legend=:none,ylims=(minimum(sim_data.i_sol_sup_stagrd_saved[:,ss.num_cathode_nodes+10])*1.1,maximum(sim_data.i_sol_sup_stagrd_saved[:,ss.num_cathode_nodes+10])*1.1),layout=layout1,subplot=2,xaxis = ("Cell Voltage  (V)", font(20)),yaxis = (L"\mathrm{Current \,\,\, (A \,m^{-2})}", font(20)),bottom_margin=0mm,left_margin=0mm)
                plot!([sim_data.cell_voltage_saved[k],sim_data.cell_voltage_saved[k]],[minimum(sim_data.i_sol_sup_stagrd_saved[:,ss.num_cathode_nodes+10])*1.1,maximum(sim_data.i_sol_sup_stagrd_saved[:,ss.num_cathode_nodes+10])*1.1],line=(1,:red),linecolor=:red,layout=layout1,subplot=2)
        end
        gif(anim, "produced_data/" *sim_data.data_dictionary_name[1:14] * "_theta_animation" * ".gif", fps = 15)
end

