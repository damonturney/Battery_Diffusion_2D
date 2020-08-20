############# Functions to plot results from   ###################################################
############# Al:AlCl3:EMImCl:graphite battery ###################################################
using LaTeXStrings, Dates, Printf
import Statistics.mean
import PyPlot

 
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
   #fig_x, fig_y, fig_dx, fig_dy = PyPlot.get_current_fig_manager().window.geometry().getRect()  #left border location , top border location, width, height     
   #PyPlot.get_current_fig_manager().window.setGeometry(0,fig_y,fig_dx,fig_dy) 
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




function make_movie(ss,sim_data,end_frame_num)
   rm("produced_data/video_images/",recursive=true,force=true)
   rm("produced_data/"*sim_data.data_dictionary_name[1:14]*"_movie.mp4",force=true)
   mkdir("produced_data/video_images/")
   PyPlot.close("all")
   for frame_num in 1:end_frame_num
      println("frame num: "*string(frame_num))
      plot_results_time_slice(ss,sim_data,frame_num)
      PyPlot.savefig("produced_data/video_images/"*Printf.@sprintf("%04d", frame_num)*".png", dpi=300)
      PyPlot.close("all") 
   end
   bash_command ="ffmpeg -start_number 4 -i produced_data/video_images/%04d.png -vcodec libx265 -x265-params \"lossless=1\" -preset slow -vf format=yuv420p produced_data/"*sim_data.data_dictionary_name[1:14]*"_movie.mp4"
   run(`bash -c $bash_command`)
end


function plot_results_time_series(ss,sim_data)
   a=1
end


