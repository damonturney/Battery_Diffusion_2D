############# Functions to plot results from   ###################################################
############# Al:AlCl3:EMImCl:graphite battery ###################################################
using LaTeXStrings, Dates, Printf
import Statistics.mean
import PyPlot

 
function plot_results_time_slice(ss,simdata,k,concentration_color_map_minimum,current_density_colormap_range)
   xs = collect(range(0.0, 4*ss.num_x_mps*ss.dx, length=4*ss.num_x_mps))*1E6
   ys = reverse(ss.locations_y[:,1])*1E6
   sup_fac = length(simdata.current_density_saved[1,:])/ss.num_x_mps  # the multiplicative factor to convert between real surface area and superficial surface area
   conc_data_cat = cat(simdata.conc_A_saved[k,:,:],reverse(simdata.conc_A_saved[k,:,:],dims=2),simdata.conc_A_saved[k,:,:],reverse(simdata.conc_A_saved[k,:,:],dims=2);dims=2)/1000 #to convert from moles/m3 to moles/L
   plot_figure_han = PyPlot.figure(figsize=((4*ss.num_x_mps+100)/100,(ss.num_y_mps+100)/100))
   PyPlot.rc("font", size=9)
   plot_axis_han = plot_figure_han.add_axes([0.12,0.29,0.84,0.7])
   ref_electrode_axis_han = plot_figure_han.add_axes([0.00,0.345,0.25,0.025])
   contour_cb_axes_han = plot_figure_han.add_axes([0.05,0.18,0.9,0.03])
   scatt_cb_axes_han   = plot_figure_han.add_axes([0.05,0.07,0.9,0.03])
   plot_axis_han.set_facecolor((0.5, 0.5, 0.5))
   #fig_x, fig_y, fig_dx, fig_dy = PyPlot.get_current_fig_manager().window.geometry().getRect()  #left border location , top border location, width, height     
   #PyPlot.get_current_fig_manager().window.setGeometry(0,fig_y,fig_dx,fig_dy) 
   levels = collect(range(concentration_color_map_minimum,stop=1,length=50))
   contourf_han = plot_axis_han.contourf(xs,ys,conc_data_cat,levels,cmap="plasma", vmin=concentration_color_map_minimum, vmax=1.0, zorder=1)
   contour_han = plot_axis_han.contour(xs,ys,conc_data_cat,levels,linewidths=0.5, colors="k", vmin=concentration_color_map_minimum, vmax=1, alpha=0.25, zorder=2)
   ylim=plot_axis_han.get_ylim()
   xlim=plot_axis_han.get_xlim()
   contour_cb_han = PyPlot.colorbar(contourf_han, cax=contour_cb_axes_han, orientation="horizontal", ticks=collect(range(concentration_color_map_minimum,step=0.1,stop=1.0)))
   contour_cb_han.set_label("concentration (moles/L)", labelpad=1)
   plot_axis_han.set_xlabel("microns", labelpad=1)
   plot_axis_han.set_ylabel("microns", labelpad=-1)
   #plot_axis_han.text(0.1,0.85,                 " time         ave. current density      overvoltage")
   #plot_axis_han.text(0.1,0.15,Printf.@sprintf("%3.2fs, %8.0f mAh/cm2, %8i mV", simdata.time_real_saved[k], abs(mean(simdata.current_density_saved[k,:])/10) , abs(simdata.electrode_voltage_saved[k]*1000)),family="monospace")
   plot_axis_han.annotate("working\nelectrode",                                                                   (0.50,0.93), horizontalalignment="center", xycoords="axes fraction", fontsize=9, weight="bold")
   plot_axis_han.annotate("voltage"                                                                           ,   (0.50,0.82), horizontalalignment="center", xycoords="axes fraction", fontsize=8)
   plot_axis_han.annotate(Printf.@sprintf("%i mV" ,      simdata.electrode_voltage_saved[k]*1000),                (0.56,0.79), horizontalalignment="right",  xycoords="axes fraction", fontsize=8)
   plot_axis_han.annotate("space ave.",                                                                           (0.17,0.96), horizontalalignment="center", xycoords="axes fraction", fontsize=8)
   plot_axis_han.annotate(Printf.@sprintf("%4.2f mA/cm", mean(simdata.current_density_saved[k,:])/10)*L"^2",      (0.29,0.93), horizontalalignment="right", xycoords="axes fraction", fontsize=8)
   plot_axis_han.annotate("time-space ave.",                                                                      (0.85,0.96), horizontalalignment="center", xycoords="axes fraction", fontsize=8)
   plot_axis_han.annotate(Printf.@sprintf("%4.2f mA/cm", simdata.superficial_cd_time_average[1]/sup_fac/10)*L"^2",(0.98,0.93), horizontalalignment="right", xycoords="axes fraction", fontsize=8)
   plot_axis_han.annotate("time"                                                              ,                   (0.88,0.04), horizontalalignment="left"  , xycoords="axes fraction", fontsize=8)
   plot_axis_han.annotate(Printf.@sprintf("%2.2fs",     simdata.time_saved[k] ),                                  (0.99,0.01), horizontalalignment="right"  , xycoords="axes fraction", fontsize=8)
   plot_axis_han.annotate("superficial",                                                                          (0.63,0.04), horizontalalignment="center", xycoords="axes fraction", fontsize=8)
   plot_axis_han.annotate(Printf.@sprintf("%4.2f mA/cm", simdata.superficial_cd_saved[k]/10 )*L"^2",              (0.77,0.01), horizontalalignment="right",  xycoords="axes fraction", fontsize=8)
   plot_axis_han.annotate("t.s. ave. superficial",                                                                (0.19,0.04), horizontalalignment="center", xycoords="axes fraction", fontsize=8)
   plot_axis_han.annotate(Printf.@sprintf("%4.2f mA/cm", simdata.superficial_cd_time_average[1]/10 )*L"^2",       (0.34,0.01), horizontalalignment="right",  xycoords="axes fraction", fontsize=8)
   ref_electrode_axis_han.set_facecolor((0.5, 0.5, 0.5))
   ref_electrode_axis_han.get_yaxis().set_visible(false); ref_electrode_axis_han.get_xaxis().set_visible(false); ref_electrode_axis_han.spines["left"].set_visible(false); ref_electrode_axis_han.spines["right"].set_visible(false); ref_electrode_axis_han.spines["top"].set_visible(false); ref_electrode_axis_han.spines["bottom"].set_visible(false);
   ref_electrode_axis_han.annotate("ref electrode", (0.04, 0.2), xycoords="axes fraction", weight="bold")
   surface_xs = [collect(1:ss.spike_num_x_mps)*ss.dx; repeat([ss.spike_num_x_mps*ss.dx],inner=ss.spike_num_y_mps) ; collect(ss.spike_num_x_mps:ss.num_x_mps)*ss.dx ; ]*1E6
   surface_ys = ys[[repeat([ss.num_y_mps - ss.spike_num_y_mps],inner=ss.spike_num_x_mps) ; collect(ss.num_y_mps - ss.spike_num_y_mps:ss.num_y_mps) ; repeat([ss.num_y_mps],inner=ss.num_x_mps-ss.spike_num_x_mps)]]
   surface_xs4 = [surface_xs ; -reverse(surface_xs) .+ 2*surface_xs[end] ; surface_xs .+ 2*surface_xs[end] ; -reverse(surface_xs) .+ 4*surface_xs[end] ] #Since the lhs and rhs are periodic boundaries... I simply repeat the results 4 times side by side
   surface_ys4 = [surface_ys ; reverse(surface_ys) ; surface_ys ; reverse(surface_ys)]                                                                   #Since the lhs and rhs are periodic boundaries... I simply repeat the results 4 times side by side
   current_density4 = [simdata.current_density_saved[k,:] ; reverse(simdata.current_density_saved[k,:]) ; simdata.current_density_saved[k,:] ; reverse(simdata.current_density_saved[k,:])]/10 #divide by 10 to convert to mA/cm2     #Since the lhs and rhs are periodic boundaries... I simply repeat the results 4 times side by side
   
   # Make a user-defined colormap.   see https://matplotlib.org/3.1.0/tutorials/colors/colormap-manipulation.html and https://github.com/JuliaPy/PyPlot.jl
   bwr_cmap = PyPlot.cm.get_cmap("bwr", 256)
   colors = bwr_cmap(collect(range(0, step=1, stop=255)));
   newcolors = 1.0*colors
   newcolors[:,3]=colors[:,2]; newcolors[:,2]=colors[:,3];
   newcmp = PyPlot.ColorMap("gwr",newcolors);

   #Scatter the current density onto the surface as color circles
   scatt_han = plot_axis_han.scatter(surface_xs4, surface_ys4, c=current_density4, s=7,  edgecolors="none", cmap=newcmp, vmin=-current_density_colormap_range, vmax=current_density_colormap_range, zorder=3) 
   plot_axis_han.set_xlim(xlim)
   plot_axis_han.set_ylim([0,ys[end,1]*1.1])
   scatter_cb_han = PyPlot.colorbar(scatt_han, cax=scatt_cb_axes_han, orientation="horizontal", alpha=0.5, ticks=collect(range(-current_density_colormap_range,step=current_density_colormap_range/4,stop=current_density_colormap_range)))
   scatter_cb_han.set_label("current density "*L"(mA/cm^2)", labelpad=1)
end




function make_movie(ss,simdata,start_frame_num,end_frame_num,concentration_color_map_minimum, current_density_colormap_range)
   rm("produced_data/video_images_"*simdata.data_dictionary_name[1:14]*"/",recursive=true,force=true)
   rm("produced_data/"*simdata.data_dictionary_name[1:14]*"_movie.mp4",force=true)
   mkdir("produced_data/video_images_"*simdata.data_dictionary_name[1:14]*"/")
   PyPlot.close("all")
   image_num=0
   for frame_num in start_frame_num:end_frame_num
      image_num=image_num+1
      println("frame num: "*string(frame_num))
      plot_results_time_slice(ss,simdata,frame_num,concentration_color_map_minimum, current_density_colormap_range)
      PyPlot.savefig("produced_data/video_images_"*simdata.data_dictionary_name[1:14]*"/"*Printf.@sprintf("%04d", image_num)*".png", dpi=300)
      PyPlot.close("all") 
   end
   println("frame num: "*string(end_frame_num))
   plot_results_time_slice(ss,simdata,end_frame_num,concentration_color_map_minimum, current_density_colormap_range)
   PyPlot.savefig("produced_data/video_images_"*simdata.data_dictionary_name[1:14]*"/"*Printf.@sprintf("%04d", image_num+1)*".png", dpi=300)
   PyPlot.close("all") 
   bash_command ="ffmpeg -start_number 1 -i produced_data/video_images_"*simdata.data_dictionary_name[1:14]*"/%04d.png -vcodec libx265 -x265-params \"lossless=1\" -preset slow -vf format=yuv420p produced_data/"*simdata.data_dictionary_name[1:14]*"_movie.mp4"
   run(`bash -c $bash_command`)
   rm("produced_data/video_images_"*simdata.data_dictionary_name[1:14]*"/",recursive=true,force=true)
end


function plot_results_time_series(ss,simdata)
   a=1
end


