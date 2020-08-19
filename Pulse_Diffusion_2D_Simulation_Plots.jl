############# Functions to plot results from   ###################################################
############# Al:AlCl3:EMImCl:graphite battery ###################################################
using Plots, LaTeXStrings, Plots.PlotMeasures, Dates, Printf
import Statistics.mean
import PyPlot

 
function plot_results_time_slice(ss,sim_data,k)
   xs = collect(range(0.0, 4*ss.num_x_mps*ss.dx, length=4*ss.num_x_mps))*1E6
   ys = reverse(ss.locations_y[:,1])*1E6
   conc_data_cat = cat(sim_data.conc_saved[k,:,:],reverse(sim_data.conc_saved[k,:,:],dims=2),sim_data.conc_saved[k,:,:],reverse(sim_data.conc_saved[k,:,:],dims=2);dims=2)
   plot_figure_han = PyPlot.figure(figsize=((4*ss.num_x_mps+300)/100,(ss.num_y_mps+100)/100))
   #plot_axis_han=plot_figure_han.add_subplot(1,1,1)
   plot_axis_han = plot_figure_han.add_axes([0.1,0.4,0.7,0.57])
   contour_cb_axes_han = plot_figure_han.add_axes([0.1,0.2,0.7,0.08])
   scatt_cb_axes_han   = plot_figure_han.add_axes([0.1,0.05,0.7,0.08])
   plot_axis_han.set_facecolor((0.5, 0.5, 0.5))
   fig_x, fig_y, fig_dx, fig_dy = PyPlot.get_current_fig_manager().window.geometry().getRect()  #left border location , top border location, width, height     
   PyPlot.get_current_fig_manager().window.setGeometry(0,fig_y,fig_dx,fig_dy) 
   levels = collect(range(600,stop=1000,length=50))
   contourf_han = plot_axis_han.contourf(xs,ys,conc_data_cat,levels,cmap="plasma", zorder=1)
   contour_han = plot_axis_han.contour(xs,ys,conc_data_cat,levels,linewidths=0.5, colors="k", alpha=0.25, zorder=2)
   ylim=plot_axis_han.get_ylim()
   xlim=plot_axis_han.get_xlim()
   contour_cb_han = PyPlot.colorbar(contourf_han, cax=contour_cb_axes_han, orientation="horizontal")
   plot_axis_han.set_xlabel("microns")
   plot_axis_han.set_ylabel("microns")
   surface_xs1 = [collect(1:ss.spike_num_x_mps)*ss.dx; repeat([ss.spike_num_x_mps*ss.dx],inner=ss.spike_num_y_mps) ; collect(ss.spike_num_x_mps:ss.num_x_mps)*ss.dx ; ]*1E6
   surface_ys1 = ys[[repeat([ss.num_y_mps - ss.spike_num_y_mps],inner=ss.spike_num_x_mps) ; collect(ss.num_y_mps - ss.spike_num_y_mps:ss.num_y_mps) ; repeat([ss.num_y_mps],inner=ss.num_x_mps-ss.spike_num_x_mps)]]
   surface_xs4 = [surface_xs1 ; -reverse(surface_xs1) .+ 2*surface_xs1[end] ; surface_xs1 .+ 2*surface_xs1[end] ; -reverse(surface_xs1) .+ 4*surface_xs1[end] ]
   surface_ys4 = [surface_ys1 ; reverse(surface_ys1) ; surface_ys1 ; reverse(surface_ys1)]
   Charge_Passed4 = [sim_data.Charge_Passed ; reverse(sim_data.Charge_Passed) ; sim_data.Charge_Passed ; reverse(sim_data.Charge_Passed)]
   #for i in 1:length(surface_xs)
   #   sleep(0.1)
   #   plot_axis_han.text(surface_xs[i], surface_ys[i], c="b")#-sim_data.Charge_Passed[i])
   #end
   scatt_han = plot_axis_han.scatter(surface_xs4, surface_ys4, c=Charge_Passed4/minimum(Charge_Passed4)*255, s=7, cmap="winter", zorder=3) 
   plot_axis_han.set_xlim(xlim)
   plot_axis_han.set_ylim([0,ys[end,1]*1.1])
   scatter_cb_han = PyPlot.colorbar(scatt_han, cax=scatt_cb_axes_han, orientation="horizontal")
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

