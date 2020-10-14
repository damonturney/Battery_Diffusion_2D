########## Description of the Physics Strategy and Discretization Strategy ########################################################################################################################

##############################################################################################################################################################################################




######### Verions ##############################################################################################################################################################################################
v1.1
v1.2
v1.3
v1.4
v1.5
v1.6
v1.7
v1.8
#######################################################################################################################################################################################################################



######### Work Flow ##############################################################################################################################################################################################
include("Pulse_Diffusion_2D_Simulation_Main.jl");
ss = Diffusion_2D.create_system_state(); 
ss,simdata = Diffusion_2D.pulse_current(ss , 0.1, 1E-5, 0.01, -200, 10, -200 , 10);  #units are always mks.  pulse_voltage(ss, simulation_duration, dt_biggest, saved_dt_spacing, superficial_current_density_target1, current_density1_ontime, superficial_current_density_target2, current_density2_ontime) 
#Diffusion_2D.plot_results_time_slice(ss,sim_data,5)
Diffusion_2D.make_movie(ss,simdata,505)
#using FileIO
#cv_data = get(load("produced_data/20200803102000/20200803102016_dictionary_results.jld2"), "cyclicV_data",0);
#ss =      get(load("produced_data/20200803102000/20200803102016_dictionary_results.jld2"), "system_state"          ,0);
#ss =      get(load(data_file_name), "system_state",0);
#cv_data = get(load(data_file_name), "cyclicV_data",0);
#Diffusion_2D.plot_results_time_series(ss,cv_data)

########################################################################################################################################################################################################






######### To Do ##################################################################################################################################################################################################################
1) I encountered "overshoot" problems for mesh points along the electrode interface when the concentration gets low (near ~5.0 mole / m^3 or lower) because the electrode voltage starts to change by large amounts in order to maintain 
   the same current density, thus the overvoltages get large and the current remains high, but since concentration is low (near ~5 mole/m3 or lower) it can shoot below zero during a dt step or it can pass the equilibrium concentration 
   (which in reality wouldn't happen because the rate of change of concentration should asymptote to the equilibrium concentration), and also the equilibrium voltage is VERY sensitive to concentration for low concentration thus small 
   inaccuracies in the time-evolution of concentraiton lead to large inaccuracies in the overvoltage calculation and it spirals out of control. You have two options to fix this:
1a) Put the calculations of Fickian time-evolution of concentration at the interfacial mesh points inside the while loop that chooses electrode voltage.  This way we can detect when the concentration overshoots and fix it to a value 
    near it's equilibrium value (using artificial "semi-physical" approximations), and then we can calculate what the true current density and electrode voltage should be.
1b) Don't do 1a but, rather, create dynamic timestep length in which you choose the timestep length to be small enough to disallow overshooting at interfacial locations.  The time duration between changes of the electrode voltage 
    would be held constant at 1E-5 seconds (or something like that) but the timesteps for calculating Fickian changes to concentration would be smaller (like 1E-6 or 1E-7 seconds).
1 conclusion)  I opted to do option 1a but I can alway redo it to be 1b.

2) Read through Julia's performance optimization tips and apply them
##############################################################################################################################################################################################







######### ALPHABETICAL LIST OF VIM COMMANDS ######################################################################################################################################################
Movement
h j k l
Basic movement keys. A step up from the cursor keys simply because they are already under your fingers. Most useful when prefixed with a number (e.g. if you need to move down by about 10 lines, hit “10j” instead of just holding j until you get there).
b w B W
Move back by token/forward by token/back by word/forward by word.  (A token is a sequence of letters, digits, and underscores. For the capital letter variations, a word consists of anything that’s not whitespace.) Faster than holding down a simple directional key.
0 ^ $
Jump to first column/first non-whitespace character/end of line, like Home and End. Faster than moving by words if you’re trying to get to the opposite end of the line.
ctrl+u ctrl+d
Basically Page Up and Page Down, but moves by half a screenful and doesn’t lose your cursor position.
<line-number>G
Jump directly to a specific line number. Most helpful if you also have line numbering enabled (:set number).
H M L
Move to the top/middle/bottom of the screen (i.e. High/Middle/Low). A good first step in getting approximately to where you want to go.
# *
Find the previous/next occurrence of the token under the cursor.
n N
Repeat the last find command forward/backward.
“
(That’s two back-ticks). Jump back to where you just were. This will jump back and forth between the same two locations if you keep pressing it.
ctrl+o ctrl+i
Move backward/forward through the jump history. Useful if you have followed a chain of method calls and need to get back to where you were.

Editing
Most editing commands may optionally be preceded by a number in order to apply it more than once (e.g. to delete three lines, press 3dd).
i a I A
Enter insert mode (insert at cursor/append after cursor/insert at beginning of line/append to end of line). Press Esc to exit insert mode and return to normal mode.  It’s rarely useful to precede one of these commands with a number, but it can come in handy. Need a comma-separated list of eight 1s? Just hit “8i1, <esc>” then delete the trailing comma.
o O
Open new line (below the current line/above the current line). A quick “o<esc>” will add a blank line below the current line, no matter where your cursor is.
cw cW
Correct the token(s)/word(s) following the cursor. Basically combines delete and insert into one step.
cc
Correct line(s) by clearing and then entering insert mode. Starts inserting at the current indent level.
dd
Delete line(s). Quickly rearrange lines by deleting them, moving to the new location, and pasting with “p”.
ct cf ci ca
dt df di da
Correct/delete up to or including specific characters. Since there are many variations, I break it down in the section below about correcting text.
s
Delete character(s) at the cursor and then enter insert mode.  cw is usually faster if you want to correct an entire word, but this is useful for correcting a fixed number of characters (e.g. “5s” will correct the next five characters).
yy
Copy line(s).  The “y” is for “yank.”
yw yW
Copy token(s)/word(s).
p P
Paste the last thing that was deleted or copied before/after cursor (for more advanced usage, you can precede it with a register specification, but that’s a topic for another day).
u ctrl+r
Undo and redo.
.
(That’s a period). Repeat the previous edit command. I use this all the time. Did you just add a line (e.g. using “o” or “O”) that you need to duplicate five more times with only slight modifications?  Hit “5.” to repeat that operation, then make your modifications; no copy/paste needed.

Correcting Text
If you need to correct some text that doesn’t fall neatly onto a token or word boundary. Fortunately, operations like “c” (correct) and “d” (delete) have a number of operators that may be applied to them:
t<char> – exclusive match: continue up to (but not including) the next <char> on this line
f<char> – inclusive match: continue up to (and including) the next <char> on this line
i<char> – exclusive inner match: apply to text bounded by <char>, where <char> is from a limited set of characters that come in pairs, like quotes, parentheses, brackets, etc.
a<char> – inclusive inner match: same as above, except it includes <char> on both ends
Say that I have the code below, and I want to completely replace the code inside the map() with something else:
signal.map(count -> String.format(“%d cookies, ah ah ah”, count));
The first thing I need to do is get my cursor anywhere inside the parentheses belonging to map(). The exact command I use depends on where I end up:
If the cursor is just inside the open paren for map(), I could use “cf)”. This corrects all text up to and including the next “)” on this line.
Say the cursor is on the word “format”. I would do “ci(”. Vim will search backward to find the first open paren, then search forward to find its match, and correct the text between (but not including) those characters.
Maybe my cursor was already closest to the word “cookies.” To break out of the inner parentheses, I would need to add a count and do “2ci(”. This is almost identical to the last example, except that the 2 is needed due to the nested parentheses.

######################################################################################################################################################################################################################################






########## ALPHABETICAL GIT COMMANDS (MASTER LIST IS IN BACKEDUP/Software/DamonWrittenSoftware/) ##################################################################################################################################
cat .git/HEAD                        see what the HEAD is current pointing to
git add .                            add all changed files to the staging area
git branch                           shows the names of all local branches
git branch --all                     shows the names of all local AND remote branches
git branch --help                    show help for git branch command
git branch -v                        shows the last commit on each branch
git branch --verbose                 same as above
git branch name                      create a new branch named "name" at the HEAD of the current branch
git branch name startpoint           create a new branch named "name" at startpoint of the current branch, startpoint can be a branch name, a commit SHA (id number), or a tag name, or other ...
git checkout branchname              use "git switch branchname" instead, switches the folder contents and HEAD to be the last commit of branchname
git commit -m 'better phi function'  commit all staged files to the current branch
git init                             initiate a git repository in the current directory
git pull                             first does a fetch from the origin then integrates commits
git pull origin master               first does a fetch from the origin then integrates commits
git pull --tags
git push                             pushes all commits from the current branch to the default origin
git push --tags                      push all lightweight tags to the origin
git push origin HEAD --tags          push all commits from the current branch to origin, including lightweight tags
git push origin master               pushes committs from master to origin
git remote -v                        list the URLS of remote copies of the repository that send (fetch command) or recieve (push command) updates
git remote add                       git remote add origin https://github.com/damonturney/20190925_AlCl4_EMImCl_Battery_Simulation.git
git reset HEAD                       erases unstaged changes
git stash                            stashes unstaged changes
git status                           reports if there are uncommitted changes and if the local repository is ahead of the origin
git switch branchname                switches the folder contents and HEAD to be the last commit of branchname
git switch --help                    show help for git switch command
git tag                              show all tags on the current branch
git tag tagname                      add a lightweight tag to the current HEAD
git tag tagname commit               add a lightweight tag to a specific commit, where you identify the commit by a commit SHA (id number)
git tag -d tagname                   delete a tag locally (you need to delete it on the remote too)
git push --delete origin tagname     delete a tag on origin


SEE COMPLETE HISTORY OF THE GIT REPOSITORY
git log --all --graph --decorate --oneline --date=short --pretty


WITHIN THE .gitignore FILE
*                                    * means "everyhting", thus this line means "ignore everything"
!.gitattributes                      ! means "don't",  thus this line means "don't ignore .gitattributes"
!.gitignore                          means "don't ignore .gitignore"
!*.py                                don't ignore all files that end with .py
!versions.txt                        don't ignore versions.txt
!*/                                  don't ignore all subdirectories, but inside each subdirectory all other gitignore rules are applied


SETUP GIT ON LOCAL COMPUTER AND ENABLE TRANSFERS TO GITHUB WITHOUT PASSWORDS
git config --global user.name "Damon Turney"
git config --global user.email "damonturney@gmail.com"
ssh-keygen -f ~/.ssh/id_rsa -t rsa -C "damonturney@gmail.com"   And enter a password when it prompts
copy the rsa key from ~/.ssh/id_rsa.pub with pbcopy < ~/.ssh/id_rsa.pub and paste it into the SSH key entry in GitHub webaccount "Settings --> SSH and GPG keys"


INITIATE A GIT REPOSITORY ON LOCAL COMPUTER
git init
vim .gitignore
git add .
git commit -m "first commit"
git tag v0.0
git remote add origin https://github.com/damonturney/20190925_AlCl4_EMImCl_Battery_Simulation.git


COMMIT CHANGES THEN PUSH CHANGES TO REMOTE
git add .                             add all changed files to the staging area
git commit -m 'better phi function'   commit all staged files to the current branch
git push origin HEAD --tags           push all commits from current branch to origin, including lightweight tags


VERSIONING
git tag tagname commit               add a lightweight tag to a specific commit, where you identify the commit by a commit SHA (id number)


PUSH ALL BRANCHES AND ALL TAGS TO REMOTE
git push --all
git push --tags


CREATE A BRANCH THEN SWITCH TO IT
git branch name                        create a new branch named "name" at the HEAD of the current branch
git switch name                        switches the folder contents and HEAD to be the last commit of branchname






######## Old code ##############################################################################################################################################################################################


      #@printf("loop.%5.0i   superficial_current_density:%+0.7e   electrode_voltage:%+0.7e   target_cd:%+0.7e   superficial_cd_error:%+0.7e\n", simulation_loop_iteration, superficial_current_density , ss.electrode_voltage[1] , superficial_current_density_target , superficial_current_density_error)
      #println(ss.conc_A[end-150,29:35])
      #println(ss.conc_A[end-149,29:35])
      #println(ss.conc_A[end-148,29:35])
      #println("conc_A: ", conc_A_along_surface[29:31])
      #println("conc_A_eq: ", conc_A_eq_along_surface)
      #println("V_eq: ", voltage_eq_along_surface[29:31])
      #println("V: ", ss.electrode_voltage[1])
      #println("molar_flux: ", molar_flux[29:31])
      #println("", findall(indices_above_eq))





while (abs(superficial_current_density_error) > 0.1) & (abs(ss.electrode_voltage[1]) < 0.5)
   # FOR DISCRETIZATION HELP SEE 2014 DISSERTATION BY ZANGANA !! IF there’s a flux boundary condition at x = 0 then: du/dt = (2*D*u_i+1,j - 2*D*u_i,j )/deltax^2 - 2*D*g/deltax,  u_i,j+1 = u_i,j + 2r*u_i+1,j - 2r*u_i,j - 2r*g*deltax   where g is du/dx on the boundary at x =0 and the 2 comes from the fact that we divide by 0.5deltax, not full deltax (because of boundary at x=0).  If there's a LHS "half boundary" due to a corner on the spike then:  du/dt = ( D*( u_i+1,j - u_i,j )/deltax - (0.5*D*g + 0.5*D*( u_i,j - u_i-1,j )/deltax) )/0.75/deltax,  du/dt = 4/3*D*( u_i+1,j - u_i,j )/deltax^2 - 2/3*D*( u_i,j - u_i-1,j )/deltax^2 - 2/3*D*g/deltax ,   u_i,j+1 = u_i,j + 4/3*r*( u_i+1,j - u_i,j ) - 2/3*r*( u_i,j - u_i-1,j ) - 2/3*r*g*deltax = 2/3*r*u_i-1,j + (1 - 2*r)*u_i,j + 4/3*r*u_i+1,j - 2/3*r*g*deltax       and if the boundary condition is a RHS "half boundary" then this is: u_i,j+1 = 4/3*r*u_i-1,j + (1 - 2*r)*u_i,j + 2/3*r*u_i+1,j + 2/3*r*g*deltax
   # Calculate the Fickian change in concentration due to y-direction gradients
   for i_x in 1:ss.spike_num_x_mps - 1
      conc_next_timestep[short_y_num+1, i_x]                        = 2*r_y*ss.conc_A[short_y_num      , i_x]                      + (1 - 2*r_y)*ss.conc_A[short_y_num+1, i_x]                                                                                         +   2*r_y*ss.dy*molar_flux[i_x]/ss.Diffusivity
   end
   i_x = ss.spike_num_x_mps
   conc_next_timestep[short_y_num+1, i_x]                           =4/3*r_y*ss.conc_A[short_y_num, i_x]                           + (1 - 2*r_y)*ss.conc_A[short_y_num+1, i_x]                          +2/3*r_y*ss.conc_A[short_y_num+2, i_x]                         + 2/3*r_y*ss.dy*molar_flux[ss.spike_num_x_mps]/ss.Diffusivity  #This line is for the corner mesh point.  The 0.5 factor is because this meshpoint only has electrode interface on half it's area.
   conc_next_timestep[ss.num_y_mps,  i_x]                           = 2*r_y*ss.conc_A[ss.num_y_mps - 1, i_x]                       + (1 - 2*r_y)*ss.conc_A[ss.num_y_mps, i_x]                                                                                          +   2*r_y*ss.dy*molar_flux[ss.spike_num_x_mps + ss.spike_num_y_mps+1]/ss.Diffusivity  #This line is for the corner mesh point.  The 0.5 factor is because this meshpoint only has electrode interface on half it's area.
   for i_x in ss.spike_num_x_mps+1:ss.num_x_mps                                           
      conc_next_timestep[ss.num_y_mps      , i_x]                   =2*r_y*ss.conc_A[ss.num_y_mps - 1  , i_x]                      + (1 - 2*r_y)*ss.conc_A[ss.num_y_mps, i_x]                                                                                          +   2*r_y*ss.dy*molar_flux[i_x + ss.spike_num_y_mps + 1]/ss.Diffusivity
   end
   conc_increment_dy[:,:] = conc_next_timestep[:,:] - ss.conc_A[:,:]             #Calculate how much the concentration increases due to y-direction gradients         
   # Calculate the Fickian change in concentration due to x-direction gradients
   i_y = short_y_num+1
   conc_next_timestep[i_y, ss.spike_num_x_mps]                      =2/3*r_x*ss.conc_A[i_y , ss.spike_num_x_mps-1]                 + (1 - 2*r_x)*ss.conc_A[i_y , ss.spike_num_x_mps]                    +4/3*r_x*ss.conc_A[i_y , ss.spike_num_x_mps+1]                 - 2/3*r_x*ss.dx*(-molar_flux[ss.spike_num_x_mps+1]/ss.Diffusivity)
   for i_y in short_y_num+2:ss.num_y_mps-1
      conc_next_timestep[i_y , ss.spike_num_x_mps]                  = 2*r_x*ss.conc_A[i_y , ss.spike_num_x_mps + 1]                + (1 - 2*r_x)*ss.conc_A[i_y , ss.spike_num_x_mps]                                                                                   -  2*r_x*ss.dx*(-molar_flux[ss.spike_num_x_mps + i_y - short_y_num]/ss.Diffusivity)
   end
   i_y = ss.num_y_mps
   conc_next_timestep[i_y , ss.spike_num_x_mps]                     = 2*r_x*ss.conc_A[i_y , ss.spike_num_x_mps+1]                  + (1 - 2*r_x)*ss.conc_A[i_y , ss.spike_num_x_mps]                                                                                   -  2*r_x*ss.dx*(-molar_flux[ss.spike_num_x_mps + ss.spike_num_y_mps + 1]/ss.Diffusivity)
   conc_increment_dx[:,:] = conc_next_timestep[:,:] - ss.conc_A[:,:]  #Calculate how much the concentration increases due to y-direction gradients
   conc_next_timestep[:,:] = ss.conc_A[:,:] + conc_increment_dy[:,:] + conc_increment_dx[:,:] # Add the Fickian change in concentration due to x-direction and y-direction gradients together!
   conc_increment_dy[:,:] .= 0.0
   conc_increment_dx[:,:] .= 0.0

   # Now let's see if the concentration overshoots the equilibrium concentration, or overshoots 0.0000 moles/m^3, and then let's calculate the more accurate current density
   #Now check if concentrations crossed the equilibrium concentration due to the discretized Fickian diffusion. If so, set conc_A equal to the equilibrium, and calculate current density (the real current density on the corner of the spike will be ~2x higher than the Butler-Volmer boundary condition, which was calculated 40 lines above, because the concentration of the corner meshpoint gets reduced by the x and y direction boundary conditions simultaneously)
   conc_A_along_surface_nexttimestep[1:ss.spike_num_x_mps]                                       = conc_next_timestep[short_y_num+1,1:ss.spike_num_x_mps]
   conc_A_along_surface_nexttimestep[ss.spike_num_x_mps+1:ss.spike_num_x_mps+ss.spike_num_y_mps] = conc_next_timestep[short_y_num+1:ss.num_y_mps,ss.spike_num_x_mps]
   conc_A_along_surface_nexttimestep[ss.spike_num_x_mps+ss.spike_num_y_mps+1:end]                = conc_next_timestep[ss.num_y_mps,ss.spike_num_x_mps:ss.num_x_mps]
   conc_B_along_surface_nexttimestep[:]                                                          = ss.total_conc .- conc_A_along_surface_nexttimestep[:]

   while length(findall(indices_crossed_eq)) > 0
      ##### The line below is not physics-based.  You need to come up with a way to fix this problem based on physical equations.
      conc_A_along_surface[indices_crossed_eq]                      = conc_A_along_surface[indices_crossed_eq] .+ (1 - exp(-abs.(conc_A_along_surface_nexttimestep[indices_crossed_eq] .- conc_A_along_surface[indices_crossed_eq])./abs.(conc_A_eq_along_surface .- conc_A_along_surface_nexttimestep[indices_crossed_eq])))*( conc_A_eq_along_surface .- conc_A_along_surface[indices_crossed_eq] ) 
      conc_B_along_surface[indices_crossed_eq]                      = ss.total_conc .- conc_A_along_surface[indices_crossed_eq]
      ss.conc_A[short_y_num+1,1:ss.spike_num_x_mps]                 = conc_A_along_surface[1:ss.spike_num_x_mps]
      ss.conc_A[short_y_num+1:ss.num_y_mps,ss.spike_num_x_mps]      = conc_A_along_surface[ss.spike_num_x_mps+1:ss.spike_num_x_mps+ss.spike_num_y_mps]
      ss.conc_A[ss.num_y_mps,ss.spike_num_x_mps:ss.num_x_mps]       = conc_A_along_surface[ss.spike_num_x_mps+ss.spike_num_y_mps+1:end]
      voltage_eq_along_surface[indices_crossed_eq] = V_eq.(conc_A_along_surface[indices_crossed_eq], conc_B_along_surface[indices_crossed_eq], ss.conc_A[1,50], ss.total_conc - ss.conc_A[1,50])
      overvoltage[indices_crossed_eq] = ss.electrode_voltage[1] .- voltage_eq_along_surface[indices_crossed_eq]
      current_density = -96500*ss.reaction_k .*sqrt.(conc_A_along_surface[:].*conc_B_along_surface[:]).* ( exp.(-(1.0 .- ss.Beta)*96500/8.3/300 .*overvoltage[:] ) .-  exp.(ss.Beta*300/8.3/300 .*overvoltage[:]) )  #(A/m2)
   else
      ss.electrode_voltage[1] = ss.electrode_voltage[1] - superficial_current_density_error*1E-7
      voltage_eq_along_surface[indices_crossed_eq] = V_eq.(conc_A_along_surface[indices_crossed_eq], conc_B_along_surface[indices_crossed_eq], ss.conc_A[1,50], ss.total_conc - ss.conc_A[1,50])
      overvoltage[indices_crossed_eq] = ss.electrode_voltage[1] .- voltage_eq_along_surface[indices_crossed_eq]
      current_density = -96500*ss.reaction_k .*sqrt.(conc_A_along_surface[:].*conc_B_along_surface[:]).* ( exp.(-(1.0 .- ss.Beta)*96500/8.3/300 .*overvoltage[:] ) .-  exp.(ss.Beta*300/8.3/300 .*overvoltage[:]) )  #(A/m2)
   end
   
   superficial_current_density       = mean(current_density[:])*length(current_density[:])/ss.num_x_mps
   superficial_current_density_error = superficial_current_density - superficial_current_density_target
   
   #@printf("loop:%5.0i   superficial_current_density:%+0.7e   electrode_voltage:%+0.7e   target_cd:%+0.7e   superficial_cd_error:%+0.7e\n", simulation_loop_iteration, superficial_current_density , ss.electrode_voltage[1] , superficial_current_density_target , superficial_current_density_error)
   #sleep(0.5)
end


      indices_above_eq = simdata.conc_A_along_surface .> conc_A_eq_along_surface




      # See if the Fickian diffusion caused "overshooting" wherein conc_A erroneously crosses over it's equilibrium value of conc_A_eq_along_surface
      indices_below_eq = conc_A_along_surface .< conc_A_eq_along_surface
      indices_crossed_eq = (indices_below_eq .& indices_above_eq) .| (.~indices_below_eq .& .~indices_above_eq)
      if length(findall(indices_crossed_eq)) > 0
         # The line below is not physics-based.  You need to come up with a way to fix this problem based on physical equations.
         conc_A_along_surface[indices_crossed_eq]                   = conc_A_along_surface_previous[indices_crossed_eq] .+ (1 - exp.(-abs.(conc_A_along_surface_previous[indices_crossed_eq] .- conc_A_along_surface[indices_crossed_eq])./abs.(conc_A_eq_along_surface .- conc_A_along_surface_previous[indices_crossed_eq])))*( conc_A_eq_along_surface .- conc_A_along_surface_previous[indices_crossed_eq] ) 
         conc_B_along_surface[indices_crossed_eq]                   = ss.total_conc .- conc_A_along_surface[indices_crossed_eq]
         ss.conc_A[short_y_num+1,1:ss.spike_num_x_mps]              = conc_A_along_surface[1:ss.spike_num_x_mps]
         ss.conc_A[short_y_num+1:ss.num_y_mps,ss.spike_num_x_mps]   = conc_A_along_surface[ss.spike_num_x_mps+1:ss.spike_num_x_mps+ss.spike_num_y_mps]
         ss.conc_A[ss.num_y_mps,ss.spike_num_x_mps:ss.num_x_mps]    = conc_A_along_surface[ss.spike_num_x_mps+ss.spike_num_y_mps+1:end]
         current_density[indices_crossed_eq]                        = current_density[indices_crossed_eq] .* exp.(-abs.(conc_A_along_surface_previous[indices_crossed_eq] .- conc_A_along_surface[indices_crossed_eq])./abs.(conc_A_eq_along_surface .- conc_A_along_surface_previous[indices_crossed_eq]))
      end




