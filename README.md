"Migration model..." is the code for the simulation and creating figures. Contains all necessary functions and packages. The simulation program "sim" is on line 1208.

The other "stats" files are saved simulations at equilibrium - can be used to create a variety of plots from existing runs. But not for using a saved population - those files are too big to upload, so need to run yourself. (3-5 minutes)

## Workflow ##

1: Run all the preamble and simulation function (I usually just run everything until line 1630.). Set your working directory.
2: Run a large simulation (~5000 years x 5 reps takes 3-4 minutes) and save the population (L 1638-39).
3: Run statstmp and save (L 1640-41)
If you want to plot, run popplots, either from local files (use arguments "local=statstmp,mainfile=test") or saved (use argument "saved=Control stats.R")

For running a climate change scenario (from line 1650), you need to first have a metapopulation at equilibrium, either locally or saved. (point 2-3 above)
Then: 
1) Set up your climate change scenario (using arguments shock.time and shock.size in sim) and run.
2) Run functions create.alldb, polish.alldb, and create.evoldyn on your output (takes 1-2 minutes)
