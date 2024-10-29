"Migration model..." is the code for the simulation and creating figures. Contains all necessary functions and packages. The simulation program "sim" is on line 1208.

The other "stats" files are saved simulations at equilibrium - can be used to create a variety of plots from existing runs. But not for using a saved population - those files are too big to upload, so need to run yourself. (3-5 minutes)

## Workflow ##

1) Run all the preamble and simulation function (I usually just run everything until line 1635). Set your working directory.
2) Run a large simulation (~5000 years x 5 reps takes 3-4 minutes) and save the population (L 1638). See line 1138-1207 for info on sim function.
3) Run statstmp and save (L 1639).
If you want to plot, run popplots, either from local files (e.g L 1640, uses arguments "local=statstmp,mainfile=test") or saved (e.g. L 1647, uses argument "saved=Control stats.R"). Vectors 'shorter' and 'allbutpopsize' (L 839) give examples of what you can plot.

For running a climate change scenario (from line 1652), you need to first have a metapopulation at equilibrium, either locally or saved. (point 2-3 above)
Then: 
1) Set up your climate change scenario (using arguments shock.time and shock.size in sim) and run.
2) Run functions statstmp, create.alldb, polish.alldb, and create.evoldyn on your output (L1692-1702, takes 1-2 minutes). Save files if you want. Note that default climate change scenarios 1 and 2 last for 200 years, whereas scenarios 3-5 last for only 100 years. In create.evoldyn you need to specify (using argument 'years') which years you want to output - use years=26:200 for CC1 and CC2, years=26:100 for CC3-5.
3) Take a look at databases alldb and evol if you want (L 1695, 1700).
4) Now the code for plotting climate change simulations should work, either using function evolplots (L 1784-1794, but remember to change all '175' to '75' depending on which climate change scenario you are doing!) or example code on L 1705-1781.

Some more code for things to do with alldb data is provided from line 1800 and out.
Good luck and come find me if you have questions! (Room D1-215, just across the hall from Irja's office)
