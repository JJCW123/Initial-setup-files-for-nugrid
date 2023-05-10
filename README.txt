Attached in the GitHub are the initial setup files that were altered. Every other file that was used was left with default settings and therefore will not be left here.

To activate the flux output files that are described in the output section, a file named 'ppn\_frame.input' contains a parameter called 'iplot_flux_option =' which would need to be set 
to value 1, with its corresponding option below of 'i_flux_integrated ='. These two options allow for integrated flux at each a chosen time-step determined by the final parameter 
'print_cycles =' which sets the interval at which flux files are created. This file also contains parameters 'ini_filename =', 'trajectry_fn =' and 'nsource =' these in turn control
the initial abundances file name which is in the above paragraph, the trajectory input file name which was always labelled 'trajectory.dat' and lastly, an option to input whether a source
file would be provided or it should run at a constant temperature and density chosen by the user.


To alter specific reaction rates the physics package had to be opened and a file named 'physics_knobs.f90' accessed to change the two commands 'rate_index(:) =' and 'rate_factor(:) ='. 
The first command identifies the specific reaction rate that you are targeting given by an integer from the 'networksetup.txt' file, then the factor command can be changed to vary the 
reaction rate for the simulation. The network setup text had the 15O(α, γ)19Ne at index number : 398 and the 22Mg(α, p)25Al rate index was: 1321.
