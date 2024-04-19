pm_v4_air.m
Calculates bending angle profile of one occultation
Called in dop2alpha_pm_v4.m

read_ar_v1_test.m	
Reads excess phase file
Called in dop2alpha_pm_v4.m

dop2alpha_pm_v4.m
Takes excess phase file and calculates and plots the bending angle profile of one occultation
Called in loop_dop2alpha.m

loop_dop2alpha.m
Takes directory of excess phase files and calculates the bending angle profiles of one IOP

graph_from_outputs_v2.m
Plots comparative bending angle profiles from PM, GO, and forward model on one plot for one occultation
Called in graph_wrapper.m

graph_wrapper.m
Plots comparative bending angle profiles from PM, GO, and forward model on one plot for one IOP

binned_std.m
Calculates the mean, standard deviation, and edges for y based on bins of x
Called in plot_all_v1.m

plot_all_v1.m
Plots OmB or OmO for two of the three: PM, GO, forward model
For one IOP

To get from excess phase files to OmB for entire IOP:
loop_dop2alpha.m: excess phase files -> bending angle files (PM)
plot_all_v1.m: 2 types of bending angle files -> OmB plots





binned_std.m				x needs to be commented					
plot_all_v1.m				x
pm_v4_air.m				x
dop2alpha_pm_v4.m			x
read_ar_v1_test.m			x needs to be commented
graph_from_outputs_v2.m			x	
graph_wrapper.m				x
loop_dop2alpha.m			x needs to be commented
transform.m