%function graph_wrapper(path)
%
% Plots bending angle vs impact height for results from phase matching,
% geometric optics, and forward model on one graph, for every occultation
% in a directory
%  
%% Inputs: 
%   path: a string that is the path to the directory
%
%% Instructions: 
%   User must change the hard-coded section to the appropriate directories
%   for phase matching, geometric optics, and forward model data
%     
%% Dependencies
%  graph_from_outputs_v2.m
%
%% Example:
%
% graph_wrapper('output')
%

function graph_wrapper(path)

%hard-coded
file_start = 'Output_alpha';
occ_index_start = 18;
occ_index_end = 24;
norp = 1;

%rest of file
len = length(file_start);
directory = dir(path);

for file = directory'
    
    if length(file.name) > 2

        if file.name(1:len) == file_start

            occ_code = file.name(occ_index_start:occ_index_end);
            graph_from_outputs(occ_code, norp)
    
        end

    end

end