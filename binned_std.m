% x is the independent variable - bins are based off of x
% y is the dependent variable - bins for x are applied to y

function [means, stds, edges] = binned_std(x,y,num_bins)

[bins, edges] = discretize(x, num_bins);
means = splitapply(@nanmean,y,bins);
stds = splitapply(@nanstd,y,bins);