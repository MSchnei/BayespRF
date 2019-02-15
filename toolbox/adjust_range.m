function x_adj = adjust_range(x,d_min_old,d_max_old, d_min_new, d_max_new)
% Adjust latent variable for new range of neural parameters
%
% x         - latent variable whose range is to be adjusted
% d_min_old - old minimum value of neural parameters
% d_max_old - old minimum value of neural parameters
% d_min_new - new minimum value of neural parameters
% d_max_new - new minimum value of neural parameters
%
% x_adj - latent variable with adjusted range

% Convert latent variable to neural parameter in old space
d_range = d_max_old - d_min_old;
x2      = (d_range .* normcdf(x,0,1)) + d_min_old;

% Convert to new range
d2_range = d_max_new - d_min_new;
x3 = (x2 - d_min_new) / d2_range;

% Convert to latent variable
x_adj = norminv(x3, 0, 1);
