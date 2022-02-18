function disp(this_lft)
%% DISP overloaded function for Ulft object
%     disp(this_lft) 
%
%     Variables:
%     ---------
%       Input:
%         this_lft : Ulft object :: the lft to be displayed
%
%     See also Ulft.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Get basic properties of LFT at first time-instant
[dim_out, dim_in] = size(this_lft);
dim_out = dim_out(1);
dim_in  = dim_in(1);

if isempty(this_lft.delta.names)
    dim_states  = 0;
    number_uncs = 0;
    dim_unc_out = 0;
    dim_unc_in  = 0;
else
    delta_state_indices = strcmp(this_lft.delta.types, 'DeltaDelayZ')...
                          | strcmp(this_lft.delta.types, 'DeltaIntegrator');
    dim_states = sum(this_lft.delta.dim_outs(delta_state_indices', 1));

    number_uncs = length(this_lft.delta.names) - sum(delta_state_indices);
    dim_unc_out = sum(this_lft.delta.dim_outs(:,1)) - dim_states;
    dim_unc_in  = sum(this_lft.delta.dim_ins(:,1)) - dim_states;
end    
%% Print LFT data
s1 = sprintf(['A [%2d, %2d]-eventually periodic LFT '...
             'with %3d outputs, %3d inputs, %3d states, and\n'...
             '%3d uncertainties forming a %3d x %3d uncertainty block\n'],...
             this_lft.horizon_period(1), this_lft.horizon_period(2),...
             dim_out, dim_in, dim_states,...
             number_uncs, dim_unc_out, dim_unc_in);
fprintf(s1);

% Print Delta, Disturbance, and Performance data
if (~isempty(this_lft.delta.names) > 0); disp(this_lft.delta); end
notDefault = @(sig) ~(length(sig.names) == 1 &&...
                      strcmp(sig.names, 'default_l2'))...
                    &&...
                    ~isempty(sig.names);
if notDefault(this_lft.disturbance);  disp(this_lft.disturbance);  end
if notDefault(this_lft.performance);  disp(this_lft.performance);  end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)