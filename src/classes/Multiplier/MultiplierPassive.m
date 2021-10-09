classdef MultiplierPassive < MultiplierDelta
%% MUTLIPLIERPASSIVE class. Used for the passive performance metric (PerformancePassive)
%  Extends the base class MultiplierPerformance
%
%  extended methods:
%    MultiplierPassive(performance) :: Constructor
%
%  See also MultiplierPassive.MultiplierPassive

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

methods
function this_mult = MultiplierPassive(delta)
%% MULTIPLIERPASSIVE constructor
%
%  this_mult = MultiplierPassive(delta)
%
%  Variables:
%  ---------
%    Input:
%      delta : DeltaPassive object :: delta defining this multiplier
%                                           
%    Output:
%      this_mult : MultiplierPassive object

% Inputs
validateattributes(delta, 'DeltaPassive', {'nonempty'}, mfilename)


this_mult.name           = delta.name;
this_mult.horizon_period = delta.horizon_period;

% Define filter
total_time = sum(delta.horizon_period);
filter.a     = cell(1, total_time);
filter.b1    = cell(1, total_time);
filter.b2    = cell(1, total_time);
filter.c1    = cell(1, total_time);
filter.c2    = cell(1, total_time);
filter.d11   = cell(1, total_time);
filter.d12   = cell(1, total_time);
filter.d21   = cell(1, total_time);
filter.d22   = cell(1, total_time);
quad.q11     = cell(1, total_time);
quad.q12     = cell(1, total_time);
quad.q21     = cell(1, total_time);
quad.q22     = cell(1, total_time);
for i = 1:total_time
    % Define filter matrices
    filter.a{i}   = [];
    filter.b1{i}  = zeros(0, delta.dim_out(i));
    filter.b2{i}  = zeros(0, delta.dim_in(i));
    filter.c1{i}  = zeros(delta.dim_out(i), 0);
    filter.c2{i}  = zeros(delta.dim_in(i), 0);
    filter.d11{i} = eye(delta.dim_out(i));
    filter.d12{i} = zeros(delta.dim_out(i), delta.dim_in(i));
    filter.d21{i} = zeros(delta.dim_in(i), delta.dim_out(i));
    filter.d22{i} = eye(delta.dim_in(i));

    % Define quad
    quad.q12{i} = eye(delta.dim_out(i));
    quad.q21{i} = quad.q12{i}';
    quad.q11{i} = zeros(delta.dim_out(i), delta.dim_in(i));
    quad.q22{i} = zeros(delta.dim_in(i), delta.dim_out(i));    
end       
this_mult.filter = filter;
this_mult.quad   = quad;
end    
end

end

%%  CHANGELOG
% Oct. 7, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)