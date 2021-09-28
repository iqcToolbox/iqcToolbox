classdef MultiplierL2 < MultiplierDisturbance
%% MULTIPLIERL2 class for L2 disturbance signals (DisturbanceL2).
%  Extends the base class MultiplierDisturbance.
%
%  extended methods:
%    MultiplierDisturbance(disturbance, dim_in_lft) :: Constructor
%
%  extended properties:
%    chan_in : cell array of column vectors of naturals :: The input channels of the LFT pertaining 
%                                                          to this performance
%    dim_in : row of naturals :: The input dimension of the LFT pertaining to this performance
%
%  See also MultiplierL2.MultiplierL2

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    chan_in (1, :) cell
    dim_in (1, :) double
end
    
methods
function this_mult = MultiplierL2(disturbance, dim_in_lft)
%% MULTIPLIERL2 constructor
%
%  this_mult = MultiplierL2(disturbance, dim_in_lft)
%
%  Variables:
%  ---------
%    Input:
%      disturbance : DisturbanceL2 object :: disturbance defining this multiplier
%      dim_in_lft : row of naturals :: input dimension of LFT associated with disturbance
%    Output:
%      this_mult : MultiplierL2 object
%
%  See also MultiplierL2
input_parser = inputParser;
addRequired(input_parser,...
            'disturbance',...
            @(dis) validateattributes(dis,...
                                      'DisturbanceL2',...
                                      {'nonempty'},...
                                      mfilename));
addRequired(input_parser,...
            'dim_in_lft',...
            @(dim) validateattributes(dim,...
                                      'numeric',...
                                      {'integer', 'row', 'positive'},...
                                      mfilename));
                              
parse(input_parser, disturbance, dim_in_lft)
disturbance = input_parser.Results.disturbance;
dim_in_lft  = input_parser.Results.dim_in_lft;

consistent_dims = cellfun(@(chan, dim) maxEmpty(chan) <= dim,...
                          disturbance.chan_in,...
                          num2cell(dim_in_lft));
assert(all(consistent_dims),...
       'MultiplierL2:MultiplierL2',...
       'Disturbance channels exceeds the input dimensions of lft')

this_mult.name           = disturbance.name;
this_mult.horizon_period = disturbance.horizon_period;
this_mult.chan_in        = disturbance.chan_in;
this_mult.dim_in         = dim_in_lft;

% Define filter and quad
total_time = sum(this_mult.horizon_period);
filter.a = cell(1, total_time);
filter.b = cell(1, total_time);
filter.c = cell(1, total_time);
filter.d = cell(1, total_time);
quad.q   = cell(1, total_time);
for i = 1:total_time
    filter.a{i} = zeros(0);
    filter.b{i} = zeros(0, this_mult.dim_in(i));
    filter.c{i} = zeros(this_mult.dim_in(i), 0);
    filter.d{i} = eye(this_mult.dim_in(i));
    quad.q{i} = zeros(this_mult.dim_in(i));
end
this_mult.filter = filter;
this_mult.quad = quad;
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)