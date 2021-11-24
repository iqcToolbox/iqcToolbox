classdef MultiplierTimeWindow < MultiplierDisturbance
%% MULTIPLIERTIMEWINDOW class for time-windowed disturbance signals (DisturbanceTimeWindow).
%  Extends the base class MultiplierDisturbance.
%
%  extended methods:
%    MultiplierTimeWindow(disturbance, varargin) :: Constructor
%
%  extended properties:
%    chan_in : cell array of column vectors of naturals :: The input channels of the LFT pertaining 
%                                                          to this disturbance
%    window : array of integers :: the timesteps (indexing from 0) in which the disturbance
%                                  is non-zero
%
%  See also MultiplierTimeWindow.MultiplierTimeWindow

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    chan_in (1, :) cell
    dim_in (1, :) double
    window (1, :) double {mustBeInteger, mustBeNonnegative}
end

properties (SetAccess = immutable)
    quad_time_varying logical
end
    
methods
function this_mult = MultiplierTimeWindow(disturbance, dim_in_lft, varargin)
%% MULTIPLIERTIMEWINDOW constructor
%
%  multiplier = disturbanceToMultiplier(disturbance, dim_in_lft, 'quad_time_varying', true)
%  multiplier = disturbanceToMultiplier(disturbance, dim_in_lft) assumes the input above
%
%  Variables:
%  ---------
%     Input:
%       disturbance : DisturbanceTimeWindow object
%       dim_in_lft : row of naturals :: input dimension of LFT associated with disturbance
%       quad_time_varying : logical :: Whether or not the decision variable
%                            quad is time-varying (if true, SDP will possibly
%                            be less conservative, but more computationally
%                            challenging)
%     Output:
%       multiplier : MultiplierTimeWindow object
%
%  See also DisturbanceTimeWindow
input_parser = inputParser;
addRequired(input_parser,...
            'disturbance',...
            @(dis) validateattributes(dis,...
                                      'DisturbanceTimeWindow',...
                                      {'nonempty'},...
                                      mfilename));
addRequired(input_parser,...
            'dim_in_lft',...
            @(dim) validateattributes(dim,...
                                      'numeric',...
                                      {'integer', 'row', 'positive'},...
                                      mfilename));
addParameter(input_parser,...
             'quad_time_varying',...
             true,...
             @(quad) validateattributes(quad, 'logical', {'nonempty'}))

parse(input_parser, disturbance, dim_in_lft, varargin{:})
disturbance       = input_parser.Results.disturbance;
quad_time_varying = input_parser.Results.quad_time_varying;

consistent_dims = cellfun(@(chan, dim) maxEmpty(chan) <= dim,...
                          disturbance.chan_in,...
                          num2cell(dim_in_lft));
assert(all(consistent_dims),...
       'MultiplierTimeWindow:MultiplierTimeWindow',...
       'Disturbance channels exceeds the input dimensions of lft')
   
% Define multiplier
this_mult.name              = disturbance.name;
this_mult.horizon_period    = disturbance.horizon_period;
this_mult.chan_in           = disturbance.chan_in;
this_mult.dim_in            = dim_in_lft;
this_mult.window            = disturbance.window;
this_mult.quad_time_varying = quad_time_varying;

% Define filter and quad
total_time    = sum(this_mult.horizon_period);
filter.a      = cell(1, total_time);
filter.b      = cell(1, total_time);
filter.c      = cell(1, total_time);
filter.d      = cell(1, total_time);
quad.q        = cell(1, total_time);
ct = [];
if this_mult.quad_time_varying
    decision_vars = cell(1, total_time - length(this_mult.window));
else
    decision_vars = cell(1);
end

ind_dec_var = 1;
for i = 1:total_time
    % Filter ss matrices
    filter.a{i} = zeros(0);
    filter.b{i} = zeros(0, this_mult.dim_in(i));
    filter.c{i} = zeros(this_mult.dim_in(i), 0);
    if ismember(i, this_mult.window + 1)
        filter.d{i} = zeros(this_mult.dim_in(i));
        quad.q{i}   = zeros(this_mult.dim_in(i));
    else
        % Create decision variables
        first_idx_outside_window = ...
            find(~ismember(1:total_time, this_mult.window + 1), 1, 'first');
                                                  
        if (i == first_idx_outside_window) || ... 
           this_mult.quad_time_varying
            p = sdpvar(1);
            ct = ct + ((p >= 0):['Time Window Multiplier, ',...
                                 this_mult.name,...
                                 ', p >= 0, #',...
                                 num2str(i)]);                                  %#ok<BDSCA>
            decision_vars{ind_dec_var} = p;
            ind_dec_var = ind_dec_var + 1;
        end
        filter.d{i} = eye(this_mult.dim_in(i));
        if ~isempty(this_mult.chan_in{i})
            quad.q{i} = -p * double(diag(ismember(1:this_mult.dim_in(i),...
                                                  this_mult.chan_in{i}')));
        else
            quad.q{i} = -p * eye(this_mult.dim_in(i));
        end
    end
end
this_mult.filter        = filter;
this_mult.quad          = quad;
this_mult.decision_vars = decision_vars;
this_mult.constraints   = ct;
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)