classdef MultiplierConstantWindow < MultiplierDisturbance
%% MULTIPLIERCONSTANT class for disturbance signals which are constant in given
%  windows of time (DisturbanceConstantWindow).
%  Extends the base class MultiplierDisturbance.
%
%  extended methods:
%    MultiplierConstantWindow(disturbance, varargin) :: Constructor
%
%  extended properties:
%    chan_in : cell array of column vectors of naturals :: The input channels of the LFT pertaining 
%                                                          to this disturbance
%    window : array of integers :: the set of multiple, contiguous timesteps (indexing from 0) in which 
%                                  the disturbance is constant
%
%  See also MultiplierConstantWindow.MultiplierConstantWindow

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
function this_mult = MultiplierConstantWindow(disturbance, dim_in_lft, varargin)
%% MULTIPLIERCONSTANTWINDOW constructor
%
%  multiplier = disturbanceToMultiplier(disturbance, 'dim_in_lft', dim_in_lft, 'quad_time_varying', true)
%  multiplier = disturbanceToMultiplier(disturbance, 'dim_in_lft', dim_in_lft) assumes the input above
%
%  Variables:
%  ---------
%     Input:
%       disturbance : DisturbanceConstantWindow object
%       dim_in_lft : row of naturals :: input dimension of LFT associated with disturbance
%       quad_time_varying : logical :: Whether or not the decision variable
%                            quad is time-varying (if true, SDP will possibly
%                            be less conservative, but more computationally
%                            challenging)
%     Output:
%       multiplier : MultiplierConstantWindow object
%
%  See also DisturbanceConstantWindow
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
       'MultiplierConstantWindow:MultiplierConstantWindow',...
       'Disturbance channels exceeds the input dimensions of lft')
   
% Define multiplier
this_mult.name              = disturbance.name;
this_mult.horizon_period    = disturbance.horizon_period;
this_mult.chan_in           = disturbance.chan_in;
this_mult.dim_in            = dim_in_lft;
this_mult.window            = disturbance.window;
this_mult.quad_time_varying = quad_time_varying;
this_mult.discrete          = true;

% Define filter
total_time = sum(this_mult.horizon_period);
chan_in    = this_mult.chan_in{1};
if isempty(chan_in)
    dim_state = dim_in_lft(1);
    chan_select = eye(dim_state);
else
    dim_state  = length(chan_in);
    cols = [eye(dim_state), zeros(dim_state, 1)];
    chan_select = [cols(:, chan_in), zeros(dim_state, dim_in_lft - dim_state)];
end
chan_select = repmat({chan_select}, 1, total_time);
chan_select = toLft(chan_select, this_mult.horizon_period);
delay      = DeltaDelayZ(dim_state, -1, this_mult.horizon_period);
filter_lft = delay - eye(dim_state);
filter_lft = filter_lft * chan_select;
filter.a = filter_lft.a;
filter.b = filter_lft.b;
filter.c = filter_lft.c;
filter.d = filter_lft.d;

% % % Produces (for each timestep) a vector indicating which 
% % % columns of the cols matrix should appear
% % channels = cell(1, total_time);
% % for i = 1:total_time
% %     chan_now = this_mult.chan_in{i};
% %     if isempty(chan_now)
% %         % Select all channels
% %         channels{i} = [1:dim_state]';
% %     else
% %         % Select and sort channels depending on chan_in
% %         channels{i} = cols(:, chan_now);
% %         % Add zero columns for remaining 
% %         channels{i}(end+1:dim_state) = dim_state + 1;
% %         
% % chan_select = cellfun(@(chan) cols(:, chan'), this_mult.chan_in,...
% %                       'UniformOutput', false);
% % chan_select = toLft(chan_select, this_mult.horizon_period);
% % filter_lft = filter_lft * chan_select;
% % filter.a = filter_lft.a;
% % filter.b = filter_lft.b;
% % filter.c = filter_lft.c;
% % filter.d = filter_lft.d;

% Define quad
quad.q        = cell(1, total_time);
ct = [];
if this_mult.quad_time_varying
    decision_vars = cell(1, length(this_mult.window));
else
    decision_vars = cell(1);
end

ind_dec_var = 1;
for i = 1:total_time
    if ~ismember(i, this_mult.window + 1) || i == this_mult.window(1) + 1
    % If not part of window, or just the first timestep, set q to zero
        quad.q{i}   = zeros(dim_state);
    else
        % Create decision variables
        if (i == this_mult.window(1) + 2) || this_mult.quad_time_varying
            p = sdpvar(1);
            ct = ct + ((p >= 0):['Constant Window Multiplier, ',...
                                 this_mult.name,...
                                 ', p >= 0, #',...
                                 num2str(i)]);                                  %#ok<BDSCA>
            decision_vars{ind_dec_var} = p;
            ind_dec_var = ind_dec_var + 1;
        end
        quad.q{i} = -p * eye(dim_state);
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
% Nov. 18, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)