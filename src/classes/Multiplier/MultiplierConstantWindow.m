classdef MultiplierConstantWindow < MultiplierDisturbance
%% MULTIPLIERCONSTANTWINDOW class for disturbance signals which are constant in given
%  windows of time (DisturbanceConstantWindow).
%  Extends the base class MultiplierDisturbance.
%
%  extended methods:
%    MultiplierConstantWindow(disturbance, varargin) :: Constructor
%
%  extended properties:
%    chan_in : cell array of column vectors of naturals :: The input channels of the LFT pertaining 
%                                                          to this disturbance
%    window : array of integers :: time-instances in which the disturbance signal is constant between 
%                                  the given time-instances and the immediately preceding time-instances 
%                                  (index starting from t0 = 0)
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
%  multiplier = MultiplierConstantWindow(disturbance, dim_in_lft, 'quad_time_varying', true)
%  multiplier = MultiplierConstantWindow(disturbance, dim_in_lft) assumes the input above
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
                                      'DisturbanceConstantWindow',...
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
dim_in_lft        = input_parser.Results.dim_in_lft;
quad_time_varying = input_parser.Results.quad_time_varying;

assert(sum(disturbance.horizon_period) == length(dim_in_lft),...
       'MultiplierConstantWindow:MultiplierConstantWindow',...
       'dim_in_lft must have the same length as total_time of the disturbance')
consistent_dims = cellfun(@(chan, dim) maxEmpty(chan) <= dim,...
                          disturbance.chan_in,...
                          num2cell(dim_in_lft));
assert(all(consistent_dims),...
       'MultiplierConstantWindow:MultiplierConstantWindow',...
       'Disturbance channels exceeds the input dimensions of lft')
assert(all(dim_in_lft(1) == dim_in_lft),...
       'MultiplierConstantWindow:MultiplierConstantWindow',...
       ['LFT input dimensions must be constant to use',...
        'DisturbanceConstantWindow in IQC analysis'])   % This can be relaxed in the future

   
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
    dim_state = this_mult.dim_in(1);
    chan_select = eye(dim_state);
else
    dim_state  = length(chan_in);
    cols = [eye(dim_state), zeros(dim_state, 1)];
    chan_select = [cols(:, chan_in), zeros(dim_state,...
                   this_mult.dim_in(1) - dim_state)];
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

% Define quad
quad.q = cell(1, total_time);
ct     = [];
if this_mult.quad_time_varying
    decision_vars = cell(1, length(this_mult.window));
else
    decision_vars = cell(1);
end

ind_dec_var = 1;
for i = 1:total_time
    if ismember(i, this_mult.window + 1)
    % If i is within (window + 1), set q(i) = -p I
        % Create decision variables
        if ind_dec_var == 1 || this_mult.quad_time_varying
            p = sdpvar(1);
            ct = ct + ((p >= 0):['Constant Window Multiplier, ',...
                                 this_mult.name,...
                                 ', p >= 0, #',...
                                 num2str(i)]);                                  %#ok<BDSCA>
            decision_vars{ind_dec_var} = p;
            ind_dec_var = ind_dec_var + 1;
        end
        quad.q{i} = -p * eye(dim_state);
    else
        quad.q{i} = zeros(dim_state);
    end
end
if isequal(this_mult.window, 1:total_time)
% A corner-case where the entire signal is specified constant (and therefore 0)
    p = sdpvar(1);
    ct = ct + ((p >= 0):['Constant Window Multiplier, ',...
                         this_mult.name,...
                         ', p >= 0, #',...
                         num2str(1)]);                                          %#ok<BDSCA>
    decision_vars{ind_dec_var} = p;
    quad.q{1} = -p * eye(dim_state);
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