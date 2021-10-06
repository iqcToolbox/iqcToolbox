classdef PerformanceL2Induced < Performance
%% PERFORMANCEL2INDUCED class for l2-norm performance measure, this extends 
% the base class Performance.
%
%   extended methods:
%     PerformanceL2Induced(name, chan_out, chan_in, gain, horizon_period) :: Constructor
%     disp(this_perf) :: Display method
%     matchHorizonPeriod(this_perf, total_time) :: Matches performance properties to new horizon_period
%     performanceToMultiplier(this_perf) :: Generate multiplier from performance
%
%   extended properties:
%     gain : double or empty :: l2-induced gain of lft for the given input/output channels
%                               if set to empty, the gain to verify against is not predetermined
%
%   See also PerformanceL2Induced.PerformanceL2Induced, Performance
%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    gain double {mustBePositive}
end

methods
    
function this_perf = PerformanceL2Induced(name, chan_out, chan_in, gain, horizon_period)
%% PERFORMANCEL2INDUCED constructor
%
%  p = PerformanceL2Induced(name, chan_out, chan_in, gain, horizon_period)
%  p = PerformanceL2Induced(name, chan_out, chan_in, gain) assumes horizon_period == [0, 1]
%  p = PerformanceL2Induced(name, chan_out, chan_in) also assumes gain == []
%  p = PerformanceL2Induced(name) also assumes chan_out == chan_in == {}
%
%  Variables:
%  ---------
%     Input:
%       name : char array :: unique ID of the performance spec (ex. 'pos_power')
%       chan_in : cell array of naturals :: channel of input signals pertinent to performance metric
%                                             if set to empty, chan_in is assumed to include all input signals of
%                                             any Ulft this performance is attached to
%       chan_out : cell array of naturals :: channel of output signals pertinent to performance metric
%                                             if set to empty, chan_out is assumed to include all output signals of
%                                             any Ulft this performance is attached to
%       gain : double or empty :: l2-induced gain of lft for the given input/output channels
%                                   if set to empty, the gain to verify against is not predetermined
%                                   and will be solved for in iqcAnalysis
%       horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of Performance properties  
%
%  See also PerformanceL2Induced, Performance.Performance

    % Defining defaults for missing arguments
    switch nargin
        case 0
            error('PerformanceL2Induced:PerformanceL2Induced',...
                  'Must provide at least the name of the performance')
        case 1
            chan_out       = {[]};
            chan_in        = {[]};
            gain           = [];  
            horizon_period = [0, 1];
        case 3
            gain           = [];
            horizon_period = [0, 1];
        case 4
            horizon_period = [0, 1];
        case 5
        otherwise
            error('PerformanceL2Induced:PerformanceL2Induced',...
                ['Must provide 1, 3, 4, or 5 arguments to construct',...
                 'PerformanceL2Induced objects'])
    end
    % Calling Performance constructor
    this_perf@Performance(name, chan_out, chan_in, horizon_period);
    
    % Checking inputs for specialized properties of PerformanceL2Induced
    if ~isempty(gain)
        validateattributes(gain, 'numeric', {'positive', 'real', 'row'})
        assert(all(gain(1) == gain),...
               'PerformanceL2Induced:PerformanceL2Induced',...
               'gain must be constant for PerformanceL2Induced objects')
    end
    this_perf.gain = gain;

    this_perf = matchHorizonPeriod(this_perf);
end

function disp(this_perf)
%% DISP function for PerformanceL2Induced object
%
%  disp(perf_l2_obj) (e.g., disp(PerformanceL2Induced('p')) )
%
%  Variables:
%  ---------
%     Input:
%       this_perf : PerformanceL2Induced object     
%
%  See also Ulft.disp, SequencePerformance.disp, Performance.disp

    disp@Performance(this_perf, 'L2-induced')
    if isempty(this_perf.gain)
        fprintf(' %16s with an unspecified gain \n', '')
    else
        fprintf(' %16s with a gain of %3.1e \n', '', this_perf.gain(1))
    end    
end

function this_perf = matchHorizonPeriod(this_perf, new_horizon_period)
%% MATCHHORIZONPERIOD function to ensure properties of PerformanceL2Induced
%  object match its own horizon_period, or a new_horizon_period
%
%  this_perf = matchHorizonPeriod(this_perf, new_horizon_period) will change the horizon_period and associated properties of this_perf
%  this_perf = matchHorizonPeriod(this_perf) only ensures that each pertinent properties of this_perf matches this_perf.horizon_period
%
%  Variables:
%  ---------
%     Input:
%        this_perf : PerformanceL2Induced object
%        new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence
%     Output:
%        this_perf : PerformanceL2Induced object
%
%  See also PerformanceL2Induced.

if nargin == 1
% Ensuring that this_perf.horizon_period matches with other properties 
% of this_perf

    % Assumes properties are a mix of sequences of length 1 or
    % horizon_period
    total_time = sum(this_perf.horizon_period);
    if length(this_perf.chan_out) ~= total_time
        assert(length(this_perf.chan_out) == 1,...
               'PerformanceL2Induced:matchHorizonPeriod',...
               'output channels of %s is not compatible w/ horizon_period',...
               this_perf.name);
        this_perf.chan_out = repmat(this_perf.chan_out, 1, total_time);
    end
    if length(this_perf.chan_in) ~= total_time
        assert(length(this_perf.chan_in) == 1,...
               'PerformanceL2Induced:matchHorizonPeriod',...
               'input channels of %s is not compatible w/ horizon_period',...
               this_perf.name);
        this_perf.chan_in = repmat(this_perf.chan_in, 1, total_time);
    end
else
% Changing this_perf.horizon_period and other properties of this_perf to
% a new horizon_period

    [indices, new_horizon_period] = ...
        makeNewIndices(this_perf.horizon_period, new_horizon_period);

    % Set properties according to indices
    this_perf.chan_out        = this_perf.chan_out(indices);
    this_perf.chan_in         = this_perf.chan_in(indices);
    this_perf.horizon_period  = new_horizon_period;

    
    this_perf = matchHorizonPeriod(this_perf);
end
end

function multiplier = performanceToMultiplier(performance, varargin)
%% PERFORMANCETOMULTIPLIER function to generate a multiplier from this object. 
%
%  multiplier = performanceToMultiplier(this_perf)
%
%  Variables:
%  ---------
%     Input:
%        this_perf : PerformanceL2Induced object
%     Output:
%        multiplier : MultiplierL2Induced object
%
%  See also PerformanceL2Induced

input_parser = inputParser;
addRequired(input_parser,...
            'performance',...
            @(perf) validateattributes(perf,...
                                       'PerformanceL2Induced',...
                                       {'nonempty'},...
                                       mfilename));
addParameter(input_parser,...
            'dim_out_lft',...
            [],...
            @(dim) validateattributes(dim,...
                                      'numeric',...
                                      {'integer', 'row', 'positive'},...
                                      mfilename));
addParameter(input_parser,...
            'dim_in_lft',...
            [],...
            @(dim) validateattributes(dim,...
                                      'numeric',...
                                      {'integer', 'row', 'positive'},...
                                      mfilename));
addParameter(input_parser,...
             'objective_scaling',...
             1,...
             @(os) validateattributes(os,...
                                      'numeric',...
                                      {'positive', 'real', 'nonempty'}))
addParameter(input_parser,... % This parameter is not used. Only defined for compatibility with other performanceToMultiplier calls
             'discrete',...
             true,...
             @(disc) validateattributes(disc, 'logical', {'nonempty'}))
            
parse(input_parser, performance, varargin{:})
performance = input_parser.Results.performance;
dim_out_lft = input_parser.Results.dim_out_lft;
dim_in_lft  = input_parser.Results.dim_in_lft;
objective_scaling = input_parser.Results.objective_scaling;
if isempty(dim_out_lft)
    error('PerformanceL2Induced:performanceToMultiplier',...
          ['Input arguments to performanceToMultiplier must include',...
           ' a dim_out_lft as a key/value pair. For example: \n',...
           'performanceToMultiplier(perf, ''dim_out_lft'', 1,...)'])
end
if isempty(dim_in_lft)
    error('PerformanceL2Induced:performanceToMultiplier',...
          ['Input arguments to performanceToMultiplier must include',...
           ' a dim_in_lft as a key/value pair. For example: \n',...
           'performanceToMultiplier(perf, ''dim_in_lft'', 1,...)'])
end
multiplier = MultiplierL2Induced(performance, dim_out_lft, dim_in_lft,...
                                 'objective_scaling', objective_scaling);
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)