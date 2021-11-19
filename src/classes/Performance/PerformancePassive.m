classdef PerformancePassive < Performance
%% PERFORMANCEPASSIVE class for passivity measure, this extends 
% the base class Performance.
%
%   extended methods:
%     PerformancePassive(name, chan_out, chan_in, horizon_period) :: Constructor
%     disp(this_perf) :: Display method
%     matchHorizonPeriod(this_perf, total_time) :: Matches performance properties to new horizon_period
%     performanceToMultiplier(this_perf) :: Generate multiplier from performance

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

methods
    
function this_perf = PerformancePassive(name, chan_out, chan_in, horizon_period)
%% PERFORMANCEPASSIVE constructor
%
%  p = PerformancePassive(name, chan_out, chan_in, horizon_period)
%  p = PerformancePassive(name, chan_out, chan_in) assumes horizon_period == [0, 1]
%  p = PerformancePassive(name) also assumes chan_out == chan_in == {}
%
%  Variables:
%  ---------
%     Input:
%       name : char array :: unique ID of the performance spec (ex. 'pos_power')
%       chan_in : cell array of naturals :: channel of input signals pertinent to performance metric
%                                         if set to empty, chan_in is assumed to include all input signals of
%                                         any Ulft this performance is attached to
%       chan_out : cell array of naturals :: channel of output signals pertinent to performance metric
%                                         if set to empty, chan_out is assumed to include all output signals of
%                                         any Ulft this performance is attached to
%       horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of Performance properties  
%
%  See also PerformancePassive, Performance.Performance

    % Defining defaults for missing arguments
    switch nargin
        case 0
            error('PerformancePassive:PerformancePassive',...
                  'Must provide at least the name of the performance')
        case 1
            chan_out       = {[]};
            chan_in        = {[]};
            horizon_period = [0, 1];
        case 3
            horizon_period = [0, 1];
        case 4
        otherwise
            error('PerformancePassive:PerformancePassive',...
                ['Must provide 1, 3, or 4 arguments to construct',...
                 'PerformancePassive objects'])
    end
    dim_out = cellfun(@length, chan_out);
    dim_in  = cellfun(@length, chan_in);
    assert(all(dim_out == dim_in),...
           'PerformancePassive:PerformancePassive',...
           'Passive Performance must have equal lengths of channels')
    
    % Calling Performance constructor
    this_perf@Performance(name, chan_out, chan_in, horizon_period);
    this_perf = matchHorizonPeriod(this_perf);
end

function disp(this_perf)
%% DISP function for PerformancePassive object
%
%  disp(perf_l2_obj) (e.g., disp(PerformancePassive('p')) )
%
%  Variables:
%  ---------
%     Input:
%       this_perf : PerformancePassive object     
%
%  See also Ulft.disp, SequencePerformance.disp, Performance.disp

    disp@Performance(this_perf, 'passivity')
end

function this_perf = matchHorizonPeriod(this_perf, new_horizon_period)
%% MATCHHORIZONPERIOD function to ensure properties of PerformancePassive
%  object match its own horizon_period, or a new_horizon_period
%
%  this_perf = matchHorizonPeriod(this_perf, new_horizon_period) will change the horizon_period and associated properties of this_perf
%  this_perf = matchHorizonPeriod(this_perf) only ensures that each pertinent properties of this_perf matches this_perf.horizon_period
%
%  Variables:
%  ---------
%     Input:
%       this_perf : PerformancePassive object
%       new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence
%     Output:
%        this_perf : PerformancePassive object
%
%  See also PerformancePassive.

if nargin == 1
% Ensuring that this_perf.horizon_period matches with other properties 
% of this_perf

    % Assumes properties are a mix of sequences of length 1 or
    % horizon_period
    total_time = sum(this_perf.horizon_period);
    if length(this_perf.chan_out) ~= total_time
        assert(length(this_perf.chan_out) == 1,...
               'PerformancePassive:matchHorizonPeriod',...
               'output channels of %s is not compatible w/ horizon_period',...
               this_perf.name);
        this_perf.chan_out = repmat(this_perf.chan_out, 1, total_time);
    end
    if length(this_perf.chan_in) ~= total_time
        assert(length(this_perf.chan_in) == 1,...
               'PerformancePassive:matchHorizonPeriod',...
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
%       this_perf : PerformancePassive object
%     Output:
%       multiplier : MultiplierPassive object
%
%  See also PerformancePassive

input_parser = inputParser;
addRequired(input_parser,...
            'performance',...
            @(perf) validateattributes(perf,...
                                       'PerformancePassive',...
                                       {'nonempty'}));
addParameter(input_parser,...
            'dim_out_lft',...
            [],...
            @(dim) validateattributes(dim,...
                                      'numeric',...
                                      {'integer', 'row', 'positive'}));
addParameter(input_parser,...
            'dim_in_lft',...
            [],...
            @(dim) validateattributes(dim,...
                                      'numeric',...
                                      {'integer', 'row', 'positive'}));
addParameter(input_parser,... % This parameter is not used. Only defined for compatibility with other performanceToMultiplier calls
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
if isempty(dim_out_lft)
    error('PerformancePassive:performanceToMultiplier',...
          ['Input arguments to performanceToMultiplier must include',...
           ' a dim_out_lft as a key/value pair. For example: \n',...
           'performanceToMultiplier(perf, ''dim_out_lft'', 1,...)'])
end
if isempty(dim_in_lft)
    error('PerformancePassive:performanceToMultiplier',...
          ['Input arguments to performanceToMultiplier must include',...
           ' a dim_in_lft as a key/value pair. For example: \n',...
           'performanceToMultiplier(perf, ''dim_in_lft'', 1,...)'])
end
multiplier = MultiplierPerformancePassive(performance, dim_out_lft, dim_in_lft);
end
end
end

%%  CHANGELOG
% Oct. 7, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)