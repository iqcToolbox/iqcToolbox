classdef (Abstract) Performance
%% PERFORMANCE abstract base class for performance measurements
% This class should not be directly accessed. Subclasses of Performance
% (such as PerformanceL2) extend and concretize this class.
%
%   concrete methods:
%     Performance(name, chan_out, chan_in, horizon_period) :: Constructor method
%     disp(this_perf) :: Display method
%
%   abstract methods:
%     matchHorizonPeriod(this_perf, horizon_period) :: Matches performance properties to a new horizon_period
%     performanceToMultiplier(this_dis, varargin) :: Method for constructing a multiplier from a performance
%
%   properties:
%     name : char array :: unique ID of the performance spec (ex. 'pos_power')
%     chan_in : cell array of naturals :: channel of input signals pertinent to performance metric
%     chan_out : cell array of naturals :: channel of output signals pertinent to performance metric
%
%   See also Performance.Performance

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    name (1, :) char
    chan_out (1, :) cell
    chan_in (1, :) cell
    horizon_period (1, 2) double {mustBeInteger, mustBeNonnegative}
end

methods (Abstract)        
    this_perf = matchHorizonPeriod(this_perf, total_time)
    multiplier = performanceToMultiplier(this_perf)
end

methods
    function this_perf = Performance(name, chan_out, chan_in, horizon_period)
    %% PERFORMANCE constructor, must be called by subclass constructor
    %
    %  this_perf = Performance(name, chan_out, chan_in, horizon_period)
    %
    %  Variables:
    %  ---------
    %    Input:
    %       name : char array :: unique ID of the performance spec (ex. 'pos_power')
    %       chan_in : cell array of naturals :: channel of input signals pertinent to performance metric
    %       chan_out : cell array of naturals :: channel of output signals pertinent to performance metric
    %       horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of Performance properties    
    %
    %  See also Performance.
    
        % Checking validity of inputs
        validateattributes(name, 'char', {'row'}, mfilename);
        if iscell(chan_out) && isempty(chan_out)
            chan_out = {[]};
        end
        validateattributes(chan_out, 'cell', {'row'}, mfilename);
        for i = 1:length(chan_out)
            validateattributes(chan_out{i},...
                                   'numeric',...
                                   {'integer', 'positive'},...
                                   mfilename);
            assert(size(chan_out{i}, 2) < 2,...
                   'Performance:Performance',...
                   ['Specified output channels must be cell array of empty',...
                    ' or column arrays'])
            assert(length(chan_out{i}) == length(unique(chan_out{i})),...
                   'Performance:Performance',...
                   'Specified output channels contains duplicates')
        end
        if iscell(chan_in) && isempty(chan_in)
            chan_in = {[]};
        end
        validateattributes(chan_in, 'cell', {'row'}, mfilename);
        for i = 1:length(chan_in)
            validateattributes(chan_in{i},...
                                   'numeric',...
                                   {'integer', 'positive'},...
                                   mfilename);
            assert(size(chan_in{i}, 2) < 2,...
                   'Performance:Performance',...
                   ['Specified input channels must be cell array of empty',...
                    ' or column arrays'])
            assert(length(chan_in{i}) == length(unique(chan_in{i})),...
                   'Performance:Performance',...
                   'Specified input channels contains duplicates')
        end  
        validateattributes(horizon_period, 'numeric', {'size', [1,2],...
                                                       'integer',...
                                                       'nonnegative'});
        validateattributes(horizon_period(2), 'numeric', {'positive'})                        
                       
        % Setting properties of Performance
        this_perf.name           = name;
        this_perf.chan_out       = chan_out;
        this_perf.chan_in        = chan_in;
        this_perf.horizon_period = horizon_period;
    end

    function disp(this_perf, type)
    %% DISP function for Performance object, must be called from a subclass disp method.
    %
    %  disp(this_perf, type)
    %
    %  Variables:
    %  ---------
    %    Input:
    %       this_dis : Performance object     
    %       type : string :: type of performance ('L2', etc.)
    
        fprintf('%4s %-7s is a %8s performance metric for ',...
                '',...
                this_perf.name,...
                type)
        if isempty(this_perf.chan_out{1})
            fprintf('\n %16s all output channels', '')
        else
            channel = sprintf( '%d; ', this_perf.chan_out{1});
            fprintf('\n %16s output channels [', '')
            fprintf([channel(1 : end - 2), ']'])
        end
        if isempty(this_perf.chan_in{1})
            fprintf(' and all input channels \n')
        else
            channel = sprintf( '%d; ', this_perf.chan_in{1});
            fprintf([' and input channels [', channel(1 : end - 2), '] \n'])
        end
    end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)