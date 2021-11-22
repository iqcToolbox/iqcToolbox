classdef DisturbanceTimeWindow < Disturbance
%% DISTURBANCETIMEWINDOW class for the set of L2 signals that are non-zero
% in a window of time-steps. This extends the base class Disturbance.
%
%   extended methods:
%     DisturbanceTimeWindow(name, channel, window, horizon_period) :: Constructor
%     disp(this_dis) :: Display method
%     matchHorizonPeriod(this_perf, horizon_period) :: Matches disturbance properties to new horizon_period
%     disturbanceToMultiplier(this_dis, 'dim_in_lft', dim_in_lft) :: Generate multiplier from disturbance
%
%   extended properties:
%     window : array of naturals :: time-steps in which signal is non-zero (index starting from t0 = 0)
%                 this window is assumed to repeat periodically according to the horizon_period.
%                 For example, if the signals are non-zero for all time, 
%                     and horizon_period = [0, 1], set window = 0;
%                     if horizon_period = [2, 3], set window = [0:4]. 
%                 Instead, if horizon_period = [3, 2], (i.e., it is [3, 2]-eventually periodc) 
%                     the initial timestep of a signal is nonzero, 
%                     and then switches on-and-off during the periodic portion,
%                     set window = [0, 3]. 
%
%  See also DisturbanceTimeWindow.DisturbanceTimeWindow

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    window (1, :) double {mustBeInteger, mustBeNonnegative}
end

methods
    function this_dis = DisturbanceTimeWindow(name,...
                                              chan_in,...
                                              window,...
                                              horizon_period)
    %% DISTURBANCETIMEWINDOW constructor
    %
    %     dis = DisturbanceTimeWindow(name, channel, window, horizon_period)
    %     dis = DisturbanceTimeWindow(name) assumes channel == {[]}, window == [1], and horizon_period == [0, 1].  This results in expressing the entire set of L2 signals
    %     
    %     Variables:
    %     ---------
    %       Input:
    %          name : char array :: unique ID of disturbance signal (ex. 'sensor_noise')
    %          channel : cell array of naturals :: channels of signals pertaining to disturbance class
    %          window : matrix array of naturals :: time-steps wherein signal is non-zero (index starting from t0 = 0)
    %                                        this window is assumed to repeat periodically according to the horizon_period.
    %                                        For example, if the signals are non-zero for all time, 
    %                                            and horizon_period = [0, 1], set window = 0;
    %                                            if horizon_period = [2, 3], set window = [0:4]. 
    %                                        Instead, if horizon_period = [3, 2], (i.e., it is [3, 2]-eventually periodc) 
    %                                            the initial timestep of a signal is nonzero, 
    %                                            and then switches on-and-off during the periodic portion,
    %                                            set window = [0, 3]. 
    %          horizon_period : 1 x 2 array of naturals :: horizon and period of properties of disturbance class
    %       Output:
    %          this_dis : DisturbanceTimeWindow object :: the produced disturbance object specifying the admissible set of disturbances
    %
    %     See also DisturbanceTimeWindow
    
    % Defining defaults for missing arguments
    switch nargin
        case 0
            error('DisturbanceTimeWindow:DisturbanceTimeWindow',...
                  'Must provide at least the name of the disturbance')
        case 1
            chan_in = {[]};
            window = 0;
            horizon_period = [0, 1];
        case 4
        otherwise
            error('DisturbanceTimeWindow:DisturbanceTimeWindow',...
                  ['Must provide 1 or 4 arguments to construct',...
                   'DisturbanceTimeWindow objects'])
    end    
    % Calling Disturbance constructor
    this_dis@Disturbance(name, chan_in, horizon_period);
    
    % Checking inputs for specialized properties of DisturbanceTimeWindow
    total_time = sum(this_dis.horizon_period);
    windowWithinHorizonPeriod = @(w) (max(w) + 1) <= total_time;
    assert(windowWithinHorizonPeriod(window),...
           'DisturbanceTimeWindow:DisturbanceTimeWindow',...
           ['The provided window is not consistent with the horizon_period',...
            '. Ensure all indices in the window would occur',...
            ' within the initial horizon and first period specified by',...
            ' horizon_period. Recall that time indices start at t0 = 0'])
    uniqueTimeIndices = @(w) length(w) == length(unique(w));
    assert(uniqueTimeIndices(window),...
           'DisturbanceTimeWindow:DisturbanceTimeWindow',...
           ['The provided window has duplicate time indices, ensure the',...
            'window has only unique indices'])
    this_dis.window = window;
    
    this_dis = matchHorizonPeriod(this_dis);
    end

    function disp(this_dis)
    %% DISP function for DisturbanceTimeWindow object
    %
    %  disp(dis_tw_obj) (e.g., disp(DisturbanceTimeWindow('dis')) )
    %
    %  Variables:
    %  ---------
    %    Input:
    %       this_dis : DisturbanceTimeWindow object     
    %
    %  See also Ulft.disp, SequenceDisturbance.dsip, Disturbance.disp
    
        disp@Disturbance(this_dis, 'Time-windowed')
        if length(this_dis.window) > 5
            fprintf(['%16s with a total of %3d time-instances where the',...
                     ' signal is non-zero. \n'],...
                    '',...
                    length(this_dis.window))
        else
            window = sprintf( '%3d, ', this_dis.window);
            fprintf(['%16s which is non-zero during time-instances ',...
                     '[', window(1 : end - 2), '] \n'],...
                    '')
        end  
    end

    function this_dis = matchHorizonPeriod(this_dis, new_horizon_period)
    %% MATCHHORIZONPERIOD function to ensure properties of DisturbanceTimeWindow
    %  object match its own horizon_period, or a new_horizon_period
    %
    %  matchHorizonPeriod(this_dis, new_horizon_period) will change the horizon_period and associated properties of this_dis
    %  matchHorizonPeriod(this_dis) only ensures that each pertinent property of this_dis matches this_dis.horizon_period
    %
    %  Variables:
    %  ---------
    %     Input:
    %       this_dis : DisturbanceTimeWindow object
    %       new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence
    %     Output:
    %       this_dis : DisturbanceTimeWindow object
    %
    %  See also DisturbanceTimeWindow.
    
    if nargin == 1
    % Ensuring that this_dis.horizon_period matches with other properties 
    % of this_dis
    
        % Assumes properties are a mix of sequences of length 1 or
        % horizon_period
        total_time = sum(this_dis.horizon_period);
        if length(this_dis.chan_in) ~= total_time
            assert(length(this_dis.chan_in) == 1,...
                   'DisturbanceTimeWindow:matchHorizonPeriod',...
                   'channels of %s is not compatible w/ horizon_period',...
                   this_dis.name);
            this_dis.chan_in = repmat(this_dis.chan_in, 1, total_time);
        end
    else
    % Changing this_dis.horizon_period and other properties of this_dis to
    % a new horizon_period
    
        [indices, new_horizon_period] = ...
            makeNewIndices(this_dis.horizon_period, new_horizon_period);

        % Set properties according to new indices
        this_dis.chan_in        = this_dis.chan_in(indices);
        window = [];
        for i = 1 : length(this_dis.window)
            % Find when new time indices coincide with old window indices
            window = [window, find(this_dis.window(i) == (indices - 1))];
        end
        window = sort(window - 1);
        % Shift window indices to start at t0 = 0
        this_dis.window         = window;
        this_dis.horizon_period = new_horizon_period;
        this_dis                = matchHorizonPeriod(this_dis);
    end
    end
    
    function multiplier = disturbanceToMultiplier(this_dis, varargin)
    %% DISTURBANCETOMULTIPLIER function to generate a multiplier from this object. 
    %
    %  multiplier = disturbanceToMultiplier(this_dis, 'dim_in_lft', dim_in_lft, 'quad_time_varying', true)
    %  multiplier = disturbanceToMultiplier(this_dis, 'dim_in_lft', dim_in_lft) assumes the input above
    %
    %  Variables:
    %  ---------
    %     Input:
    %       this_dis : DisturbanceTimeWindow object
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
    addParameter(input_parser,...
                'dim_in_lft',...
                [],...
                @(dim) validateattributes(dim,...
                                          'numeric',...
                                          {'integer', 'row', 'positive'},...
                                          mfilename));
    addParameter(input_parser,...
                 'quad_time_varying',...
                 true,...
                 @(quad) validateattributes(quad, 'logical', {'nonempty'}))
    addParameter(input_parser,... % This parameter is not used. Only defined for compatibility with other disturbanceToMultiplier calls
                 'discrete',...
                 true,...
                 @(disc) validateattributes(disc, 'logical', {'nonempty'}))

    parse(input_parser, this_dis, varargin{:})
    disturbance = input_parser.Results.disturbance;
    dim_in_lft  = input_parser.Results.dim_in_lft;
    if isempty(dim_in_lft)
        error('DisturbanceTimeWindow:disturbanceToMultiplier',...
              ['Input arguments to disturbanceToMultiplier must include',...
               ' a dim_in_lft as a key/value pair. For example: \n',...
               'disturbanceToMultiplier(dis, ''dim_in_lft'', 1,...)'])
    end
    multiplier = MultiplierTimeWindow(disturbance, dim_in_lft);
    end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)