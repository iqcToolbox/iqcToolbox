classdef DisturbanceConstantWindow < Disturbance
%% DISTURBANCECONSTANTWINDOW class for the set of L2 signals that are non-zero
% in a window of time-steps. This extends the base class Disturbance.
%
%   extended methods:
%     DisturbanceConstantWindow(name, channel, window, horizon_period) :: Constructor
%     disp(this_dis) :: Display method
%     matchHorizonPeriod(this_perf, horizon_period) :: Matches disturbance properties to new horizon_period
%     disturbanceToMultiplier(this_dis) :: Generate multiplier from disturbance
%
%   extended properties:
%     window : array of naturals :: set of contiguous time-steps in which signal is constant (index starting from t0 = 0)
%                 this window is assumed to repeat periodically according to the horizon_period.
%                 For example, if the signals are constant for all time, 
%                     and horizon_period = [0, 1], set window = 0;
%                     if horizon_period = [2, 3], set window = [0:4]. 
%                 Instead, if horizon_period = [3, 2], (i.e., it is [3, 2]-eventually periodc) 
%                     the initial timestep of a signal is nonzero, 
%                     and then switches on-and-off during the periodic portion,
%                     set window = [0, 3]. 
%                 Instead, if horizon_period = [3, 2], (i.e., it is [3, 2]-eventually periodic) 
%                     the first three timesteps (i.e., in the non-periodic portion) are constant
%                     then set window = [0:2];
%                 Instead, if horizon_period = [3, 3],
%                     and the first two timesteps of each period are constant,
%                     then set window = [3:4];
%     override : boolean :: indicates that user wishes to override an error thrown when
%                 window connects timesteps between the non-periodic and periodic segments    
%
%  See also DisturbanceConstantWindow.DisturbanceConstantWindow

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    window (1, :) double {mustBeInteger, mustBeNonnegative}
    override logical
end

methods
    function this_dis = DisturbanceConstantWindow(name,...
                                                  chan_in,...
                                                  window,...
                                                  horizon_period,...
                                                  override)
    %% DISTURBANCECONSTANTWINDOW constructor
    %
    %     dis = DisturbanceConstantWindow(name, channel, window, horizon_period, override)
    %     dis = DisturbanceConstantWindow(name, channel, window, horizon_period) assumes override == false
    %     dis = DisturbanceConstantWindow(name) also assumes channel == {[]}, window == [1], and horizon_period == [0, 1].  This results in expressing the entire set of L2 signals
    %     
    %     Variables:
    %     ---------
    %       Input:
    %          name : char array :: unique ID of disturbance signal (ex. 'sensor_noise')
    %          channel : cell array of columns of naturals :: channels of signals pertaining to disturbance class
    %          window : array of naturals :: set of contiguous time-steps wherein the signal is constant
    %                                        (index starting from t0 = 0)
    %                                        this window is assumed to repeat periodically according to the horizon_period.
    %                                        For example, if the signals are constant for all time, 
    %                                            and horizon_period = [0, 2], set window = [0:1];
    %                                            if horizon_period = [2, 3], set window = [0:4]. 
    %                                        Instead, if horizon_period = [3, 2], (i.e., it is [3, 2]-eventually periodc) 
    %                                            the first three timesteps (i.e., in the non-periodic portion) 
    %                                            then set window = [0:2];
    %                                        Instead, if horizon_period = [3, 3],
    %                                            and the first two timesteps of each period are constant,
    %                                            then set window = [3:4];
    %                                        Specifying that the last timestep of the non-periodic portion
    %                                            is the same as the first timestep of the periodic portion has the
    %                                            the consequence of also specifying the last timestep of the
    %                                            periodic portion is the same as the first timestep in the next period.
    %                                            Because of this, errors are thrown unless the user acknowledges such behavior
    %                                            by setting `override = true`.
    %                                            For example, if hp = [1, 2], and window = [0, 1], this specifies that
    %                                            u_0 == u_1, u_2 == u_3, u_4 == u_5, etc.
    %                                            If you wish to only specify u_0 == u_1, use a one-step-long horizon in the
    %                                            horizon_period: in the last example, hp = [2, 2], and window = [0, 1].
    %                                            original_hp = lft.horizon_period;
    %                                            new_hp = original_hp + [1, 0];
    %                                            lft = lft.matchHorizonPeriod(new_hp)
    %                                            dis = DisturbanceConstantWindow('dis', {[]}, [0, 1], new_hp;
    %                                            lft = lft.addDisturbance({dis})
    %          horizon_period : 1 x 2 array of naturals :: horizon and period of properties of disturbance class
    %          override : boolean :: indicates that the user wishes to override an error thrown by connecting timesteps
    %                                            from the non-periodic portion to the periodic portion
    %       Output:
    %          this_dis : DisturbanceConstantWindow object :: the produced disturbance set
    %
    %     See also DisturbanceConstantWindow
    
    % Defining defaults for missing arguments
    switch nargin
        case 0
            error('DisturbanceConstantWindow:DisturbanceConstantWindow',...
                  'Must provide at least the name of the disturbance')
        case 1
            chan_in = {[]};
            window = 0;
            horizon_period = [0, 1];
            override = false;
        case 4
            override = false;
        case 5
        otherwise
            error('DisturbanceConstantWindow:DisturbanceConstantWindow',...
                  ['Must provide 1, 4, or 5 arguments to construct',...
                   'DisturbanceConstantWindow objects'])
    end
    
    % Checking conditions on chan_in
    assert(length(chan_in) == 1,...
           'DisturbanceConstantWindow:DisturbanceConstantWindow',...
           ['chan_in must be constant for each timestep, ',...
            'provide chan_in as a cell of length 1'])
    % Calling Disturbance constructor
    this_dis@Disturbance(name, chan_in, horizon_period);
    
    % Checking inputs for specialized properties of DisturbanceConstantWindow
    total_time = sum(this_dis.horizon_period);
    windowWithinHorizonPeriod = @(w) (max(w) + 1) <= total_time;
    assert(windowWithinHorizonPeriod(window),...
           'DisturbanceConstantWindow:DisturbanceConstantWindow',...
           ['The provided window is not consistent with the horizon_period.',...
            ' Ensure all indices in the window would occur',...
            ' within the initial horizon and first period specified by',...
            ' horizon_period. Recall that time indices start at t0 = 0'])
    contiguousTimeIndices = @(w) isequal(w, w(1):w(end));
    assert(contiguousTimeIndices(window),...
           'DisturbanceConstantWindow:DisturbanceConstantWindow',...
           'The provided window must have continuous timesteps')
    assert(length(window) > 1,...
           'DisturbanceConstantWindow:DisturbanceConstantWindow',...
           'The provided window must specify at least two timesteps')
    % Error logic if window spans non-periodic and periodic portion
    horizonOnly = @(w, hp) max(w) < hp(1);
    periodOnly = @(w, hp) min(w) >= hp(1);
    if ~(horizonOnly(window, horizon_period) ||...
         periodOnly(window, horizon_period) || ...
         override)
        error('DisturbanceConstantWindow:DisturbanceConstantWindow',...
           ['The provided window should not specify timesteps that span',...
            ' the non-periodic and periodic portion of a horizon_period',...
            ' because the user may unintentially be specifying additional',...
            ' undesired constraints.\n\n',...
            'If you want to specify constant signals only for the first X',...
            ' timesteps (and X > hp(1)), make an LFT and disturbance whose',...
            ' hp is hp_new = hp + [X - hp(1), 0], and reconstruct this',...
            ' disturbance with hp_new and the same window.\n\n',...
            'This error may be overridden by setting override = true when',...
            ' constructing this disturbance'])
    end
    this_dis.window = window;
    
    this_dis = matchHorizonPeriod(this_dis);
    end

    function disp(this_dis)
    %% DISP function for DisturbanceConstantWindow object
    %
    %  disp(dis_cw_obj) (e.g., disp(DisturbanceConstantWindow('dis')) )
    %
    %  Variables:
    %  ---------
    %    Input:
    %       this_dis : DisturbanceConstantWindow object     
    %
    %  See also Ulft.disp, SequenceDisturbance.dsip, Disturbance.disp
    
        disp@Disturbance(this_dis, 'Constant-in-a-window ')
        if length(this_dis.window) > 5
            fprintf(['%16s with a total of %3d time-instances where the',...
                     ' signal is constant. \n'],...
                    '',...
                    length(this_dis.window))
        else
            window = sprintf( '%2d, ', this_dis.window);
            fprintf(['%16s which is constant between time-instances ',...
                     '[', window(1 : end - 2), '] \n'],...
                    '')
        end  
    end

    function this_dis = matchHorizonPeriod(this_dis, new_horizon_period)
    %% MATCHHORIZONPERIOD function to ensure properties of DisturbanceConstantWindow
    %  object match its own horizon_period, or a new_horizon_period
    %
    %  matchHorizonPeriod(this_dis, new_horizon_period) will change the horizon_period and associated properties of this_dis
    %  matchHorizonPeriod(this_dis) only ensures that each pertinent property of this_dis matches this_dis.horizon_period
    %
    %  Variables:
    %  ---------
    %     Input:
    %       this_dis : DisturbanceConstantWindow object
    %       new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence
    %     Output:
    %       this_dis : DisturbanceConstantWindow object
    %
    %  See also DisturbanceConstantWindow.
    
    if nargin == 1
    % Ensuring that this_dis.horizon_period matches with other properties 
    % of this_dis
    
        % Assumes properties are a mix of sequences of length 1 or
        % horizon_period
        total_time = sum(this_dis.horizon_period);
        if length(this_dis.chan_in) ~= total_time
            assert(length(this_dis.chan_in) == 1,...
                   'DisturbanceConstantWindow:matchHorizonPeriod',...
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
    %       this_dis : DisturbanceConstantWindow object
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
        error('DisturbanceConstantWindow:disturbanceToMultiplier',...
              ['Input arguments to disturbanceToMultiplier must include',...
               ' a dim_in_lft as a key/value pair. For example: \n',...
               'disturbanceToMultiplier(dis, ''dim_in_lft'', 1,...)'])
    end
    multiplier = MultiplierConstantWindow(disturbance, dim_in_lft);
    end
end
end

%%  CHANGELOG
% Nov. 18, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)