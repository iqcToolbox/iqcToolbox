classdef DisturbanceConstantWindow < Disturbance
%% DISTURBANCECONSTANTWINDOW class for the set of L2 signals that are constant
% in a window of time-steps. This extends the base class Disturbance.
%
%   extended methods:
%     DisturbanceConstantWindow(name, channel, window, horizon_period) :: Constructor
%     disp(this_dis) :: Display method
%     matchHorizonPeriod(this_perf, horizon_period) :: Matches disturbance properties to new horizon_period
%     disturbanceToMultiplier(this_dis, 'dim_in_lft', dim_in_lft) :: Generate multiplier from disturbance
%
%   extended properties:
%     window : array of naturals :: time-instances in which the disturbance signal is constant from the time-instances immediately 
%                                   preceding the given time-instances(index starting from t0 = 0)
%                 this window is assumed to repeat periodically according to the horizon_period.
%                 For example, if the signals are constant for all time, (i.e., d_1 = d_0, d_2 = d_1, etc.)
%                    and horizon_period = [0, 1], set window = 1;
%                 Instead, if horizon_period = [0, 2], (i.e., it is [0, 2]-eventually periodic)
%                    and d_0 = d_1, d_2 = d_3, d_4 = d_5, etc. (the signal is constant within each period, but changes between periods)
%                    set window = 1 (this specifies d_1 - d_0 = 0, and the constraint is extended every two steps because the system is [0, 2]-eventually periodic
%                 Instead, if horizon_period = [3, 3],
%                    the signal during the first three time instances (i.e. in the non-periodic portion) is constant
%                    but then may change arbitrarily afterward,
%                    set window = 1:2 (this specifices d_1 = d_0, and d_2 = d_1) 
%                 Instead, if horizon_period = [3, 3], 
%                   the signal during the first three time instances (i.e. in the non-periodic portion) is constant
%                    and d_3 = d_2, d_6 = d_5, d_9 = d_8, d_12 = d_11, etc. (i.e., the first two timesteps during every period are constant)
%                    then set window = 1:3 (this specifies d_1 = d_0, d_2 = d_1, and d_3 = d_2. Then, because d_3 pertains to the periodic portion the d_3i = d_3i - 1 constraint repeats every period)
%     override : boolean :: indicates that user wishes to override an error thrown when
%                 window connects timesteps between the non-periodic and periodic segments    
%
%  See also DisturbanceConstantWindow.DisturbanceConstantWindow

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    window (1, :) double {mustBeInteger, mustBePositive}
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
    %          window : array of naturals :: time-instances in which the disturbance signal is constant between the given time-instances and the immediately 
    %                                            preceding time-instances (index starting from t0 = 0)
    %                                            this window is assumed to repeat periodically according to the horizon_period.
    %                                            For example, if the signals are constant for all time, (i.e., d_1 = d_0, d_2 = d_1, etc.)
    %                                                and horizon_period = [0, 1], set window = 1;
    %                                            Instead, if horizon_period = [0, 2], (i.e., it is [0, 2]-eventually periodic)
    %                                                and d_0 = d_1, d_2 = d_3, d_4 = d_5, etc. (the signal is constant within each period, but changes between periods)
    %                                                set window = 1 (this specifies d_1 - d_0 = 0, and the constraint is extended every two steps because the system is [0, 2]-eventually periodic
    %                                            Instead, if horizon_period = [3, 3],
    %                                                the signal during the first three time instances (i.e. in the non-periodic portion) is constant
    %                                                but then may change arbitrarily afterward,
    %                                                set window = 1:2 (this specifices d_1 = d_0, and d_2 = d_1) 
    %                                            Instead, if horizon_period = [3, 3], 
    %                                               the signal during the first three time instances (i.e. in the non-periodic portion) is constant
    %                                                and d_3 = d_2, d_6 = d_5, d_9 = d_8, d_12 = d_11, etc. (i.e., the first two timesteps during every period are constant)
    %                                                then set window = 1:3 (this specifies d_1 = d_0, d_2 = d_1, and d_3 = d_2. Then, because d_3 pertains to the periodic portion the d_3i = d_3i - 1 constraint repeats every period)
    %                                        Specifying that the last timestep of the non-periodic portion (i.e., ismember(horizon_period(1), window))
    %                                            is the same as the first timestep of the periodic portion has 
    %                                            the consequence of also specifying the last timestep of the
    %                                            periodic portion is the same as the first timestep in the next period.
    %                                            Because of this, errors are thrown unless the user acknowledges such behavior
    %                                            by setting `override = true`.
    %                                            For example, if hp = [1, 2], and window = [1, 2], this specifies that
    %                                            d_1 == d_0, d_2 == d_1, d_4 == d_3, etc.
    %                                            If you wish to only specify u_0 == u_1, expand the non-periodic portion of the 
    %                                            horizon_period by one timestep and specify the same window:
    %                                            from the last example, this would mean changing to hp = [2, 2], and keeping window = [1, 2].
    %                                            original_hp = lft.horizon_period;
    %                                            new_hp = original_hp + [1, 0];
    %                                            lft = lft.matchHorizonPeriod(new_hp)
    %                                            dis = DisturbanceConstantWindow('dis', {[]}, [1, 2], new_hp;
    %                                            lft = lft.addDisturbance({dis})
    %          horizon_period : 1 x 2 array of naturals :: horizon and period of properties of disturbance class
    %          override : boolean :: indicates that the user wishes to override an error thrown by connecting timesteps
    %                                            from the non-periodic portion to the periodic portion
    %       Output:
    %          this_dis : DisturbanceConstantWindow object :: the produced disturbance object specifying the admissible set of disturbances
    %
    %     See also DisturbanceConstantWindow
    
    % Defining defaults for missing arguments
    switch nargin
        case 0
            error('DisturbanceConstantWindow:DisturbanceConstantWindow',...
                  'Must provide at least the name of the disturbance')
        case 1
            % This produces the trivial "constant at all time" signal (which must therefore be 0)
            % override set to true to enable this, but override should generally be set to false
            chan_in = {[]};
            window = 1;
            horizon_period = [0, 1];
            override = true;
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
    withinHorizonPeriod = @(w) max(w) < total_time;
    fillsHorizonPeriod = @(w) isequal(w, 1:total_time);
    % Only allow ismember(total_time, window) if treating corner-case "0" signal 
    assert(withinHorizonPeriod(window) || fillsHorizonPeriod(window),...
           'DisturbanceConstantWindow:DisturbanceConstantWindow',...
           ['The provided window is not consistent with the horizon_period.',...
            ' Ensure all indices in the window would occur',...
            ' within the initial horizon and first period specified by',...
            ' horizon_period. Recall that time indices start at t0 = 0'])
    % Error logic if window spans non-periodic and periodic portion
    betweenNonperiodAndPeriod = @(w, hp) ismember(hp(1), w);
    if betweenNonperiodAndPeriod(window, this_dis.horizon_period) && ~override
        error('DisturbanceConstantWindow:DisturbanceConstantWindow',...
           ['The provided window should not specify the first timestep',...
            ' of the periodic portion',...
            ' because the user may unintentially be specifying additional',...
            ' undesired constraints.\n\n',...
            'If you want to specify the signal is constant between the last',...
            ' non-periodic time instant and the first periodic time instant,',...
            ' (i.e., ismember(hp(1), window)), but ensure that the last instant of',...
            ' a period may be different from the first instant of the next period,',...
            ' make an LFT and disturbance whose',...
            ' hp is hp_new = hp + [1, 0], and reconstruct this',...
            ' disturbance with hp_new and the same window.\n\n',...
            'This error may be overridden by setting override = true when',...
            ' constructing this disturbance'])
    end
    uniqueTimeIndices = @(w) length(w) == length(unique(w));
    assert(uniqueTimeIndices(window),...
       'DisturbanceConstantWindow:DisturbanceConstantWindow',...
       ['The provided window has duplicate time indices, ensure the',...
        'window has only unique indices'])
    this_dis.window   = window;
    this_dis.override = override;
    
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
        if length(this_dis.window) > 2
            fprintf(['%16s with a total of %3d timesteps wherein the',...
                     ' signal is constant. \n'],...
                    '',...
                    length(this_dis.window))
        else
            
            window_str = sprintf('[%3d, %3d] ',...
                                 [this_dis.window-1; this_dis.window]);
            fprintf(['%16s which is constant between time-instances ',...
                    window_str(1 : end - 2), '] \n'],...
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
        window = [];                                                            %#ok<PROPLC>
        for i = 1 : length(this_dis.window)
            % Find when new time indices coincide with old window indices
            window = [window, find(this_dis.window(i) == (indices - 1))];       %#ok<AGROW,PROPLC>
        end
        window = sort(window - 1);                                              %#ok<PROPLC>
        if isequal(this_dis.window, 1:sum(this_dis.horizon_period))
        % Corner case where signal is constant through all time
            window = 1:sum(new_horizon_period);                                 %#ok<PROPLC>
        end
        % Shift window indices to start at t0 = 0
        this_dis.window         = window;                                       %#ok<PROPLC>
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
                                          {'DisturbanceConstantWindow'},...
                                          {'nonempty'},...
                                          mfilename));
    addParameter(input_parser,...
                'dim_in_lft',...
                [],...
                @(dim) validateattributes(dim,...
                                          {'numeric'},...
                                          {'integer', 'row', 'positive'},...
                                          mfilename));
    addParameter(input_parser,...
                 'quad_time_varying',...
                 true,...
                 @(quad) validateattributes(quad, {'logical'}, {'nonempty'}))
    addParameter(input_parser,... % This parameter is not used. Only defined for compatibility with other disturbanceToMultiplier calls
                 'discrete',...
                 true,...
                 @(disc) validateattributes(disc, {'logical'}, {'nonempty'}))

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