classdef DisturbanceBandedWhite < Disturbance
%% DISTURBANCEBANDEDWHITE class for the subset of L2 signals that are white
% within the given frequency band [-omega, omega]. This extends the base class 
% Disturbance.
%
%   extended methods:
%     DisturbanceBandedWhite(name, chan_in, omega, horizon_period) :: Constructor
%     disp(this_dis) :: Display method
%     matchHorizonPeriod(this_perf, horizon_period) :: Matches disturbance properties to new horizon_period
%     disturbanceToMultiplier(this_dis, 'dim_in_lft', dim_in_lft) :: Generate multiplier from disturbance
%
%  See also DisturbanceBandedWhite.DisturbanceBandedWhite, Disturbance

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    omega (1, 1) double {mustBeNonnegative}
end

methods
    function this_dis = DisturbanceBandedWhite(name, chan_in, omega, horizon_period)
    %% DISTURBANCEBandedWhite constructor
    %
    %  dis = DisturbanceBandedWhite(name, chan_in, omega, horizon_period)
    %  dis = DisturbanceBandedWhite(name, chan_in, omega) assumes horizon_period == [0, 1]
    %  dis = DisturbanceBandedWhite(name, chan_in) also assumes omega == pi
    %  dis = DisturbanceBandedWhite(name) also assumes chan_in == {[1]}
    %
    %  Variables:
    %  ---------
    %     Input:
    %       name : char array :: unique ID of disturbance signal (ex. 'sensor_noise')
    %       chan_in : cell array of natural :: channel of signals pertaining to disturbance class.
    %                                          Unlike most disturbances, this must satisfy:
    %                                          length(chan_in) == 1 && length(chan_in{1}) == 1
    %       omega : positive double :: frequency defining the band [-omega, omega] wherein the signal is white (for discrete-time systems, omega <= pi)
    %       horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of Disturbance properties
    %     Output:
    %       this_disturbance : DisturbanceBandedWhite object 
    %
    %     See also DisturbanceBandedWhite, Disturbance.Disturbance.
    
        % Defining defaults for missing arguments
        switch nargin
            case 0
                error('DisturbanceBandedWhite:DisturbanceBandedWhite',...
                      'Must provide at least the name of the disturbance')
            case 1
                chan_in = {1};
                omega = pi;
                horizon_period = [0, 1];
            case 2
                omega = pi;
                horizon_period = [0, 1];
            case 3
                horizon_period = [0, 1];
            case 4
        end
        
        % Checking conditions on chan_in
        assert(length(chan_in) == 1,... 
               'DisturbanceBandedWhite:DisturbanceBandedWhite',...
               ['chan_in must be constant for each timestep, ',...
                'provide chan_in as a cell of length 1'])
        assert(length(chan_in{1}) == 1,...
               'DisturbanceBandedWhite:DisturbanceBandedWhite',...
               ['chan_in must only specify one channel, ',...
                'provide chan_in such that chan_in{1} is a scalar'])
        % Calling Disturbance constructor
        this_dis@Disturbance(name, chan_in, horizon_period);
        
        % Setting specialized properties of DisturbanceBandedWhite
        validateattributes(omega, 'numeric', {'scalar', 'positive'})
        this_dis.omega = omega;
        
        this_dis = matchHorizonPeriod(this_dis);
    end

    function disp(this_dis)
    %% DISP function for DisturbanceBandedWhite object
    %
    %     disp(dis_bw_obj) (e.g., disp(DisturbanceBandedWhite('dis')) )
    %
    %     Variables:
    %     ---------
    %       Input:
    %          this_dis : DisturbanceBandedWhite object    
    %
    %     See also Ulft.disp, SequenceDisturbance.disp, Disturbance.disp
        disp@Disturbance(this_dis, 'Banded white')
        fprintf('%16s which is white between the band [-%3.1f, %3.1f]\n',...
                '',...
                this_dis.omega, this_dis.omega)
    end

    function this_dis = matchHorizonPeriod(this_dis, new_horizon_period)
    %% MATCHHORIZONPERIOD function to ensure properties of DisturbanceBandedWhite
    %  object match its own horizon_period, or a new_horizon_period
    %
    %  this_dis = matchHorizonPeriod(this_dis, new_horizon_period) will change the horizon_period and associated properties of this_dis
    %  this_dis = matchHorizonPeriod(this_dis) only ensures that each pertinent property of this_dis matches this_dis.horizon_period
    %
    %  Variables:
    %  ---------
    %     Input:
    %       this_dis : DisturbanceBandedWhite object
    %       new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence
    %     Output:
    %       this_dis : DisturbanceBandedWhite object
    %
    %  See also DisturbanceBandedWhite.
    
    if nargin == 1
    % Ensuring that this_dis.horizon_period matches with other properties 
    % of this_dis
    
        % Assumes properties are a mix of sequences of length 1 or
        % horizon_period
        total_time = sum(this_dis.horizon_period);
        if length(this_dis.chan_in) ~= total_time
            assert(length(this_dis.chan_in) == 1,...
                   'DisturbanceBandedWhite:matchHorizonPeriod',...
                   'channels of %s is not compatible w/ horizon_period',...
                   this_dis.name);
            this_dis.chan_in = repmat(this_dis.chan_in, 1, total_time);
        end
    else
    % Changing this_dis.horizon_period and other properties of this_dis to
    % a new horizon_period
    
        [indices, new_horizon_period] = ...
            makeNewIndices(this_dis.horizon_period, new_horizon_period);

        % Set properties according to indices
        this_dis.chan_in        = this_dis.chan_in(indices);
        this_dis.horizon_period = new_horizon_period;
        this_dis                = matchHorizonPeriod(this_dis);
    end
    end
    
    function multiplier = disturbanceToMultiplier(this_dis, varargin)
    %% DISTURBANCETOMULTIPLIER function to generate a multiplier from this object. 
    %
    %  multiplier = disturbanceToMultiplier(this_dis, 'dim_in_lft', dim_in_lft)
    %
    %  Variables:
    %  ---------
    %     Input:
    %       this_dis : DisturbanceBandedWhite object
    %      dim_in_lft : row of naturals :: input dimension of LFT associated with disturbance
    %     Output:
    %       multiplier : MultiplierBandedWhite object
    %
    %  See also DisturbanceBandedWhite
    input_parser = inputParser;
    addRequired(input_parser,...
                'disturbance',...
                @(dis) validateattributes(dis,...
                                          'DisturbanceBandedWhite',...
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
                 'discrete',...
                 true,...
                 @(disc) validateattributes(disc, 'logical', {'nonempty'}))
    addParameter(input_parser,...
                 'poles',...
                 -0.5,...
                 @(poles) validateattributes(poles,...
                                             'numeric',...
                                             {'real', 'vector'}));
    parse(input_parser, this_dis, varargin{:})
    disturbance = input_parser.Results.disturbance;
    dim_in_lft  = input_parser.Results.dim_in_lft;
    poles       = input_parser.Results.poles;
    discrete    = input_parser.Results.discrete;
    if isempty(dim_in_lft)
        error('DisturbanceBandedWhite:disturbanceToMultiplier',...
              ['Input arguments to disturbanceToMultiplier must include',...
               ' a dim_in_lft as a key/value pair. For example: \n',...
               'disturbanceToMultiplier(dis, ''dim_in_lft'', 1,...)'])
    end
    if ismember('discrete', input_parser.UsingDefaults)
        warning('DisturbanceBandedWhite:disturbanceToMultiplier',...
                ['Input arguments do not specify if the multiplier is discrete',...
                 '-time. Assuming as a default that the multiplier is discrete',...
                 'time. This can be specified with a discrete key/value pair.',...
               '\n For example: disturbanceToMultiplier(dis, ''discrete'', true,...)'])
    end
    multiplier = MultiplierBandedWhite(disturbance, dim_in_lft, discrete,...
                                       'poles', poles);
    end
end
end

%%  CHANGELOG
% Nov. 23, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)