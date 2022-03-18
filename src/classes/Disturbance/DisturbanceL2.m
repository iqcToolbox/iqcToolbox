classdef DisturbanceL2 < Disturbance
%% DISTURBANCEL2 class for the set of L2 signals. This extends the base 
% class Disturbance.
%
%   extended methods:
%     DisturbanceL2(name, chan_in, horizon_period) :: Constructor
%     disp(this_dis) :: Display method
%     matchHorizonPeriod(this_perf, horizon_period) :: Matches disturbance properties to new horizon_period
%     disturbanceToMultiplier(this_dis, 'dim_in_lft', dim_in_lft) :: Generate multiplier from disturbance
%
%  See also DisturbanceL2.DisturbanceL2, Disturbance

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

methods
    function this_dis = DisturbanceL2(name, chan_in, horizon_period)
    %% DISTURBANCEL2 constructor
    %
    %  dis = DisturbanceL2(name, chan_in, horizon_period)
    %  dis = DisturbanceL2(name, chan_in) assumes horizon_period == [0, 1]
    %  dis = DisturbanceL2(name) also assumes chan_in == {[]}
    %
    %  Variables:
    %  ---------
    %     Input:
    %       name : char array :: unique ID of disturbance signal (ex. 'sensor_noise')
    %       chan_in : cell array of naturals :: channels of signals pertaining to disturbance class
    %       horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of Disturbance properties
    %     Output:
    %       this_disturbance : DisturbanceL2 object 
    %
    %     See also DisturbanceL2, Disturbance.Disturbance.
    
        % Defining defaults for missing arguments
        switch nargin
            case 1
                chan_in = {[]};
                horizon_period = [0, 1];
            case 2
                horizon_period = [0, 1];
        end
        % Calling Disturbance constructor
        this_dis@Disturbance(name, chan_in, horizon_period);
        this_dis = matchHorizonPeriod(this_dis);
    end

    function disp(this_dis)
    %% DISP function for DisturbanceL2 object
    %
    %     disp(dis_l2_obj) (e.g., disp(DisturbanceL2('dis')) )
    %
    %     Variables:
    %     ---------
    %       Input:
    %          this_dis : DisturbanceL2 object    
    %
    %     See also Ulft.disp, SequenceDisturbance.disp, Disturbance.disp
        disp@Disturbance(this_dis, 'L2')        
    end

    function this_dis = matchHorizonPeriod(this_dis, new_horizon_period)
    %% MATCHHORIZONPERIOD function to ensure properties of DisturbanceL2
    %  object match its own horizon_period, or a new_horizon_period
    %
    %  this_dis = matchHorizonPeriod(this_dis, new_horizon_period) will change the horizon_period and associated properties of this_dis
    %  this_dis = matchHorizonPeriod(this_dis) only ensures that each pertinent property of this_dis matches this_dis.horizon_period
    %
    %  Variables:
    %  ---------
    %     Input:
    %       this_dis : DisturbanceL2 object
    %       new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence
    %     Output:
    %       this_dis : DisturbanceL2 object
    %
    %  See also DisturbanceL2.
    
    if nargin == 1
    % Ensuring that this_dis.horizon_period matches with other properties 
    % of this_dis
    
        % Assumes properties are a mix of sequences of length 1 or
        % horizon_period
        total_time = sum(this_dis.horizon_period);
        if length(this_dis.chan_in) ~= total_time
            assert(length(this_dis.chan_in) == 1,...
                   'DisturbanceL2:matchHorizonPeriod',...
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
        this_dis.chan_in       = this_dis.chan_in(indices);
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
    %       this_dis : DisturbanceL2 object
    %      dim_in_lft : row of naturals :: input dimension of LFT associated with disturbance
    %     Output:
    %       multiplier : MultiplierL2 object
    %
    %  See also DisturbanceL2
    input_parser = inputParser;
    addRequired(input_parser,...
                'disturbance',...
                @(dis) validateattributes(dis,...
                                          {'DisturbanceL2'},...
                                          {'nonempty'},...
                                          mfilename));
    addParameter(input_parser,...
                'dim_in_lft',...
                [],...
                @(dim) validateattributes(dim,...
                                          {'numeric'},...
                                          {'integer', 'row', 'positive'},...
                                          mfilename));
    addParameter(input_parser,... % This parameter is not used. Only defined for compatibility with other disturbanceToMultiplier calls
                 'discrete',...
                 true,...
                 @(disc) validateattributes(disc, {'logical'}, {'nonempty'}))

    parse(input_parser, this_dis, varargin{:})
    disturbance = input_parser.Results.disturbance;
    dim_in_lft  = input_parser.Results.dim_in_lft;
    if isempty(dim_in_lft)
        error('DisturbanceL2:disturbanceToMultiplier',...
              ['Input arguments to disturbanceToMultiplier must include',...
               ' a dim_in_lft as a key/value pair. For example: \n',...
               'disturbanceToMultiplier(dis, ''dim_in_lft'', 1,...)'])
    end
    multiplier = MultiplierL2(disturbance, dim_in_lft);
    end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)