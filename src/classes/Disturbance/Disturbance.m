classdef (Abstract) Disturbance
%% DISTURBANCE abstract base class for specifying disturbance classes
% This class should not be directly accessed. Subclasses of Disturbance
% (such as DisturbanceL2) extend and concretize this class.
%
%   concrete methods:
%     Disturbance(name, channel, horizon_period) :: Constructor method
%     disp(this_dis) :: Display method
%
%   abstract methods:
%     matchHorizonPeriod(this_dis, horizon_period) :: Matches disturbance properties to a new horizon_period
%     disturbanceToMultiplier(this_dis, varargin) ::  Method for constructing a multiplier from a disturbance
%
%   properties:
%     name : char array :: unique ID of disturbance signal (ex. 'sensor_noise')
%     channel : cell array of naturals :: channels of signals pertaining to disturbance class
%     horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) 
%                                                 of Disturbance properties
%
%   See also Disturbance.Disturbance

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    name (1, :) char
    chan_in (1, :) cell
    horizon_period (1, 2) double {mustBeInteger, mustBeNonnegative}
end

methods (Abstract)
    this_perf = matchHorizonPeriod(this_dis, horizon_period)
    multiplier = disturbanceToMultiplier(this_dis, varargin)
end

methods
    function this_dis = Disturbance(name, chan_in, horizon_period)
    %% DISTURBANCE constructor, must be called from a subclass constructor
    %
    %     this_dis = Disturbance(name, dim_in)
    %
    %     Variables:
    %     ---------
    %       Input:
    %           name : char array :: unique ID of disturbance signal (ex. 'sensor_noise')
    %           channel : cell array of naturals :: channels of signals for disturbance set
    %           horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of Disturbance properties
    %       Output:
    %           this_disturbance : Disturbance object :: the base class Disturbance
    %
    %     See also Disturbance.
    validateattributes(name, {'char'}, {'row'}, mfilename);
    validateattributes(chan_in, {'cell'}, {'row'}, mfilename);
    for i = 1:length(chan_in)
        validateattributes(chan_in{i},...
                               {'numeric'},...
                               {'integer', 'positive'},...
                               mfilename);
        assert(size(chan_in{i}, 2) < 2,...
               'Disturbance:Disturbance',...
               ['Specified input channels must be cell array of empty',...
                ' or column arrays'])
        assert(length(chan_in{i}) == length(unique(chan_in{i})),...
               'Disturbance:Disturbance',...
               'Specified channels contains duplicates')
    end        
    validateattributes(horizon_period, {'numeric'}, {'size', [1,2],...
                                                   'integer',...
                                                   'nonnegative'});
    validateattributes(horizon_period(2), {'numeric'}, {'positive'}) 
    
    % Setting properties of Disturbance
    this_dis.name           = name;
    this_dis.chan_in        = chan_in;
    this_dis.horizon_period = horizon_period;
    end

    function disp(this_dis, type)
    %% DISP function for Disturbance object, must be called from a subclass disp method.
    %
    %     disp(this_dis, type)
    %
    %     Variables:
    %     ---------
    %       Input:
    %           this_dis : Disturbance object     
    %           type : char array :: type of disturbance class ('banded', etc.)
    %
    %     See also Ulft.disp, SequenceDisturbance.disp

    fprintf('%4s %-7s is a %8s disturbance pertaining to',...
            '',...
            this_dis.name,...
            type)
    if isempty(this_dis.chan_in{1})
        fprintf('\n %16s all input channels \n', '')
    else
        chan = sprintf( '%d; ', this_dis.chan_in{1});
        fprintf(['\n %16s input channels [', chan(1 : end - 2), '] \n'],'')
    end
    end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)