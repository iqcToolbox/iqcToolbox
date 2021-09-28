classdef (Abstract) MultiplierDisturbance < matlab.mixin.Heterogeneous
%% MULTIPLIERDISTURBANCE abstract base class for IQC multipliers of disturbances
%  This class should not be directly accessed. Subclasses of
%  MultiplierDisturbance (such as MultiplierL2, etc.) extend
%  and concretize this class. This class inherits from matlab.mixin.Heterogeneous 
%  to allow for object arrays of this class (and mixtures of subclasses)
%
%  methods:
%     getDefaultScalarElement() :: Necessary method for mixin.Heterogeneous subclasses
%  
%  properties:
%     name : char array :: unique ID of multiplier (same as Disturbance
%            signal)
%     filter : struct with fields a, b, c, d, each a cell of matrices 
%            :: state-space realization for Multiplier's filter
%     decision_vars : cell array of sdpvar objects (yalmip) ::
%                     decision_vars defining Multiplier
%     constraints : lmi object (yalmip) :: list of constraints on decision_vars
%     quad : struct with field q, a cell of sdpvar objs
%          :: structured collection of decision_vars, such that the IQC 
%             multiplier == filter' * quad * filter
%     horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) 
%                                                 of Multiplier properties
%
%  See also MultiplierPerformance, MultiplierDelta

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    name 
    filter struct
    decision_vars (1, :) cell
    constraints 
    quad struct
    horizon_period (1, 2) double {mustBeNonnegative, mustBeInteger}
end

methods (Static, Sealed, Access = protected)
    function default_multiplier = getDefaultScalarElement
        default_multiplier = MultiplierDisturbanceDefault();
    end
end   

end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)