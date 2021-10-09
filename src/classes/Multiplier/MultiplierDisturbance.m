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
%     filter_lft : Ulft object :: convenience getter property to express the filter in
%                                 less structured terms (a Ulft object)
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
    discrete logical
end

properties (Dependent)
   filter_lft
end

methods (Static, Sealed, Access = protected)
    function default_multiplier = getDefaultScalarElement
        default_multiplier = MultiplierDisturbanceDefault();
    end
end

methods
function filter_lft = get.filter_lft(this_mult)
filt = this_mult.filter;
if all(cellfun(@isempty, filt.a))
    del = SequenceDelta();
else
    if this_mult.discrete
        % timestep may be wrong, will need to correct this to concatenate
        timestep = -1; 
        del = DeltaDelayZ(cellfun(@(a) size(a, 2), filt.a),...
                          timestep,...
                          this_mult.horizon_period);
    else
        del = DeltaIntegrator(size(filt.a{1}, 2));
    end
end
filter_lft = Ulft(filt.a, filt.b, filt.c, filt.d, del,...
                  'horizon_period', this_mult.horizon_period);
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)
