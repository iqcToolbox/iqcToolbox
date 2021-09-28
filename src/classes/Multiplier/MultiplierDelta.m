classdef (Abstract) MultiplierDelta < matlab.mixin.Heterogeneous
%% MULTIPLIERDELTA abstract base class for IQC multipliers of uncertainties
%  This class should not be directly accessed. Subclasses of MultiplierDelta
%  (such as MultiplierSlti, MultiplierBounded, etc.) extend and concretize
%  this class. This class inherits from matlab.mixin.Heterogeneous to allow for
%  object arrays of this class (and mixtures of subclasses)
%
%  methods:
%     getDefaultScalarElement() :: Necessary method for mixin.Heterogeneous subclasses
%
%  properties:
%     name : char array :: unique ID of multiplier (same as Delta uncertainty)
%     filter : struct with fields a, b1, b2, c1, c2, d11, d12, d21, d22,
%              each a cell of matrices 
%            :: state-space realization for Multiplier's filter
%     decision_vars : cell array of sdpvar objects (yalmip) ::
%                     decision_vars defining Multiplier
%     constraints : lmi object (yalmip) :: list of constraints on decision_vars
%     quad : struct with fields q11, q12, q21, q22, each a cell of sdpvar objs
%          :: structured collection of decision_vars, such that the IQC 
%             multiplier == filter' * quad * filter
%     horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) 
%                                                 of Multiplier properties
%
%  See also MultiplierPerformance, MultiplierDisturbance

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    name 
    filter
    decision_vars
    constraints
    quad
    horizon_period
end

methods (Static, Sealed, Access = protected)
    function default_multiplier = getDefaultScalarElement
    %% GETDEFAULTSCALARELEMENT method for filling unspecified elements in an object array
    %
    %  default_mult = MultiplierDelta.getDefaultScalarElement()
    %
    %  Variables:
    %  ---------
    %    Output:
    %       default_multiplier : MultiplierDeltaDefault object :: placeholder for undefined multiplier
    %
    %  See also MultiplierDelta
        default_multiplier = MultiplierDeltaDefault();
    end
end        

end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)