classdef (Abstract) MultiplierPerformance < matlab.mixin.Heterogeneous
%% MULTIPLIERPERFORMANCE abstract base class for IQC multipliers of uncertainties
%  This class should not be directly accessed. Subclasses of MultiplierPerformance
%  (such as MultiplierL2Induced, etc.) extend and concretize
%  this class. This class inherits from matlab.mixin.Heterogeneous 
%  to allow for object arrays of this class (and mixtures of subclasses)
%
%  methods:
%     getDefaultScalarElement() :: Necessary method for mixin.Heterogeneous subclasses
%  
%   properties:
%     name : char array :: unique ID of multiplier (same as Performance class)
%     filter : struct with fields a, b1, b2, c1, c2, d11, d12, d21, d22,
%              each a cell of matrices 
%            :: state-space realization for Multiplier's filter
%     decision_vars : cell array of sdpvar objects (yalmip) ::
%                     decision_vars defining Multiplier
%     constraints : lmi object (yalmip) :: list of constraints on decision_vars
%     quad : struct with fields q11, q12, q21, q22, each a cell of sdpvar objs
%          :: structured collection of decision_vars, such that the IQC 
%             multiplier == filter' * quad * filter
%     horizon_period : 1 x 2 array of naturals :: [horizon, period] (in
%                    timesteps) of MultiplierPerformance properties
%     objective : double or sdpvar (yalmip) :: Objective to be minimized for this Performance
%     objective_scaling : double scalar :: scaling factor for object to scalarize
%                                          the vector optimization problem that would arise when
%                                          combining multiple performances
%     filter_lft : Ulft object :: convenience getter property to express the filter in
%                                 less structured terms (a Ulft object)
%
%  See also MultiplierDisturbance, MultiplierDelta
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
    objective
    objective_scaling double {mustBeReal, mustBeNonnegative}
    discrete logical
end

properties (Dependent)
    filter_lft
end

methods (Static, Sealed, Access = protected)
    function default_multiplier = getDefaultScalarElement
        default_multiplier = MultiplierPerformanceDefault();
    end
end        

methods
function filter_lft = get.filter_lft(this_mult)
filt = this_mult.filter;
a = filt.a;

% b and c matrices are trickier to get because some may be empty
total_time = sum(this_mult.horizon_period);
b = cell(1, total_time);
c = cell(1, total_time);
for i = 1:total_time
    % Commented block was necessary before changing combineAllMultipliers, 
    % maintaining old code in case of regressions
    if isempty(filt.b1{i}) && ~isempty(filt.b2{i})
        filt.b1{i} = zeros(size(filt.a{i}, 1), size(filt.b1{i}, 2));
    elseif isempty(filt.b2{i}) && ~isempty(filt.b1{i})
        filt.b2{i} = zeros(size(filt.a{i}, 1), size(filt.b2{i}, 2));
    end
    b{i} = [filt.b1{i}, filt.b2{i}];

    % Commented block was necessary before changing combineAllMultipliers, 
    % maintaining old code in case of regressions
    if isempty(filt.c1{i}) && ~isempty(filt.c2{i})
        filt.c1{i} = zeros(size(filt.c1{i}, 1), size(filt.a{i}, 2));
    elseif isempty(filt.c2{i}) && ~isempty(filt.c1{i})
        filt.c2{i} = zeros(size(filt.c2{i}, 1), size(filt.a{i}, 2));
    end
    c{i} = [filt.c1{i}; filt.c2{i}];
end

d = cellfun(@(d11, d12, d21, d22) [d11, d12; d21, d22],...
            filt.d11, filt.d12, filt.d21, filt.d22,...
            'UniformOutput', false);
if all(cellfun(@isempty, a))
    del = SequenceDelta();
else
    if this_mult.discrete
        % timestep may be wrong, will need to correct this to concatenate
        timestep = -1; 
        del = DeltaDelayZ(cellfun(@(a) size(a, 2), a),...
                          timestep,...
                          this_mult.horizon_period);
    else
        del = DeltaIntegrator(size(a{1}, 2));
    end
end
filter_lft = Ulft(a, b, c, d, del, 'horizon_period', this_mult.horizon_period);
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)