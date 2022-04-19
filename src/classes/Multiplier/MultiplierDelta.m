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
%     filter_lft : Ulft object :: convenience getter property to express the filter in
%                                 less structured terms (a Ulft object)
%     exponential : scalar double :: the scalar used to define the exponential
%                                     IQC (rho in discrete-time, alpha in continuous-time
%                                     per the literature). If no argument is given for
%                                     the default constructor, "exponential" is set
%                                     to [], which signals the default value:
%                                     0 for continuous-time or 1 for discrete-time.
%
%
%  See also MultiplierDelta/MultiplierDelta, MultiplierPerformance, MultiplierDisturbance

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

properties (SetAccess = protected)
    exponential double {mustBeNonnegative}
end

properties (Dependent)
    filter_lft
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

methods

function this_mult = MultiplierDelta(exponential)
%% MUTLIPLIERDELTA base constructor
%  This does not construct the entire multiplier. It only fills the exponential field
%  (for rho- or alpha-IQC multipliers). If not explicitly called in a subclass's
%  constructor, this contructor will be called before anything else in the subclass
%  constructor.  This design approach allows developers unaware of exponential IQCs
%  to design multiplier classes without this generalization in mind.
%
%  this_mult = this_mult@MultiplierDelta(exponential) sets this_mult.exponential to the provided value
%  this_mult = this_mult@MultiplierDelta() sets this_mult.exponential to [], which indicates 
%                                          the "default" value (1 for disrete-time, 0 for continuous-time)
%
%  Variables:
%  ---------
%    Input:
%       exponential : double scalar :: scalar indicating exponential rate of exponential IQCs
%                                      this is "rho" for discrete-time IQCs (rho ^ -2k),
%                                      and "alpha" for continuous-time IQCs (exp ^ -2*alpha*t)
%
%    Output:
%       this_mult : MultiplierDelta object :: the base class MultiplierDelta
%     
%  See also MultiplierDelta, 
%           MultiplierSlti and MultiplierSltv (examples on exponential-agnostic multipliers)
%           MultiplierDlti (example of a mult that rejects non-default exponential-rates)
%           MultiplierConstantDelay and MultiplierConstantDelay2 (examples of mults whose cnstraints depend on exponential-rates)

if nargin == 0
    this_mult.exponential = [];
elseif nargin == 1
    this_mult.exponential = exponential;
end
end

function this_mult = shiftMultiplier(this_mult, exponential)
%% SHIFTMULTIPLIER method for MultiplierDelta objects.
%  This standardizes how all multiplier shifting is done for any subclass of MultiplierDelta.
%  It is expected that this method would be called near the end of the subclass's
%  constructor method (after this_mult.filter and this_mult.discrete are assigned in that constructor).
%  This will shift the multiplier's filter depending on the multiplier being discrete- or continuous-time
%  After shifting, this method checks that the multiplier is still stable.
%
%  this_mult = shiftMultiplier@MultiplierDelta(this_mult, exponential) % Shifts state-space matrices of this_mult.filter, and sets this_mult.exponential
%  this_mult = shiftMultiplier@MultiplierDelta(this_mult) % Shifts state-space matrices according to pre-defined this_mult.exponential
%
%  Variables:
%  ---------
%    Input:
%       exponential : double scalar :: scalar indicating exponential rate of exponential IQCs
%                                      this is "rho" for discrete-time IQCs (rho ^ -2k),
%                                      and "alpha" for continuous-time IQCs (exp ^ 2*alpha*t)
%
%    Output:
%       this_mult : MultiplierDelta object :: the shifted MultiplierDelta* object
%
%  See also MultiplierDelta, MultiplierDelta.MultiplierDelta

    % First check if multiplier has continuous/discrete-time dynamics, quit if not
    if isempty(this_mult.discrete)
        return
    end
    % If exponential is not given, use the pre-existing field in the object
    if nargin == 1
        exponential = this_mult.exponential;
    end
    % If exponential is empty (denoting default), set to 1 for discrete-time, 0 for continuous-time
    if isempty(exponential) 
        if this_mult.discrete
            exponential = 1;
        else
            exponential = 0;
        end
    end
    % Shift A-matrices of filter
    if this_mult.discrete
        a_mat = cellfun(@(a) a / exponential, this_mult.filter.a,...
                        'UniformOutput', false);
    else
        a_mat = cellfun(@(a) a + eye(size(a)) * exponential, this_mult.filter.a,...
                        'UniformOutput', false);
    end
    this_mult.filter.a = a_mat;
    % Shift B-matrices of filter
    if this_mult.discrete
        b1_mat = cellfun(@(b) b / exponential, this_mult.filter.b1,...
                        'UniformOutput', false);
        b2_mat = cellfun(@(b) b / exponential, this_mult.filter.b2,...
                        'UniformOutput', false);
        this_mult.filter.b1 = b1_mat;
        this_mult.filter.b2 = b2_mat;
    end
    % Check that filter remains stable
    if sum(this_mult.horizon_period) == 1
    % Time-invariant systems
        if this_mult.discrete
            stable = max(abs(eig(this_mult.filter.a{1}))) < 1;
        else
            stable = max(real(eig(this_mult.filter.a{1}))) < 0;
        end
%     else
% Uncomment this when functionality for exponential stability is enabled for LTV systems
%         filt_lft = addPerformance(this_mult.filter_lft,...
%                                     {PerformanceStable(this_mult.horizon_period)});
%         result = iqcAnalysis(filt_lft);
%         stable = result.valid;
    end
    assert(stable, 'MultiplierDelta:shiftMultiplier',...
           ['Shifting Multiplier ', this_mult.name, ' by the exponential rate ',...
            exponential, ' produces an unstable multiplier.  You must construct',...
            ' a Multiplier and provide an exponential rate such that the',...
            ' shifted Multiplier remains stable'])
    this_mult.exponential = exponential;
end

function filter_lft = get.filter_lft(this_mult)
filt = this_mult.filter;
a = filt.a;
b = cellfun(@(b1, b2) [b1, b2], filt.b1, filt.b2, 'UniformOutput', false);
c = cellfun(@(c1, c2) [c1; c2], filt.c1, filt.c2, 'UniformOutput', false);
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