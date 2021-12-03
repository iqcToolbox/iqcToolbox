classdef (InferiorClasses = {?ss}) Ulft
%% ULFT class for (upper) linear fractional transformation objects.
%
%     methods:
%       lft = Ulft(a, b, c, d, delta, varargin) :: Constructor
%       disp(this_lft) :: Display method
%       lft = addDisturbance(lft, disturbance) :: Appends the given cell-array of disturbances to lft
%       lft = addPerformance(lft, performance) :: Appends the given cell-array of performances to lft
%       lft = blkdiag(lft1, lft2, lft3, ...) :: block-diagonally concatenate multiple lfts
%       lft = gatherLft(lft) :: Gathering multiple copies of the same Delta/Disturbance/Performance blocks into one (typically back-end)
%       lft = horzcat(lft1, lft2, lft3, ...) :: Horizontally concatenate multiple lfts
%       lft = interconnect(upper_lft, lower_lft) :: Form one lft from two lfts connected in an upper lft interconnection
%       lft = inv(lft) :: Invert lft (if invertible)
%       lft = matchHorizonPeriod(lft, new_horizon_period) :: Ensures horizon_period of lft is consistent with its properties (typically back-end)
%       lft = minus(left_lft, right_lft) :: Subtraction operation between lfts
%       lft = mpower(lft, exponent) :: Power operation for lfts
%       lft = mrdivide(left_lft, right_lft) :: Division operation between lfts
%       lft = mtimes(left_lft, right_lft) :: Multiplication operation between lfts
%       lft = normalizeLft(lft) :: normalizes the Deltas of a lft
%       lft = plus(left_lft, right_lft) :: Addition operation between lfts
%       lft = power(lft, exponent) :: .^ Power operation for lfts
%       lft = random(varargin) :: Constructs a random / arbitrary lft
%       lft = rdivide(left_lft, right_lft) :: ./ Division operation between lfts
%       lft = removeUncertainty(lft, delta) :: Remove specified deltas from lft
%       lft = removeDisturbance(lft, disturbance) :: Remove specified disturbances from lft
%       lft = removePerformance(lft, performance) :: Remove specified performances from lft
%       lft = reorderLftDelta(lft, new_order) :: Reordering deltas within lft
%       lft = reorderLftDisturbance(lft, new_order) :: Reordering disturbances within lft
%       lft = reorderLftPerformance(lft, new_order) :: Reordering performances within lft
%       lft = sampleDeltas(lft, deltas, values, varargin) :: realizes the specified deltas of lft to specified or random values
%       [output, time, state] = simulate(lft, input, time, initial_state)
%       [dim_out, dim_in] = size(lft) :: Provide input/output dimensions of lft
%       lft = times(left_lft, right_lft) :: .* Multiplication operation between lfts
%       lft = uminus(lft) :: Negative operator for an lft
%       lft = uplus(lft) :: Positive operator for an lft
%       lft = vertcat(lft1, lft2, lft3, ...) :: Vertically concatenate multiple lfts
%
%     properties:
%       a : cell array of double matrices :: sequence of state/delta matrices
%       b : cell array of double matrices :: sequence of input matrices
%       c : cell array of double matrices :: sequence of output matrices
%       d : cell array of double matrices :: sequence of feedthrough matrices
%       delta : SequenceDelta object :: sequence of uncertainties in LFT
%       horizon_period : 1 x 2 array of naturals :: [horizon, period] (in
%                       timesteps) of LFT a, b, c, and d matrices
%       disturbance : SequenceDisturbance object :: sequence of LFT disturbances
%       performance : SequencePerformance object :: sequence of LFT performances
%       uncertain : boolean :: Indicator if LFT is uncertain or not
%       timestep : row of doubles (possibly empty) :: Timestep of LFT. Empty if memoryless,
%                                                                      0 if continuous-time,
%                                                                      -1 if discrete-time and unspecified timestep,
%                                                                      Positive scalars express the discrete-time system's timestep
%
%     See also Ulft.Ulft, toLft.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    a (1,:) cell
    b (1,:) cell
    c (1,:) cell
    d (1,:) cell
    delta SequenceDelta
    horizon_period (1, 2) double {mustBeNonnegative, mustBeInteger}
    disturbance SequenceDisturbance
    performance SequencePerformance
end

properties (Dependent)
    uncertain
    timestep
end

methods (Static)
    this_lft = random(varargin)
end

methods
this_lft = gatherLft(this_lft)
this_lft = matchHorizonPeriod(this_lft, new_horizon_period)
this_lft = reorderLftDelta(this_lft, new_order)
this_lft = reorderLftDisturbance(this_lft, new_order)
this_lft = reorderLftPerformance(this_lft, new_order)
this_lft = interconnect(upper_lft, lower_lft)
this_lft = addDisturbance(this_lft, disturbance)
this_lft = addPerformance(this_lft, performance)
this_lft = removeUncertainty(this_lft, delta)
this_lft = removeDisturbance(this_lft, disturbance)
this_lft = removePerformance(this_lft, performance)
this_lft = normalizeLft(this_lft)
this_lft = sampleDeltas(this_lft, deltas, values, varargin)
[output, time, state] = simulate(lft, input, time, initial_state)

% Overloading methods
this_lft = uplus(this_lft)
this_lft = plus(left_lft, right_lft)
this_lft = uminus(this_lft)
this_lft = minus(left_lft, right_lft)
this_lft = inv(this_lft)
this_lft = mrdivide(left_lft, right_lft)
this_lft = rdivide(left_lft, right_lft)
this_lft = times(left_lft, right_lft)
this_lft = mtimes(left_lft, right_lft)
this_lft = mpower(this_lft, exponent)
this_lft = power(this_lft, exponent)
this_lft = horzcat(varargin)
this_lft = vertcat(varargin)
this_lft = blkdiag(varargin)
disp(this_lft)

function this_lft = Ulft(a, b, c, d, delta, varargin)
%% ULFT constructor
%
%     this_lft = Ulft(a, b, c, d, delta, 'horizon_period', horizon_period, 'disturbance', disturbance, 'performance', performance)
%     this_lft = Ulft(a, b, c, d, delta) assumes horizon_period == [0, 1],
%                                           performance == l2-induced performance (default)
%                                           disturbance == l2 signals (default).
%     If specifying any but not all of the three aforementioned properties,
%        the unspecified properties assume the default values.
%        e.g., this_lft = Ulft(a, b, c, d, delta, 'horizon_period', [2, 3])
%                    creates an lft with defaults assumed for 'performance' and 'disturbance'
%
%     Variables:
%     ---------
%       Input:
%          a : double matrix or cell array of double matrices :: sequence of state/delta matrices
%          b : double matrix or cell array of double matrices :: sequence of input matrices
%          c : double matrix or cell array of double matrices :: sequence of output matrices
%          d : double matrix or cell array of double matrices :: sequence of feedthrough matrices
%          delta : Delta object, or cell array of Delta objs, or SequenceDelta object :: seq of uncertainties in LFT
%          horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of LFT a,b,c,d matrices
%          disturbance : Disturbance object, or cell array of Disturbance objs, or SequenceDisturbance object :: seq of LFT disturbances
%          performance : Performance object, or cell array of Performance objs, or SequencePerformance object :: seq of LFT performance
%       Output:
%          this_lft : Ulft object :: the produced lft
%
%       See also Ulft, toLft.

%% PARSE INPUTS
input_parser = inputParser;

% a, b, c, d matrices
isMatrix = @(m) isnumeric(m) && ismatrix(m);
isCellOfMatrices = @(m) iscell(m) && all(cellfun(isMatrix, m)) && size(m, 1)==1;
validABCD = @(m) isMatrix(m) || isCellOfMatrices(m);
addRequired(input_parser, 'a', validABCD);
addRequired(input_parser, 'b', validABCD);
addRequired(input_parser, 'c', validABCD);
addRequired(input_parser, 'd', validABCD);

% delta
isDelta = @(del) isa(del, 'Delta');
isCellOfDeltas = @(del) iscell(del) && all(cellfun(isDelta,del));
validDelta = @(del) isDelta(del)...
                    || isCellOfDeltas(del)...
                    || isa(del, 'SequenceDelta');
addRequired(input_parser, 'delta', validDelta);

% horizon_period
validHorizonPeriod = @(hp) validateattributes(hp,...
                                              'numeric',...
                                              {'size',[1,2],...
                                               'integer',...
                                               'nonnegative'});
addParameter(input_parser, 'horizon_period', [0, 1], validHorizonPeriod);

% disturbance
isDisturbance = @(dis) isa(dis, 'Disturbance');
isCellOfDisturbances = @(dis) iscell(dis) && all(cellfun(isDisturbance,dis));
validDisturbance = @(dis) isDisturbance(dis)...
                          || isCellOfDisturbances(dis)...
                          || isa(dis, 'SequenceDisturbance');
addParameter(input_parser,...
             'disturbance',...
             DisturbanceL2('default_l2'),...
             validDisturbance);

% performance
isPerformance = @(perf) isa(perf, 'Performance');
isCellOfPerformances = @(perf) iscell(perf) && ...
                               all(cellfun(isPerformance,perf));
validPerformance = @(perf) isPerformance(perf)...
                           || isCellOfPerformances(perf)...
                           || isa(perf, 'SequencePerformance');
addParameter(input_parser,...
             'performance',...
             PerformanceL2Induced('default_l2'),...
             validPerformance);

% Parse everything
parse(input_parser, a, b, c, d, delta, varargin{:})

horizon_period = input_parser.Results.horizon_period;

%% ENSURE THAT A, B, C, D ARE CELLS
% a, b, c, d
if ~iscell(input_parser.Results.a)
    a = {input_parser.Results.a};
else
    a = input_parser.Results.a;
end
if ~iscell(input_parser.Results.b)
    b = {input_parser.Results.b};
else
    b = input_parser.Results.b;
end
if ~iscell(input_parser.Results.c)
    c = {input_parser.Results.c};
else
    c = input_parser.Results.c;
end
if ~iscell(input_parser.Results.d)
    d = {input_parser.Results.d};
else
    d = input_parser.Results.d;
end

%% ENSURE THAT DELTA, PERFORMANCE, AND DISTURBANCE ARE SEQUENCE*
% delta, performance, disturbance
if ~isa(input_parser.Results.delta, 'SequenceDelta')
    delta = SequenceDelta(input_parser.Results.delta);
else
    delta = input_parser.Results.delta;
end
if ~isa(input_parser.Results.performance, 'SequencePerformance')
    performance = SequencePerformance(input_parser.Results.performance);
else
    performance = input_parser.Results.performance;
end
if ~isa(input_parser.Results.disturbance, 'SequenceDisturbance')
    disturbance = SequenceDisturbance(input_parser.Results.disturbance);
else
    disturbance = input_parser.Results.disturbance;
end

%% CHECK CONSISTENCY BETWEEN INPUTS

% horizon_period against a, b, c, d
total_time = sum(horizon_period);
validateattributes(a, 'cell', {'numel', total_time})
validateattributes(b, 'cell', {'numel', total_time})
validateattributes(c, 'cell', {'numel', total_time})
validateattributes(d, 'cell', {'numel', total_time})

% horizon_period against delta, performance, disturbance
delta = matchHorizonPeriod(delta, horizon_period);
performance = matchHorizonPeriod(performance, horizon_period);
disturbance = matchHorizonPeriod(disturbance, horizon_period);

% delta against a
cellOfEmpty = @(a) all(cellfun('size', a, 1) == 0);
if cellOfEmpty(a) && isempty(delta.deltas)
elseif cellOfEmpty(a)
    error('Ulft:Ulft',...
          ['A matrices are empty,',...
           ' which is inconsistent with a non-empty Delta'])
elseif isempty(delta.deltas)
    error('Ulft:Ulft',...
          ['Delta is empty,',...
           ' which is inconsistent with a non-empty A matrix'])
else
    aOutMatchesDelIn = all(cellfun('size',a, 1) == sum(delta.dim_ins, 1));
    aInMatchesDelOut = all(cellfun('size',a, 2) == sum(delta.dim_outs, 1));
    aConsistentWithDel = aOutMatchesDelIn && aInMatchesDelOut;
    assert(aConsistentWithDel,...
           'Ulft:Ulft',...
           'Dimensions of A matrices arent consistent with Delta')
end

% disturbance against b and d
allChannelsEmpty = @(chan) all(cellfun(@isempty, chan), 'all');
matInputContainsChannels = @(chan, mat) ...
    all(cellfun('size', mat, 2) >= max(cellfun(@maxEmpty, chan), [], 1));
bConsistentWithDis = allChannelsEmpty(disturbance.chan_ins)...
                     ||...
                     matInputContainsChannels(disturbance.chan_ins, b);
dConsistentWithDis = allChannelsEmpty(disturbance.chan_ins)...
                     ||...
                     matInputContainsChannels(disturbance.chan_ins, d);
assert(bConsistentWithDis && dConsistentWithDis,...
       'Ulft:Ulft',...
       ['Disturbance channels are not consistent with'...
        ' dimensions of B matrices or D matrices'])

% performance against b, c and d
matOutputContainsChannels = @(chan, mat) ...
    all(cellfun('size', mat, 1) >= max(cellfun(@maxEmpty, chan), [], 1));

bConsistentWithPerf = allChannelsEmpty(performance.chan_ins)...
                      || ...
                      matInputContainsChannels(performance.chan_ins, b);
cConsistentWithPerf = allChannelsEmpty(performance.chan_outs)...
                      || ...
                      matOutputContainsChannels(performance.chan_outs, c);
dConsistentWithPerfOut = allChannelsEmpty(performance.chan_outs)...
                         || ...
                         matOutputContainsChannels(performance.chan_outs, d);
dConsistentWithPerfIn = allChannelsEmpty(performance.chan_ins)...
                        || ...
                        matInputContainsChannels(performance.chan_ins, d);
assert(bConsistentWithPerf && cConsistentWithPerf && ...
       dConsistentWithPerfOut && dConsistentWithPerfIn,...
       'Ulft:Ulft',...
       ['Dimensions of performances are not consistent with'...
        ' dimensions of C matrices or D matrices'])

% b against a and d
bConsistentWithA = all(cellfun('size', b, 1) == cellfun('size', a, 1));
bConsistentWithD = all(cellfun('size', b, 2) == cellfun('size', d, 2));
assert(bConsistentWithA && bConsistentWithD ,...
       'Ulft:Ulft',...
       ['Dimensions of B matrices are not consistent with'...
        ' dimensions of A matrices or D matrices'])

% c against a and d
cConsistentWithA = all(cellfun('size', c, 2) == cellfun('size', a, 2));
cConsistentWithD = all(cellfun('size', c, 1) == cellfun('size', d, 1));
assert(cConsistentWithA && cConsistentWithD,...
       'Ulft:Ulft',...
       ['Dimensions of C matrices are not consistent with'...
        ' dimensions of A matrices or D matrices'])

% d against disturbance and performance

%% SET PROPERTIES
this_lft.a = a;
this_lft.b = b;
this_lft.c = c;
this_lft.d = d;
this_lft.horizon_period = horizon_period;
this_lft.delta = delta;
this_lft.disturbance = disturbance;
this_lft.performance = performance;

%% FINAL CONDITIONING OF LFT
this_lft = gatherLft(this_lft);
end

%% Getter methods for dependent properties

function timestep = get.timestep(this_lft)
    if isempty(this_lft.delta.deltas)
        timestep = [];
    else
        if isa(this_lft.delta.deltas{1}, 'DeltaIntegrator')
            timestep = 0;
        elseif isa(this_lft.delta.deltas{1}, 'DeltaDelayZ')
            timestep = this_lft.delta.deltas{1}.timestep;
        else
            timestep = [];
        end
    end
end

function uncertain = get.uncertain(this_lft)
    if ~isempty(this_lft.timestep)
    % this_lft has DeltaDelayZ or DeltaIntegrator
        uncertain = length(this_lft.delta.deltas) > 1;
    else
        uncertain = ~isempty(this_lft.delta.deltas);
    end
end

end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Added sampleDeltas method and .uncertain and .timestep properties - Jason Nezvadovitz and Micah Fry
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)