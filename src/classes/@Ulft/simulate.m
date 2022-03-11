function [output, time, state] = simulate(this_lft, input, time, initial_state)
%% SIMULATE method for simulating a certain Ulft (discrete-time LTV, or continuous-time LTI) with 
%    given input signals and initial conditions
%
%    [output, time, state] = simulate(this_lft, input, time, initial_state) 
%    output = simulate(this_lft, input, time) assumes initial_state == 0.  time must be defined for continuous-time systems
%    output = simulate(this_lft, input)  for discrete-time systems (time should NOT be specified)
%    output = simulate(this_lft, input, [], initial_state) for specifying the initial_state of a discrete-time Ulft
%
%    Note: The calling syntax for this function is designed to mostly match lsim
%
%    Variables:
%    ---------
%      Inputs:
%        input : matrix or cell array of doubles :: the input signal to be fed into the Ulft
%        time : array of doubles :: the vector of time indices to which the input signal pertains. 
%                                   This should only be specified for continuous-time systems
%        initial_state : n x 1 array of doubles :: the initial state of the Ulft
%
%      Output:
%        output : array of doubles :: the output signal from the Ulft
%        time : array of doubles :: the vector of time indices to which the input signal pertains
%        state : array of doubles :: the state of the Ulft during simulation
%
%    See also Ulft.sampleDeltas

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Process inputs
switch nargin
    case 1
        error('Ulft:simulate', 'Must provide at least an input trajectory')
    case 2
        time = [];
        initial_state = [];
    case 3
        initial_state = [];
end

%% Validate inputs
assert(isempty(time) || length(time) == size(input, 2),...
       'Ulft:simulate',...
       ['The given time sequence must match the length of the input signal, '...
        'or may be empty for memoryless/discrete-time systems'])

assert(~this_lft.uncertain,...
       'Ulft:simulate',...
       ['Cannot simulate an uncertain system. Consider sampling uncertainties',...
       ' first by using Ulft.sampleDeltas'])

% Check that input matches lft dimensions
dim_in = size(this_lft, 2);
if all(dim_in == dim_in(1))
    % If lft dimensions are constant
    validateattributes(input, {'numeric'}, {'2d', 'nrows', dim_in(1)})
else
    % If lft dimensions are time-varying, need to map sequence to new indices
    total_time = size(input, 2);
    horz_per = [total_time - 1, 1];
    horz_per = commonHorizonPeriod([this_lft.horizon_period; horz_per]);
    time_indices = makeNewIndices(this_lft.horizon_period, horz_per);
    time_indices = time_indices(1 : total_time);
    % Check matching dimensions
    validateattributes(input, {'cell'}, {'nonempty'})
    dim_match = all(dim_in(time_indices) == cellfun(@length, input));
    assert(dim_match, 'Ulft:simulate', 'Dims of input and Ulft do not match')
end

%% Call different sub-functions if system is LTI, Discrete-time LTV, or Memoryless LTV
if isequal(this_lft.horizon_period, [0, 1])
    [output, time, state] = timeInvariantSim(this_lft, input, time, initial_state);
else
    if ~isempty(this_lft.timestep)
        [output, time, state] = discreteTimeVaryingSim(this_lft, input, time, initial_state);    
    else
        [output, time, state] = noStateTimeVaryingSim(this_lft, input);    
    end
end
    

end

function [output, time, state] = timeInvariantSim(this_lft, input, time, initial_state)
%% Simulate an LTI system
sys = lftToSs(this_lft);
% If lft is memoryless, have the ss specified as discrete-time
if isempty(this_lft.timestep)
    sys.Ts = -1;
end
[output, time, state] = lsim(sys, input, time, initial_state);
% Match format of argout with other lft types
output = output';
time = time';
state = state';
end

function [output, time, state] = discreteTimeVaryingSim(this_lft, input, time, initial_state)
%% Simulate a Discrete-time LTV system

% Obtain time_indices which map this_lft properties to input time
total_time = size(input, 2);
horizon_period = [total_time - 1, 1];
horizon_period = commonHorizonPeriod([this_lft.horizon_period; horizon_period]);
time_indices = makeNewIndices(this_lft.horizon_period, horizon_period);
time_indices = time_indices(1 : total_time);

% Calculate correct time, assert the given time matches
if any(this_lft.timestep == -1)
    correct_time = 0 : total_time - 1;
else
    correct_time = cumsum(this_lft.timestep(time_indices)) - this_lft.timestep(1);
end
assert(isempty(time) || all(abs(time - correct_time) < 1e-8) ,...
       'Ulft:simulate',...
       ['The given time sequence is incorrect. For discrete-time time-varying',...
        ' systems, this may be given as an empty ([])'])
time = correct_time;

% Check that initial_state is correct
dim_state = this_lft.delta.deltas{1}.dim_out;
if isempty(initial_state)
    initial_state = zeros(dim_state(1), 1);
end
assert(size(initial_state, 1) == dim_state(1),...
       'Ulft:simulate',...
       'The given initial_state has the wrong size')

% Convert input to cell to allow for uniform treatment
if ~iscell(input)
    input = mat2cell(input, size(input, 1), ones(1, total_time));
end
output = cell(1, total_time);
state = cell(1, total_time);
state{1} = initial_state;
a = this_lft.a(time_indices);
b = this_lft.b(time_indices);
c = this_lft.c(time_indices);
d = this_lft.d(time_indices);
% Though time starts at zero, matlab doesn't allow zero-indexed arrays, so we start at i = 1
for i = 1:total_time - 1
    state{1, i + 1} = a{i} * state{1, i} + b{i} * input{1, i};
    output{1, i}    = c{i} * state{1, i} + d{i} * input{1, i};
end
output{1, end} = c{end} * state{1, end} + d{end} * input{1, end};

% Bring output and state to matrices if not needed as cell arrays
dim_out = size(this_lft, 1);
if all(dim_out(1) == dim_out)
    output = cell2mat(output);
end
if all(dim_state(1) == dim_state)
    state = cell2mat(state);
end
end

function [output, time, state] = noStateTimeVaryingSim(this_lft, input)
%% Simulate a Memoryless LTV system

% Obtain time_indices which map this_lft properties to input time
total_time = size(input, 2);
horizon_period = [total_time - 1, 1];
horizon_period = commonHorizonPeriod([this_lft.horizon_period; horizon_period]);
time_indices = makeNewIndices(this_lft.horizon_period, horizon_period);
time_indices = time_indices(1 : total_time);

% Time default defined because this is memoryless
time = 0 : total_time - 1;

% Convert input to cell to allow for uniform treatment
if ~iscell(input)
    input = mat2cell(input, size(input, 1), ones(1, total_time));
end
output = cell(1, total_time);
d = this_lft.d(time_indices);
% Though time starts at zero, matlab doesn't allow zero-indexed arrays, so we start at i = 1
for i = 1:total_time
    output{1, i} = d{i} * input{1, i};
end

% Bring output to matrices if not needed as cell arrays
dim_out = size(this_lft, 1);
if all(dim_out(1) == dim_out)
    output = cell2mat(output);
end

state = double.empty(0, total_time);
end

%%  CHANGELOG
% Oct. 6, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)