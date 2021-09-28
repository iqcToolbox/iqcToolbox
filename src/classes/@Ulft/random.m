function this_lft = random(varargin)
%% RANDOM method (static) for constructing a random / arbitrary Ulft object.
%
%    Usages:
%      lft_out = Ulft.random()
%      lft_out = Ulft.random('dim_in', dim_in, 'dim_out', dim_out, 'num_deltas', num_deltas, 'req_deltas', req_deltas, 'horizon_period', horizon_period)
%
%    Inputs:
%      dim_in : positive integer :: specifically desired dimensionality of the LFT's input
%      dim_out : positive integer :: specifically desired dimensionality of the LFT's output
%      num_deltas : nonnegative integer :: specifically desired total number of deltas in the LFT's delta sequence
%      req_deltas : cell array of Delta objects or strings (Delta type names) :: specifically desired deltas to be included in the LFT's delta sequence
%      horizon_period : 1 x 2 array of naturals :: specifically desired [horizon, period] (in timesteps) of LFT a,b,c,d matrices
%
%    Outputs:
%      lft_out : Ulft object :: the randomly generated lft
%
%    Notes:
%      All unspecified inputs are randomized.
%      The deltas in 'req_deltas' count toward 'num_deltas', but are guaranteed to appear even if `length(req_deltas) > num_deltas`.
%      To ensure that the LFT defines a continuous / discrete dynamic, include a DeltaIntegrator / DeltaDelayZ in 'req_deltas'.
%      If the LFT defines a continuous / discrete dynamic, the state-subblock of the A-matrix will be Hurwitz / Schur for all time.
%      If the LFT defines a continuous dynamic, the horizon period will be [0,1] regardless of the input argument 'horizon_period'.
%
%    See also Ulft, toLft.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% General ranges for randomization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (for randomly naming deltas)
ALPHABET = ['A':'Z', 'a':'z'];
% (for scaling random numbers)
N = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parser = inputParser();
parser.addParameter('dim_in',...
                    randi([1,N]),...
                    @(arg) validateattributes(arg,...
                                              'numeric',...
                                              {'integer', 'positive'}));
parser.addParameter('dim_out',...
                    randi([1,N]),...
                    @(arg) validateattributes(arg,...
                                              'numeric',...
                                              {'integer', 'positive'}));
parser.addParameter('num_deltas',...
                    randi([0,N]),...
                    @(arg) validateattributes(arg,...
                                              'numeric',...
                                              {'integer', 'nonnegative'}));
cellOfDelsOrChars = @(arr) all(cellfun(@(a) isa(a,'Delta') || ischar(a), arr));
parser.addParameter('req_deltas',...
                    {},...
                    @(arg) iscell(arg) && cellOfDelsOrChars(arg));
parser.addParameter('horizon_period',...
                    [],...
                    @(arg) validateattributes(arg,...
                                              'numeric',...
                                              {'integer', 'nonnegative',...
                                               'size', [1,2]}));

parser.parse(varargin{:});
dim_in = parser.Results.dim_in;
dim_out = parser.Results.dim_out;
num_deltas = max(parser.Results.num_deltas, length(parser.Results.req_deltas));
req_deltas = parser.Results.req_deltas;
horizon_period = parser.Results.horizon_period;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine temporal mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The horizon-period can be implied by the user through req_deltas
if isempty(horizon_period) && ~isempty(req_deltas)
    req_horizon_periods = repmat([0, 1], length(req_deltas), 1);
    found = false;
    for i = 1 : length(req_deltas)
        if isa(req_deltas{i}, 'Delta')
            req_horizon_periods(i,:) = req_deltas{i}.horizon_period;
            found = true;
        end
    end
    if found
        horizon_period = commonHorizonPeriod(req_horizon_periods);
    end
end

tmode = 'unspecified';
can_be_continuous = isempty(horizon_period) || isequal(horizon_period, [0,1]);

% The user can specify continuous/discrete by including an integrator or delay in req_deltas
for i = 1 : length(req_deltas)
    if isa(req_deltas{i}, 'DeltaIntegrator') ...
       || strcmp(req_deltas{i}, 'DeltaIntegrator')
        if ~can_be_continuous
            error('Continuous-time LFTs must have a horizon_period of [0,1]');
        end
        tmode = 'continuous';
        break;
    elseif isa(req_deltas{i}, 'DeltaDelayZ') ...
           || strcmp(req_deltas{i}, 'DeltaDelayZ')
        tmode = 'discrete';
        break;
    end
end

% If they didn't specify, then we have to randomly decide whether to make it dynamical
if strcmp(tmode, 'unspecified')
    if num_deltas == length(req_deltas)
        % (no room to randomly include an integrator or delay)
        tmode = 'stateless';
    else
        tmodes = {'continuous', 'discrete', 'stateless'};
        if can_be_continuous
            tmode = tmodes{randi([1, 3])};
        else
            tmode = tmodes{randi([2, 3])};
        end
        switch tmode
            case 'continuous'
                req_deltas{end+1} = random_delta_integrator(N);
            case 'discrete'
                req_deltas{end+1} = random_delta_delayz(N);
        end
    end
end

% If horizon_period wasn't implied by req_deltas, then this one will be used
if strcmp(tmode, 'continuous')
    horizon_period = [0,1];
elseif isempty(horizon_period)
    horizon_period = random_horizon_period(N, 0.3);
end
total_time = sum(horizon_period);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate deltas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deltas = cell(1, num_deltas);

% The first few deltas must be the required ones (randomly generated if only the type was specified)
for i = 1 : length(req_deltas)
    if isa(req_deltas{i}, 'Delta')
        deltas{i} = req_deltas{i};
    else
        name = ['req_', ALPHABET(randi(length(ALPHABET), [1,N]))];
        switch req_deltas{i}
            case 'DeltaIntegrator'
                deltas{i} = random_delta_integrator(N);
            case 'DeltaDelayZ'
                deltas{i} = random_delta_delayz(N);
            case 'DeltaBounded'
                deltas{i} = random_delta_bounded(name, N, horizon_period);
            case 'DeltaDlti'
                deltas{i} = random_delta_dlti(name, N);
            case 'DeltaSlti'
                deltas{i} = random_delta_slti(name, N);
            case 'DeltaSltv'
                deltas{i} = random_delta_sltv(name, N, horizon_period);
            case 'DeltaSltvRateBnd'
                deltas{i} = random_delta_sltvratebnd(name, N, horizon_period);
            otherwise
                error(['Required Delta-type ''', req_deltas{i},...
                       ''' has no randomizer.']);
        end
    end
end

% The rest are randomly chosen
for i = length(req_deltas)+1 : num_deltas
    name = ALPHABET(randi(length(ALPHABET), [1,N]));
    switch randi(5)
        case 1
            deltas{i} = random_delta_bounded(name, N, horizon_period);
        case 2
            deltas{i} = random_delta_dlti(name, N);
        case 3
            deltas{i} = random_delta_slti(name, N);
        case 4
            deltas{i} = random_delta_sltv(name, N, horizon_period);
        case 5
            deltas{i} = random_delta_sltvratebnd(name, N, horizon_period);
        otherwise
            error('Switch-case mismatch (this is a backend bug)');
    end
end

% Now that we have our final delta list, count up the important dimensions
dim_deltas_in = 0;
dim_deltas_out = 0;
dim_state = 0;
for i = 1 : num_deltas
    dim_deltas_in = dim_deltas_in + deltas{i}.dim_in(1);
    dim_deltas_out = dim_deltas_out + deltas{i}.dim_out(1);
    if isa(deltas{i}, 'DeltaIntegrator') || isa(deltas{i}, 'DeltaDelayZ')
        dim_state = deltas{i}.dim_out(1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate LFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare a stable dynamic via the state-subblock of the A-matrix
if ~strcmp(tmode, 'stateless')
    [basis, poles] = eig(N*(2*rand(dim_state, dim_state)-1));
    poles = diag(poles);
    switch tmode
        % (hurwitz)
        case 'continuous'
            unstable = (real(poles) > 0);
            poles(unstable) = -poles(unstable);
        % (schur)
        case 'discrete'
            mags = abs(poles);
            unstable = (mags > 0.98);
            poles(unstable) = ...
                poles(unstable) .* tanh(mags(unstable)/10) ./ mags(unstable);
    end
    stable = real(basis * diag(poles) * inv(basis));
end

% Form abcd-matrices
a = cell(1, total_time);
b = cell(1, total_time);
c = cell(1, total_time);
d = cell(1, total_time);
for t = 1 : total_time
    a{t} = N*(2*rand(dim_deltas_in, dim_deltas_out)-1);
    b{t} = N*(2*rand(dim_deltas_in, dim_in        )-1);
    c{t} = N*(2*rand(dim_out,       dim_deltas_out)-1);
    d{t} = N*(2*rand(dim_out,       dim_in        )-1);
    % Stabilize the dynamic
    if ~strcmp(tmode, 'stateless')
        % (assumes the integrator / delay is always first in Ulft.delta)
        % (state-subblock held constant to guarantee LTV stability)
        a{t}(1:dim_state, 1:dim_state) = stable;
    end
end

this_lft = Ulft(a, b, c, d, deltas, 'horizon_period', horizon_period);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Randomization Helpers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function horizon_period = random_horizon_period(N, p)
    if rand() < p
        horizon_period = [0,1];
    else
        horizon_period = [randi([0,N]), randi([1,N])];
    end
end

function delta = random_delta_integrator(N)
    delta = DeltaIntegrator(randi([1,N]));
end

function delta = random_delta_delayz(N)
    delta = DeltaDelayZ(randi([1,N]), rand());
end

function delta = random_delta_bounded(name, N, horizon_period)
    delta = DeltaBounded(name, randi([1,N]), randi([1,N]), N*rand(), horizon_period);
end

function delta = random_delta_dlti(name, N)
    delta = DeltaDlti(name, randi([1,N]), randi([1,N]), N*rand());
end

function delta = random_delta_slti(name, N)
    bound = N*rand();
    delta = DeltaSlti(name, randi([1,N]), -bound, bound);
end

function delta = random_delta_sltv(name, N, horizon_period)
    bound = N*rand();
    delta = DeltaSltv(name, randi([1,N]), -bound, bound, horizon_period);
end

function delta = random_delta_sltvratebnd(name, N, horizon_period)
    bound = N*rand();
    rate_bound = N*rand();
    delta = DeltaSltvRateBnd(name, randi([1,N]), -bound, bound,...
        -rate_bound, rate_bound, horizon_period, randi([2,N]));
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Jason Nezvadovitz (jason.nezvadovitz@ll.mit.edu)