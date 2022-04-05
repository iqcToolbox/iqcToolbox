function lft_out = toLft(varargin)
%% TOLFT for converting double, Delta, and ss to Ulft objects.
%
%     lft_out = toLft(3*eye(n))
%     lft_out = toLft(DeltaSlti('d'))
%     lft_out = toLft(ss(a, b, c, d))
%     lft_out = toLft(ss(a, b, c, d, timestep))
%     lft_out = toLft(a, b, c, d) (where a, b, c, d, may be a double or a Ulft)
%     lft_out = toLft(a, b, c, d, timestep) (where a,b,c,d may be a double or Ulft)
%
%     Variables:
%     ---------
%       Input:
%         varargin : May be a number of different types and arguments.
%             - If giving a single argument, varargin may be a:
%                  double, cell of doubles, Delta object, continuous-time ss, or discrete-time ss
%             - If giving four arguments, each argument must be a:
%                  double or Ulft object, which represent state-space matrices (with uncertainties if a Ulft object) describing the dynamical system
%             - If giving five arguments,
%                  the first four arguments must be double or Ulft objects (representing discrete-time state-space matrices)
%                  the fifth argument must be a double (representing the sampling time of the discrete-time system)
%                      if the sampling time is unknown, provide -1 or []
%       Output:
%         lft_out : Ulft object :: the resultant lft
%
%     See also Ulft, Ulft.Ulft.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Calling appropriate functions for different input arguments
if nargin == 1
    v1 = varargin{1};
    if isa(v1, 'numeric')
        lft_out = numericToLft(v1);
    elseif isa(v1, 'ss')
        lft_out = ssToLft(v1);
    elseif isa(v1, 'Delta')
        lft_out = deltaToLft(v1);
    elseif isa(v1, 'Ulft')
        lft_out = v1;
    else
        error('toLft:toLft', ['When providing a single argument for toLft, ',...
              'that input must be a numeric, ss, or Delta']);
    end
elseif nargin == 2
    lft_out = numericCellToLft(varargin{:});
elseif nargin == 4 || nargin == 5
    lft_out = abcdToLft(varargin{:});
else
    error('toLft:toLft', 'Must provide 1, 2, 4, or 5 arguments for toLft')
end

end

function lft_out = numericToLft(num)
%% NUMERICTOLFT for converting numeric objects (double, int, etc.) to Ulft
validateattributes(num, {'numeric'}, {'nonempty', 'finite'})

[dim_out, dim_in] = size(num);
lft_out = Ulft([], zeros(0, dim_in), zeros(dim_out, 0), num, SequenceDelta());
end

function lft_out = numericCellToLft(sequence, horizon_period)
%% NUMERICCELLTOLFT for converting sequences of numeric objects to Ulft
validateattributes(horizon_period, {'numeric'}, {'size', [1,2],...
                                               'integer',...
                                               'nonnegative'});
validateattributes(horizon_period(2), {'numeric'}, {'positive'})  
total_time = sum(horizon_period);
validateattributes(sequence, {'cell'}, {'numel', total_time})
for i = 1:total_time
    validateattributes(sequence{i}, {'numeric'}, {'nonempty', 'finite', 'nonnan'})
end

dim_out = cellfun(@(c) size(c, 1), sequence);
dim_in  = cellfun(@(c) size(c, 2), sequence);

a = cell(1, total_time);
b = cell(1, total_time);
c = cell(1, total_time);
d = cell(1, total_time);

for i = 1:length(sequence)
    a{i} = [];
    b{i} = zeros(0, dim_in(i));
    c{i} = zeros(dim_out(i), 0);
    d{i} = sequence{i};
end

lft_out = Ulft(a, b, c, d, SequenceDelta(), 'horizon_period', horizon_period);
end

function lft_out = ssToLft(g)
%% SSTOLFT for converting ss objects to Ulft
validateattributes(g, {'ss'}, {'nonempty'})

if g.Ts == 0
    lft_out = abcdToLft(g.a, g.b, g.c, g.d);
else
    lft_out = abcdToLft(g.a, g.b, g.c, g.d, g.Ts);
end
end

function lft_out = deltaToLft(del)
%% DELTATOLFT for converting Delta objects to Ulft
validateattributes(del, {'Delta'}, {'nonempty'})

dim_in  = del.dim_in;
dim_out = del.dim_out;
total_time = sum(del.horizon_period);

a = cell(1, total_time);
b = cell(1, total_time);
c = cell(1, total_time);
d = cell(1, total_time);
for i = 1:total_time
    a{i} = zeros(dim_in(i), dim_out(i));
    b{i} = eye(dim_in(i));
    c{i} = eye(dim_out(i));
    d{i} = zeros(dim_out(i), dim_in(i));
end
lft_out = Ulft(a, b, c, d, del, 'horizon_period', del.horizon_period);
end

function lft_out = abcdToLft(varargin)
%% ABCDTOLFT for converting uncertain SS (represented with a, b, c, d 
% matrices) to Ulft

% Ensuring first four arguments are Ulft objects with admissible uncertainties
goodDelta = @(del) isa(del, 'DeltaSlti') || ...
                   isa(del, 'DeltaSltv') || ...
                   isa(del, 'DeltaSltvRateBnd') || ...
                   isa(del, 'DeltaSltvRepeated') || ...
                   isa(del, 'DeltaSectorBounded');
for i = 1:4
    varargin{i} = toLft(varargin{i});
    validateattributes(varargin{i}, {'Ulft'}, {'nonempty'})
    for j = 1:length(varargin{i}.delta.deltas)
        assert(goodDelta(varargin{i}.delta.deltas{j}),...
               'toLft:toLft',...
               'Cannot make a ss-based lft where matrices have dynamic operators')
    end
end
% a, b, c, d are Ulft objects
a = varargin{1};
b = varargin{2};
c = varargin{3};
d = varargin{4};

if nargin == 5
    timestep = varargin{5};
    assert(~isempty(timestep) && (all(timestep == -1) || all(timestep > 0)),...
           'toLft',...
           'timestep for discrete LFTs must be -1 (unspecified) or positive')
else
    timestep = 0;
end

% Check that continuous-time systems have LTI matrices
if timestep == 0
    assert(all(a.horizon_period == [0, 1]),...
           'toLft:abcdToLft',...
           'Cannot use time-varying a-matrix for continuous-time systems')
    assert(all(b.horizon_period == [0, 1]),...
           'toLft:abcdToLft',...
           'Cannot use time-varying b-matrix for continuous-time systems')
    assert(all(c.horizon_period == [0, 1]),...
           'toLft:abcdToLft',...
           'Cannot use time-varying c-matrix for continuous-time systems')
    assert(all(d.horizon_period == [0, 1]),...
           'toLft:abcdToLft',...
           'Cannot use time-varying d-matrix for continuous-time systems')
end

% Construct matrices for combined, dynamic lft
dim_states = cellfun(@(a) size(a, 2), a.d);
abcd = [a, b;
        c, d];
total_time = sum(abcd.horizon_period);
m11 = cell(1, total_time);
m12 = cell(1, total_time);
m21 = cell(1, total_time);
m22 = cell(1, total_time);
for i = 1:total_time
    a11 = abcd.d{i}(1:dim_states, 1:dim_states);
    a12 = abcd.c{i}(1:dim_states, :);
    b1  = abcd.d{i}(1:dim_states, dim_states+1:end);
    a21 = abcd.b{i}(:, 1:dim_states);
    a22 = abcd.a{i};
    b2  = abcd.b{i}(:, dim_states+1:end);
    c1  = abcd.d{i}(dim_states+1:end, 1:dim_states);
    c2  = abcd.c{i}(dim_states+1:end, :);
    d   = abcd.d{i}(dim_states+1:end, dim_states+1:end);
    m11{i} = [a11, a12; a21, a22];
    m12{i} = [b1; b2];
    m21{i} = [c1, c2];
    m22{i} = d;
end

% Construct combined delta and lft
if timestep == 0
    delta_state = DeltaIntegrator(dim_states);
else
    delta_state = DeltaDelayZ(dim_states, timestep, abcd.horizon_period);    
end
delta = SequenceDelta(delta_state, abcd.delta.deltas);
lft_out = Ulft(m11, m12, m21, m22, delta, 'horizon_period', abcd.horizon_period);
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Changed behavior when converting from cell arrays so that horizon_period is specified - Micah Fry
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)