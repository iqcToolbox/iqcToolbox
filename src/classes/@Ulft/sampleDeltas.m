function this_lft = sampleDeltas(this_lft, deltas, values, varargin)
%% SAMPLEDELTAS method for realizing some or all deltas to specified or randomly sampled values.
%
%    Example Usage:
%      this_lft = this_lft.sampleDeltas({'somedelta', 'otherdelta'}, {{}, otherlft}, 'override', true, 'time_invariant', true);
%
%    Variables:
%    ---------
%      Inputs:
%        deltas : cell array of strings :: names of the deltas to sample
%        values : cell array of Ulft objects or empties :: corresponding samples for each delta; empty cells request random generation
%        override : boolean :: if true, then only the dimensions and horizon_period of each value will be validated (DEFAULT: false)
%        time_invariant : boolean :: if true, then all randomly generated delta samples will be time-invariant (DEFAULT: false)
%
%      Output:
%        this_lft : Ulft object :: the LFT with the deltas sampled
%
%    See also removeUncertainty, random.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

% Validate arguments (except for delta-specific validation which happens later)
parser = inputParser;
parser.addRequired('deltas',...
                   @(arg) iscell(arg) && all(cellfun(@(a) ischar(a), arg)));
isCellOfUlft = @(arg) iscell(arg) &&...
                     all(cellfun(@(a) isa(a, 'Ulft') || isempty(a), arg));
parser.addRequired('values', isCellOfUlft);
parser.addParameter('override',...
                    false,...
                    @(o) validateattributes(o, 'logical', {'scalar'}));
parser.addParameter('time_invariant',...
                    false,...
                    @(ti) validateattributes(ti, 'logical', {'scalar'}));
parser.parse(deltas, values, varargin{:});
assert(length(deltas) == length(values),...
       'Ulft:sampleDeltas',...
       ['Must provide exactly one value for each specified delta. ',...
        'Use {} to request a random value.']);

% Extract optional arguments
override = parser.Results.override;
time_invariant = parser.Results.time_invariant;

% Determine if this LFT is continuous-time, discrete-time, or memoryless
timestep = this_lft.timestep;

% Create diagonal LFT of (all or partially) realized deltas
realized = cell(1, length(this_lft.delta.deltas));
for i = 1 : length(this_lft.delta.deltas)
    % Check to see if this delta was named for sampling
    delta_obj = this_lft.delta.deltas{i};
    value_ind = find(strcmp(delta_obj.name, deltas));
    % Keep as-is if not
    if isempty(value_ind)
        realized{i} = toLft(delta_obj);
    else
        % Otherwise, grab the corresponding sample value
        value = values{value_ind(1)};
        % Empty value requests random sample
        if isempty(value)
            value = delta_obj.sample(timestep);
            % Reduce to time-invariant if flagged to do so
            if time_invariant
                if value.timestep
                    value.delta.deltas{1} = ...
                        DeltaDelayZ(value.delta.deltas{1}.dim_out(1),...
                                    value.delta.deltas{1}.timestep(1));
                end
                value = Ulft(value.a{1}, value.b{1}, value.c{1}, value.d{1},...
                             value.delta,...
                             'horizon_period', [0,1]);
                value = value.matchHorizonPeriod(this_lft.horizon_period);
            end
        else
            % Just check dimensions...
            if override
                assert(all(value.horizon_period == delta_obj.horizon_period),...
                       'Ulft:sampleDeltas',...
                       ['Specified value has invalid horizon_period for delta "',...
                       delta_obj.name,'".']);
                assert(isequal(size(value, 1), delta_obj.dim_out) &&...
                       isequal(size(value, 2), delta_obj.dim_in),...
                       'Ulft:sampleDeltas',...
                       ['Specified value has invalid dimensions for delta "',...
                        delta_obj.name,'".']);
            % ...or check everything
            else
                delta_obj.validateSample(value, timestep);
            end
        end
        % Sampled or provided value is saved in the array of Ulft realizations
        realized{i} = value;
    end
end
realized = blkdiag(realized{:});

% Forge the original ABCD into an LFT with just D
total_time = sum(this_lft.horizon_period);
base_a = cell(1, total_time);
base_b = cell(1, total_time);
base_c = cell(1, total_time);
base_d = cell(1, total_time);
for t = 1 : total_time
    base_d{t} = [this_lft.a{t}, this_lft.b{t}; this_lft.c{t}, this_lft.d{t}];
    [dim_out, dim_in] = size(base_d{t});
    base_c{t} = zeros(dim_out, 0);
    base_b{t} = zeros(0, dim_in);
    base_a{t} = [];
end
base = Ulft(base_a, base_b, base_c, base_d, SequenceDelta(),...
            'horizon_period', this_lft.horizon_period,...
            'disturbance', this_lft.disturbance,...
            'performance', this_lft.performance);

% Interconnect the "deltas LFT" with the "ABCD LFT" to close the loop
this_lft = interconnect(realized, base);

end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Added after v0.5.0 - Jason Nezvadovitz (jason.nezvadovitz@ll.mit.edu)