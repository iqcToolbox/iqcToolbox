function lft_out = blkdiag(varargin)
%% BLKDIAG overloaded function for block-diagonally concatenating Ulft objects.
%
%     lft_out = blkdiag(lft1, lft2, lft3, ...)
%
%     Variables:
%     ---------
%       Input:
%         varargin : Ulft or Ulft-convertible object(s) :: one or many lft(s) to be concatenated
%       Output:
%         lft_out : Ulft object :: Concatenated lft
%
%     See also vertcat, horzcat.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Check correctness and consistency of inputs
for i = 1 : nargin
    if ~isa(varargin{i}, 'Ulft')
        varargin{i} = toLft(varargin{i});
    end
    validateattributes(varargin{i}, {'Ulft'}, {'nonempty'}, mfilename)
end
horizon_periods = cell2mat(cellfun(@(lft) lft.horizon_period,...
                                   varargin',...
                                   'UniformOutput', false));
horizon_period = commonHorizonPeriod(horizon_periods);

varargin{1} = matchHorizonPeriod(varargin{1}, horizon_period);
total_time = sum(horizon_period);
for i = 1 : nargin
    if (~all(varargin{i}.horizon_period == horizon_period))
        varargin{i} = matchHorizonPeriod(varargin{i}, horizon_period);
    end
    assert(all(varargin{i}.horizon_period == horizon_period),...
           'Ulft:diag',...
           'Ulft objects must have the same horizon_period')
end

%% Form delta, performance, a, b, c, and d matrices
deltas = cellfun(@(lft) lft.delta.deltas, varargin, 'UniformOutput', false);
delta = SequenceDelta(deltas{:});

% Add pertinent input size to each disturbance.chan_in
lft_disturbance = cell(1, nargin);
dist_chan_shift = num2cell(zeros(1, total_time));
for i = 1:nargin
    disturbances = varargin{i}.disturbance.disturbances;
    for j = 1:length(disturbances)
        disturbances{j}.chan_in = cellfun(@plus,...
                                          disturbances{j}.chan_in,...
                                          dist_chan_shift,...
                                          'UniformOutput', false);
    end
    dist_chan_shift = cellfun(@plus,...
                            dist_chan_shift,...
                            num2cell(size(varargin{i}, 2)),...
                            'UniformOutput', false);
    lft_disturbance{i} = disturbances;
end
disturbance = SequenceDisturbance(lft_disturbance{:});

% Add pertinent input/output size to each performance.chan_out/chan_in
lft_performance = cell(1, nargin);
perf_chan_out_shift = num2cell(zeros(1, total_time));
perf_chan_in_shift = num2cell(zeros(1, total_time));
for i = 1:nargin
    performances = varargin{i}.performance.performances;
    for j = 1:length(performances)
        performances{j}.chan_out = cellfun(@plus,...
                                           performances{j}.chan_out,...
                                           perf_chan_out_shift,...
                                           'UniformOutput', false);
        performances{j}.chan_in = cellfun(@plus,...
                                          performances{j}.chan_in,...
                                          perf_chan_in_shift,...
                                          'UniformOutput', false);
    end
    perf_chan_out_shift = cellfun(@plus,...
                                  perf_chan_out_shift,...
                                  num2cell(size(varargin{i}, 1)),...
                                  'UniformOutput', false);
    perf_chan_in_shift = cellfun(@plus,...
                                 perf_chan_in_shift,...
                                 num2cell(size(varargin{i}, 2)),...
                                 'UniformOutput', false);
    lft_performance{i} = performances;
end
performance = SequencePerformance(lft_performance{:});

% a, b, c, d matrices
d = cell(1, total_time);
c = cell(1, total_time);
b = cell(1, total_time);
a = cell(1, total_time);
for i = 1:total_time
    d_lfts = cellfun(@(lft) lft.d{i}, varargin, 'UniformOutput', false);
    c_lfts = cellfun(@(lft) lft.c{i}, varargin, 'UniformOutput', false);
    b_lfts = cellfun(@(lft) lft.b{i}, varargin, 'UniformOutput', false);
    a_lfts = cellfun(@(lft) lft.a{i}, varargin, 'UniformOutput', false);

    d{i} = blkdiag(d_lfts{:});
    c{i} = blkdiag(c_lfts{:});
    b{i} = blkdiag(b_lfts{:});
    a{i} = blkdiag(a_lfts{:});
end

%% Construct new lft
lft_out = Ulft(a, b, c, d, delta,...
           'horizon_period', horizon_period,...
           'performance', performance,...
           'disturbance', disturbance);
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Jason Nezvadovitz (jason.nezvadovitz@ll.mit.edu)