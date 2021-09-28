function lft_out = vertcat(varargin)
%% VERTCAT overloaded function for vertically concatenating Ulft objects.
%
%     lft_out = vertcat(lft1, lft2)
%     lft_out = [lft1; lft2]
%
%     Variables:
%     ---------
%       Input:
%         varargin : Ulft or Ulft-convertible object(s) :: one or many lft(s) to be concatenated
%       Output:
%         lft_out : Ulft object :: Concatenated lft
%
%     See also Ulft, horzcat, blkdiag.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Check correctness and consistency of inputs
for i = 1 : nargin
    if ~isa(varargin{i}, 'Ulft')
        varargin{i} = toLft(varargin{i}); 
    end
    validateattributes(varargin{i}, 'Ulft', {'nonempty'}, mfilename)
end
horizon_periods = cell2mat(cellfun(@(lft) lft.horizon_period,...
                                   varargin',...
                                   'UniformOutput', false));
horizon_period = commonHorizonPeriod(horizon_periods);

varargin{1} = matchHorizonPeriod(varargin{1}, horizon_period);
[~, dim_in_lft] = size(varargin{1});
names_disturbance = varargin{1}.disturbance.names;

for i = 1 : nargin
    if (~all(varargin{i}.horizon_period == horizon_period))
        varargin{i} = matchHorizonPeriod(varargin{i}, horizon_period);
    end
    assert(all(varargin{i}.horizon_period == horizon_period),...
           'Ulft:vertcat',...
           'Ulft objects must have the same horizon_period')
    dim_in = size(varargin{i}, 2);       
    assert(all(dim_in == dim_in_lft, 'all'),...
           'Ulft:vertcat',...
           'Ulft objects must have the same number of input channels')
    assert(all(strcmp(varargin{i}.disturbance.names, names_disturbance)),...
           'Ulft:vertcat',...
           'Ulft objects must have the same disturbances')
end

%% Form delta, performance, a, b, c, and d matrices
deltas = cellfun(@(lft) lft.delta.deltas, varargin, 'UniformOutput', false);
delta = SequenceDelta(deltas{:});

% Add pertinent output size to each performance.chan_out
total_time = sum(horizon_period);
lft_performance = cell(1, nargin);
channel_shift = num2cell(zeros(1, total_time));
for i = 1:nargin
    performances = varargin{i}.performance.performances;
    for j = 1:length(performances)
        performances{j}.chan_out = cellfun(@plus,...
                                           performances{j}.chan_out,...
                                           channel_shift,...
                                           'UniformOutput', false);
    end
    channel_shift = cellfun(@plus,...
                            channel_shift,...
                            num2cell(size(varargin{i}, 1)),...
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
    
    d{i} = vertcat(d_lfts{:});
    c{i} = blkdiag(c_lfts{:});
    b{i} = vertcat(b_lfts{:});
    a{i} = blkdiag(a_lfts{:});
end

%% Construct new lft
lft_out = Ulft(a, b, c, d, delta,...
           'horizon_period', horizon_period,...
           'performance', performance,...
           'disturbance', varargin{1}.disturbance);

end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)