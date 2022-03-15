function this_lft = removePerformance(this_lft, performance)
%% REMOVEPERFORMANCE method for removing specified performances from the lft.
%
%     lft_out = removePerformance(this_lft, 2)
%     lft_out = removePerformance(this_lft, 'fuel_use')
%     lft_out = removePerformance(this_lft, [2, 6, 3])
%     lft_out = removePerformance(this_lft, {'fuel_use', 'coolness'})
%
%     Variables:
%     ---------
%       Input:
%         this_lft : Ulft object :: the lft whose specified performances are to be removed
%         performance : string/cell of strings/array of natural numbers :: the names of the performance objects to be removed
%       Output:
%         this_lft : Ulft object :: the lft whose specified performances are removed
%
%     See also Ulft, removeUncertainty, removeDisturbance.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

% Check and process inputs
if isempty(performance); return; end
validateattributes(performance, {'numeric', 'char', 'cell'}, {'nonempty'})

% Convert all input types to be indices
isChar = @(perf) isa(perf, 'char');
if isChar(performance)
    performance_ind = find(strcmp(performance, this_lft.performance.names), 1);
    if isempty(performance_ind)
        warning('Ulft:removePerformance',...
                ['The named performance, ',...
                 performance,...
                 ', does not appear in the lft'])
    end
elseif isa(performance, 'cell') && all(all(cellfun(isChar, performance)))
    performance_ind = [];
    for i = 1:length(performance)
        performance_i = performance{i};
        performance_ind = [performance_ind, find(strcmp(performance_i, this_lft.performance.names))];
        if isempty(find(strcmp(performance_i, this_lft.performance.names), 1))
            warning('Ulft:removePerformance',...
                    ['The named performance, ',...
                     performance_i,...
                     ', does not appear in the lft'])
        end
    end
elseif isnumeric(performance)
    validateattributes(performance, {'numeric'}, {'positive',...
                                                'integer',...
                                                'finite',...
                                                'nonnan'})
        bad_inds = performance(performance > length(this_lft.performance.names));
        if ~isempty(bad_inds)
            warning('Ulft:removePerformance',...
                    ['The numbered performance, ',...
                     num2str(bad_inds),...
                     ', does not appear in the lft'])
        end
    performance_ind = performance;
end

% Define the indices of performances to keep
keep_ind = setdiff(1:length(this_lft.performance.names), performance_ind);

% Define reduced list of performances
new_performance = this_lft.performance.performances(keep_ind);

% Make reduced lft
if isempty(new_performance)
    this_lft = Ulft(this_lft.a, this_lft.b, this_lft.c, this_lft.d, this_lft.delta,...
                    'horizon_period', this_lft.horizon_period,...
                    'disturbance', this_lft.disturbance);
else
    this_lft = Ulft(this_lft.a, this_lft.b, this_lft.c, this_lft.d, this_lft.delta,...
                    'horizon_period', this_lft.horizon_period,...
                    'disturbance', this_lft.disturbance,...
                    'performance', new_performance);
end

end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Jason Nezvadovitz (jason.nezvadovitz@ll.mit.edu)