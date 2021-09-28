function this_lft = removeDisturbance(this_lft, disturbance)
%% REMOVEDISTURBANCE method for removing specified disturbances from the lft.
%
%     lft_out = removeDisturbance(this_lft, 2)
%     lft_out = removeDisturbance(this_lft, 'kick')
%     lft_out = removeDisturbance(this_lft, [2, 6, 3])
%     lft_out = removeDisturbance(this_lft, {'kick', 'punch'})
%
%     Variables:
%     ---------
%       Input:
%         this_lft : Ulft object :: the lft whose specified disturbances are to be removed
%         disturbance : string/cell of strings/array of natural numbers :: the names of the disturbance objects to be removed
%       Output:
%         this_lft : Ulft object :: the lft whose specified disturbances are removed
%
%     See also Ulft, removeUncertainty, removePerformance.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

% Check and process inputs
if isempty(disturbance); return; end
validateattributes(disturbance, {'numeric', 'char', 'cell'}, {'nonempty'})

% Convert all input types to be indices
isChar = @(dist) isa(dist, 'char');
if isChar(disturbance)
    disturbance_ind = find(strcmp(disturbance, this_lft.disturbance.names), 1);
    if isempty(disturbance_ind)
        warning('Ulft:removeDisturbance',...
                ['The named disturbance, ',...
                 disturbance,...
                 ', does not appear in the lft'])
    end
elseif isa(disturbance, 'cell') && all(cellfun(isChar, disturbance), 'all')
    disturbance_ind = [];
    for i = 1:length(disturbance)
        disturbance_i = disturbance{i};
        disturbance_ind = [disturbance_ind, find(strcmp(disturbance_i, this_lft.disturbance.names))];
        if isempty(find(strcmp(disturbance_i, this_lft.disturbance.names), 1))
            warning('Ulft:removeDisturbance',...
                    ['The named disturbance, ',...
                     disturbance_i,...
                     ', does not appear in the lft'])
        end
    end
elseif isnumeric(disturbance)
    validateattributes(disturbance, 'numeric', {'positive',...
                                                'integer',...
                                                'finite',...
                                                'nonnan'})
        bad_inds = disturbance(disturbance > length(this_lft.disturbance.names));
        if ~isempty(bad_inds)
            warning('Ulft:removeDisturbance',...
                    ['The numbered disturbance, ',...
                     num2str(bad_inds),...
                     ', does not appear in the lft'])
        end
    disturbance_ind = disturbance;
end

% Define the indices of disturbances to keep
keep_ind = setdiff(1:length(this_lft.disturbance.names), disturbance_ind);

% Define reduced list of disturbances
new_disturbance = this_lft.disturbance.disturbances(keep_ind);

% Make reduced lft
if isempty(new_disturbance)
    this_lft = Ulft(this_lft.a, this_lft.b, this_lft.c, this_lft.d, this_lft.delta,...
                    'horizon_period', this_lft.horizon_period,...
                    'performance', this_lft.performance);
else
    this_lft = Ulft(this_lft.a, this_lft.b, this_lft.c, this_lft.d, this_lft.delta,...
                    'horizon_period', this_lft.horizon_period,...
                    'disturbance', new_disturbance,...
                    'performance', this_lft.performance);
end

end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Jason Nezvadovitz (jason.nezvadovitz@ll.mit.edu)