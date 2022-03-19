function this_lft = gatherLft(this_lft)
%% GATHERLFT method for gathering multiples of the same Delta/Disturbance/Performance block 
%     to be represented as a single Delta/Disturbance/Performance for each set of multiples.
%
%     gathered_lft = gatherLft(this_lft)
%
%     Variables:
%     ---------
%       Input:
%         this_lft : Ulft object :: the lft whose Deltas/Disturbances/Performances are to be gathered
%       Output:
%         this_lft : Ulft object :: the lft with gathered sets of Deltas/Disturbance/Performances
%
%     See also Ulft, reorderLftDelta

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Gathering Delta
isState = @(del) isa(del, 'DeltaDelayZ') || isa(del, 'DeltaIntegrator');
state_deltas = cellfun(isState, this_lft.delta.deltas);
state_not_in_front = any(state_deltas) && ~state_deltas(1);

number_of_deltas = length(this_lft.delta.names);
repeated_deltas = (number_of_deltas ~= length(unique(this_lft.delta.names)));
if (state_not_in_front || repeated_deltas)
    old_order = 1:number_of_deltas;
    
    %%%%%%%%%%%%%%%%%
    % Define ordering of non-state deltas to put duplicates next to each other
    %%%%%%%%%%%%%%%%%
    non_state_inds = old_order(~state_deltas);
    unsorted_names = this_lft.delta.names(non_state_inds);
    del_ns_names = this_lft.delta.names(non_state_inds);
    % Initialize re-sorted order and unsorted list (start w/ first element)
    gather_order = nan(1, length(unsorted_names));
    gather_order(1) = 1;
    unsorted_names{1} = [];
    for i = 2:length(unsorted_names)
        % Check name of last uncertainty in re-sorted order
        last_name = del_ns_names{gather_order(i - 1)};
        % Find the first match of that name in the unsorted list
        first_match = find(strcmp(last_name, unsorted_names), 1, 'first');
        % If there is no match, get the next unclaimed delta from the
        % unsorted list
        if isempty(first_match)
            first_match = ...
                find(cellfun(@(name) ~isempty(name), unsorted_names),...
                     1,...
                     'first');
        end
        gather_order(i) = first_match;
        % Remove sorted delta from unsorted list
        unsorted_names{first_match} = [];
    end
    if ~isempty(non_state_inds)
        non_state_inds = non_state_inds(gather_order);
    end
    
    % Place state deltas in front of gathered non-state deltas
    new_order = [old_order(state_deltas), non_state_inds];  
    
    % Shift a, b, c matrices
    this_lft = reorderLftDelta(this_lft, new_order);
    
    % Now start gathering together repeated deltas
    [~, unique_indices] = unique(this_lft.delta.names, 'stable');
    gathered_deltas = cell(1, length(unique_indices));
    j = 1;
    for i = unique_indices'
        delta = this_lft.delta.deltas{i};
        delta_group = strcmp(this_lft.delta.names{i}, this_lft.delta.names);
        delta.dim_in  = sum(this_lft.delta.dim_ins(delta_group,:), 1);
        delta.dim_out = sum(this_lft.delta.dim_outs(delta_group,:), 1);
        gathered_deltas{j} = delta;
        j = j + 1;
    end
    this_lft.delta = SequenceDelta(gathered_deltas);
end

total_time = sum(this_lft.horizon_period);
%% Gathering Disturbance
[~, unique_indices] = unique(this_lft.disturbance.names, 'stable');
if (length(this_lft.disturbance.names) ~= length(unique_indices))
    % Unify the same disturbances
    gathered_dis = cell(1, length(unique_indices));
    j = 1;
    for i = unique_indices'
        dis = this_lft.disturbance.disturbances{i};
        group_inds = strcmp(this_lft.disturbance.names{i},...
                               this_lft.disturbance.names);
        dis_group = ...
            SequenceDisturbance(this_lft.disturbance.disturbances{group_inds});
        chan_in = cell(1, total_time);
        for k = 1:total_time
            chan_in{k} = ...
                uniqueStableNonempty(cell2mat(dis_group.chan_ins(:, k)));
        end 
        dis.chan_in  = chan_in;
        gathered_dis{j} = dis;
        j = j + 1;
    end
    this_lft.disturbance = SequenceDisturbance(gathered_dis);
end
% Eliminate default disturbances if others exist
default_disturbances = strcmp(this_lft.disturbance.names, 'default_l2');
if length(default_disturbances) > 1
    disturbances = this_lft.disturbance.disturbances(~default_disturbances);
    this_lft.disturbance = SequenceDisturbance(disturbances);
end

%% Gathering Performance
[~, unique_indices] = unique(this_lft.performance.names, 'stable');
if (length(this_lft.performance.names) ~= length(unique_indices))
    % Unify the same performances
    gathered_perf = cell(1, length(unique_indices));
    j = 1;
    for i = unique_indices'
        perf = this_lft.performance.performances{i};
        group_inds = strcmp(this_lft.performance.names{i},...
                             this_lft.performance.names);
        perf_group = ...
            SequencePerformance(this_lft.performance.performances{group_inds});
        chan_in = cell(1, total_time);
        chan_out = cell(1, total_time);
        for k = 1:total_time
            chan_in{k} = ...
                uniqueStableNonempty(cell2mat(perf_group.chan_ins(:, k)));
            chan_out{k} = ...
                uniqueStableNonempty(cell2mat(perf_group.chan_outs(:, k)));
        end 
        perf.chan_in  = chan_in;
        perf.chan_out = chan_out;
        gathered_perf{j} = perf;
        j = j + 1;
    end
    this_lft.performance = SequencePerformance(gathered_perf);
end
% Eliminate default performances if others exist
default_performances = strcmp(this_lft.performance.names, 'default_l2');
if length(default_performances) > 1
    performances = this_lft.performance.performances(~default_performances);
    this_lft.performance = SequencePerformance(performances);
end
end

function out = uniqueStableNonempty(in)
% This block is defined to avoid changing 0 x 0 empty doubles to 
% 0 x 1 empty doubles[out, ~, inv_order] = unique(in, 'legacy');
out = in;
if ~isempty(out)
    out = unique(in, 'stable');
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)