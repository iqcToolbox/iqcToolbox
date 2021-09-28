%% Requirements:
%  1. removeDisturbance shall take as an input a user-specified list of 
%      to-be-removed disturbances in an LFT, and return as an output an LFT
%      which does not have those disturbances.
%  1.1. removeDisturbance shall take, as an input, a string specifying the 
%      name of a to-be-removed Disturbance object, a cell of strings 
%      specifying the names of multiple to-be-removed Disturbance objects, 
%      or a numerical index specifying the indices of the to-be-removed 
%      Disturbance objects.
%  1.2. removeDisturbance shall remove the specified Disturbance objects 
%      from a Ulft's disturbance property.
%  2. If the list of to-be-removed disturbances includes a disturbance 
%      not present in the LFT, removeDisturbance will report a warning, but
%      not an error.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for Ulft.
classdef testUlftRemoveDisturbance < matlab.unittest.TestCase

properties
    chars = ['A':'Z', 'a':'z', '0':'9'];
end
methods (Test)
    
function testOneDisturbanceRemovalFromManyDisturbances(testCase) 
    chars = ['A':'Z', 'a':'z', '0':'9'];
    for i = 1:10
        % Create lft and disturbances
        lft = Ulft.random();
        num_disturbance = randi([2, 5]);
        disturbances = cell(1, num_disturbance);
        for j = 1 : num_disturbance - 1
            name = chars(randi(length(chars), [1, 5]));
            disturbances{j} = DisturbanceL2(name);
        end
        disturbances{end} = DisturbanceL2('to_remove');
        lft = lft.addDisturbance(disturbances);
        lft_without_last = ...
            Ulft(lft.a, lft.b, lft.c, lft.d, lft.delta,...
                 'horizon_period', lft.horizon_period,...
                 'disturbance', lft.disturbance.disturbances(1: end - 1));
        
        % Check removal with text input
        lft_removed1 = removeDisturbance(lft, 'to_remove');
        verifyEqual(testCase, lft_removed1, lft_without_last)

        % Check removal with cell input
        lft_removed2 = removeDisturbance(lft, {'to_remove'});
        verifyEqual(testCase, lft_removed2, lft_without_last)
        
        % Check removal with numeric array input
        lft_removed3 = removeDisturbance(lft, num_disturbance);
        verifyEqual(testCase, lft_removed3, lft_without_last)
        
        % An incorrect removal
        lft_wrong_removed = removeDisturbance(lft, 1);
        verifyNotEqual(testCase, lft_wrong_removed, lft_without_last)
    end
end

function testCompleteDisturbanceRemoval(testCase) 
    chars = ['A':'Z', 'a':'z', '0':'9'];
    for i = 1:10
        % Create lft and disturbances
        lft = Ulft.random();
        num_disturbance = randi([2, 5]);
        disturbances = cell(1, num_disturbance);
        for j = 1 : num_disturbance
            name = chars(randi(length(chars), [1, 5]));
            disturbances{j} = DisturbanceL2(name);
        end
        lft = lft.addDisturbance(disturbances);
        lft_without_dis = Ulft(lft.a, lft.b, lft.c, lft.d, lft.delta,...
                           'horizon_period', lft.horizon_period);                 
        
        % Check removal with cell input
        lft_removed1 = removeDisturbance(lft, lft.disturbance.names);
        verifyEqual(testCase, lft_removed1, lft_without_dis)
        
        % Check removal with numeric array input
        lft_removed2 = removeDisturbance(lft, 1:num_disturbance);
        verifyEqual(testCase, lft_removed2, lft_without_dis)
    end
end

function testMultipleDisturbanceRemoval(testCase) 
    chars = ['A':'Z', 'a':'z', '0':'9'];
    for i = 1:10
        % Create lft and disturbances
        lft = Ulft.random();
        num_disturbance = 2 * randi([1, 3]);
        disturbances = cell(1, num_disturbance);
        % Disturbances to remain
        ind_keep = 1 : 2 : num_disturbance;
        for j = ind_keep
            name = chars(randi(length(chars), [1, 5]));
            disturbances{j} = DisturbanceL2(name);
        end
        % Disturbances to remove
        ind_rem = 2 : 2 : num_disturbance;
        for j = ind_rem
            name = ['to_be_removed_', num2str(j)];
            disturbances{j} = DisturbanceL2(name);
        end
        lft = lft.addDisturbance(disturbances);
        lft_rem_dis = Ulft(lft.a, lft.b, lft.c, lft.d, lft.delta,...
                           'horizon_period', lft.horizon_period,...
                           'disturbance', disturbances(ind_keep));
        
        % Check removal with cell input
        lft_removed1 = removeDisturbance(lft, lft.disturbance.names(ind_rem));
        verifyEqual(testCase, lft_removed1, lft_rem_dis)
        
        % Check removal with numeric array input
        lft_removed2 = removeDisturbance(lft, ind_rem);
        verifyEqual(testCase, lft_removed2, lft_rem_dis)
    end    
end

function testWarningMissingDisturbances(testCase)
    chars = ['A':'Z', 'a':'z', '0':'9'];
    % Create lft
    lft = Ulft.random();
    num_disturbance = randi([1, 5]);
    disturbances = cell(1, num_disturbance);
    for j = 1 : num_disturbance
        name = chars(randi(length(chars), [1, 5]));
        disturbances{j} = DisturbanceL2(name);
    end
    lft = lft.addDisturbance(disturbances);
    
    % Bad non-breaking calls
    bad_text = @() removeDisturbance(lft, 'non-existent');
    bad_cell = @() removeDisturbance(lft, {name, 'non-existent'});
    bad_idx_out_of_bounds = @() removeDisturbance(lft,...
                                                  [1, num_disturbance + 1]);
    
    % Check that bad calls don't break removeDisturbance
    bad_text();
    bad_cell();
    bad_idx_out_of_bounds();
    
    % Setup warning to throw error
    old_warning_state = warning('query', 'Ulft:removeDisturbance');
    warning('error', 'Ulft:removeDisturbance');
        
    % Check if error is thrown
    verifyError(testCase, @() bad_text(), ?MException)
    verifyError(testCase, @() bad_cell(), ?MException)    
    verifyError(testCase, @() bad_idx_out_of_bounds(), ?MException)

    % Return warning to previous "watch" state
    warning(old_warning_state)
end
    
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)