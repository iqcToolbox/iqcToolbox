%% Requirements:
%  1. removePerformance shall take as an input a user-specified list of 
%      to-be-removed performances in an LFT, and return as an output an LFT
%      which does not have those performances.
%  1.1. removePerformance shall take, as an input, a string specifying the 
%      name of a to-be-removed Performance object, a cell of strings 
%      specifying the names of multiple to-be-removed Performance objects, 
%      or a numerical index specifying the indices of the to-be-removed 
%      Performance objects.
%  1.2. removePerformance shall remove the specified Performance objects 
%      from a Ulft's performance property.
%  2. If the list of to-be-removed performances includes a performance 
%      not present in the LFT, removePerformance will report a warning, but
%      not an error.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for Ulft.
classdef testUlftRemovePerformance < matlab.unittest.TestCase

properties
    chars = ['A':'Z', 'a':'z', '0':'9'];
end
methods (Test)
    
function testOnePerformanceRemovalFromManyPerformances(testCase) 
    chars = ['A':'Z', 'a':'z', '0':'9'];
    for i = 1:10
        % Create lft and performances
        lft = Ulft.random();
        num_performance = randi([2, 5]);
        performances = cell(1, num_performance);
        for j = 1 : num_performance - 1
            name = chars(randi(length(chars), [1, 5]));
            performances{j} = PerformanceL2Induced(name);
        end
        performances{end} = PerformanceL2Induced('to_remove');
        lft = lft.addPerformance(performances);
        lft_without_last = ...
            Ulft(lft.a, lft.b, lft.c, lft.d, lft.delta,...
                 'horizon_period', lft.horizon_period,...
                 'performance', lft.performance.performances(1: end - 1));
        
        % Check removal with text input
        lft_removed1 = removePerformance(lft, 'to_remove');
        verifyEqual(testCase, lft_removed1, lft_without_last)

        % Check removal with cell input
        lft_removed2 = removePerformance(lft, {'to_remove'});
        verifyEqual(testCase, lft_removed2, lft_without_last)
        
        % Check removal with numeric array input
        lft_removed3 = removePerformance(lft, num_performance);
        verifyEqual(testCase, lft_removed3, lft_without_last)
        
        % An incorrect removal
        lft_wrong_removed = removePerformance(lft, 1);
        verifyNotEqual(testCase, lft_wrong_removed, lft_without_last)
    end
end

function testCompletePerformanceRemoval(testCase) 
    chars = ['A':'Z', 'a':'z', '0':'9'];
    for i = 1:10
        % Create lft and performances
        lft = Ulft.random();
        num_performance = randi([2, 5]);
        performances = cell(1, num_performance);
        for j = 1 : num_performance
            name = chars(randi(length(chars), [1, 5]));
            performances{j} = PerformanceL2Induced(name);
        end
        lft = lft.addPerformance(performances);
        lft_without_perf = Ulft(lft.a, lft.b, lft.c, lft.d, lft.delta,...
                           'horizon_period', lft.horizon_period);                 
        
        % Check removal with cell input
        lft_removed1 = removePerformance(lft, lft.performance.names);
        verifyEqual(testCase, lft_removed1, lft_without_perf)
        
        % Check removal with numeric array input
        lft_removed2 = removePerformance(lft, 1:num_performance);
        verifyEqual(testCase, lft_removed2, lft_without_perf)
    end
end

function testMultiplePerformanceRemoval(testCase) 
    chars = ['A':'Z', 'a':'z', '0':'9'];
    for i = 1:10
        % Create lft and performances
        lft = Ulft.random();
        num_performance = 2 * randi([1, 3]);
        performances = cell(1, num_performance);
        % Performances to remain
        ind_keep = 1 : 2 : num_performance;
        for j = ind_keep
            name = chars(randi(length(chars), [1, 5]));
            performances{j} = PerformanceL2Induced(name);
        end
        % Performances to remove
        ind_rem = 2 : 2 : num_performance;
        for j = ind_rem
            name = ['to_be_removed_', num2str(j)];
            performances{j} = PerformanceL2Induced(name);
        end
        lft = lft.addPerformance(performances);
        lft_rem_perf = Ulft(lft.a, lft.b, lft.c, lft.d, lft.delta,...
                           'horizon_period', lft.horizon_period,...
                           'performance', performances(ind_keep));
        
        % Check removal with cell input
        lft_removed1 = removePerformance(lft, lft.performance.names(ind_rem));
        verifyEqual(testCase, lft_removed1, lft_rem_perf)
        
        % Check removal with numeric array input
        lft_removed2 = removePerformance(lft, ind_rem);
        verifyEqual(testCase, lft_removed2, lft_rem_perf)
    end    
end

function testWarningMissingPerformances(testCase)
    chars = ['A':'Z', 'a':'z', '0':'9'];
    % Create lft
    lft = Ulft.random();
    num_performance = randi([1, 5]);
    performances = cell(1, num_performance);
    for j = 1 : num_performance
        name = chars(randi(length(chars), [1, 5]));
        performances{j} = PerformanceL2Induced(name);
    end
    lft = lft.addPerformance(performances);
    
    % Bad non-breaking calls
    bad_text = @() removePerformance(lft, 'non-existent');
    bad_cell = @() removePerformance(lft, {name, 'non-existent'});
    bad_idx_out_of_bounds = @() removePerformance(lft,...
                                                  [1, num_performance + 1]);
    
    % Check that bad calls don't break removePerformance
    bad_text();
    bad_cell();
    bad_idx_out_of_bounds();
    
    % Setup warning to throw error
    old_warning_state = warning('query', 'Ulft:removePerformance');
    warning('error', 'Ulft:removePerformance');
        
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