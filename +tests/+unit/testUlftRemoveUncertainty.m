%% Requirements:
%  1. removeUncertainty shall take as an input a user-specified list of 
%      to-be-removed deltas in an LFT, and return as an output an LFT which 
%      is not affected by the given list of deltas.  
%  1.1. removeUncertainty shall take, as an input, a string specifying the 
%      name of a to-be-removed Delta object, a cell of strings specifying 
%      the names of multiple to-be-removed Delta objects, or a numerical 
%      index specifying the indices of the to-be-removed Delta objects.
%  1.2. removeUncertainty shall remove the specified Delta objects from a 
%      Ulft object's delta property.
%  1.3. removeUncertainty shall remove the sections of a Ulft's a, b, and c 
%      properties which pertain to the to-be-removed deltas.
%  2. If the list of to-be-removed deltas includes a delta not present in 
%      the LFT, removeUncertainty will report a warning, but not an error.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for Ulft.
classdef testUlftRemoveUncertainty < matlab.unittest.TestCase
methods (Test)
    
function testOneDeltaRemoval(testCase) 
    % Create lft
    g = drss(3);
    nominal_lft = toLft(g);
    lft = nominal_lft + toLft(DeltaSltv('test')); 
    
    % Variations on arguments for removing an uncertainty
    nom_1 = removeUncertainty(lft, 'test');
    nom_2 = removeUncertainty(lft, {'test'});
    nom_3 = removeUncertainty(lft, 2);
    
    % An incorrect removal
    nom_wrong = removeUncertainty(lft, 1);
    
    verifyEqual(testCase, nom_1, nominal_lft)
    verifyEqual(testCase, nom_2, nominal_lft)
    verifyEqual(testCase, nom_3, nominal_lft)
    verifyNotEqual(testCase, nom_wrong, nominal_lft)    
end

function testMultipleDeltaRemoval(testCase) 
    % Create lft
    g = drss(3);
    nominal_lft = toLft(g);
    lft = nominal_lft + toLft(DeltaSlti('t1')) + toLft(DeltaBounded('t2')); 
    
    % Variations on arguments for removing an uncertainty
    nom_1 = removeUncertainty(lft, {'t1', 't2'});
    nom_2 = removeUncertainty(lft, [2, 3]);
    
    % An incorrect removal
    nom_bad = removeUncertainty(lft, 't1');
    
    verifyEqual(testCase, nom_1, nominal_lft)
    verifyEqual(testCase, nom_2, nominal_lft)
    verifyNotEqual(testCase, nom_bad, nominal_lft)
    
    % Cut out DeltaDelayZ as well
    d_lft = toLft(g.d);
    d_nominal1 = removeUncertainty(lft, {'Z', 't1', 't2'});
    d_nominal2 = removeUncertainty(lft, [1, 2, 3]);
    verifyEqual(testCase, d_nominal1, d_lft)
    verifyEqual(testCase, d_nominal2, d_lft)    
end

function testWarningMissingDeltas(testCase)
    % Create lft
    g = drss(3);
    nominal_lft = toLft(g);
    lft = nominal_lft + toLft(DeltaSlti('t1')) + toLft(DeltaBounded('t2')); 
    
    % Bad non-breaking calls
    bad_text = @() removeUncertainty(lft, 'non-existent');
    bad_cell = @() removeUncertainty(lft, {'t1', 'non-existent'});
    bad_idx_out_of_bounds = @() removeUncertainty(lft, [1, 5]);
    
    % Check that bad calls don't break removeUncertainty
    bad_text();
    bad_cell();
    bad_idx_out_of_bounds();
    
    % Setup to warning to throw error
    old_warning_state = warning('query', 'Ulft:removeUncertainty');
    warning('error', 'Ulft:removeUncertainty');
        
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