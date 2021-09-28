%% Requirements:
%  1. gatherLft shall ensure that duplicates of an Delta object in a Ulft
%       ("duplicates" referring to objects which are exactly the same except 
%       for their dimensions) are gathered to represent a single Delta object
%       with appropriate dimensions
%  2. gatherLft shall ensure that any instance of a DeltaDelayZ or
%       DeltaIntegrator object is listed as the first Delta object in the
%       "delta" property of a Ulft

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for Ulft.gatherLft.
classdef testUlftGatherLft < matlab.unittest.TestCase
methods (Test)
    
function testDeltaGathersDuplicates(testCase)
    % Make LFT where there are multiple duplicates
    Z = DeltaDelayZ();
    d1 = DeltaSlti('a');
    d2 = DeltaSltv('b');
    Z_3_by_3 = DeltaDelayZ(3);
    d1_2_by_2 = DeltaSlti('a', 2);
    
    lft = Ulft(randn(6), randn(6,1), randn(1,6), 0, {d1, d1, Z, d2, Z, Z});
    
    % Check if duplicates have been gathered
    verifyEqual(testCase, length(lft.delta.deltas), 3);
    verifyEqual(testCase, lft.delta.deltas{1}, Z_3_by_3);
    verifyEqual(testCase, lft.delta.deltas{2}, d1_2_by_2);
    verifyEqual(testCase, lft.delta.deltas{3}, d2);
end
    

% % % function testDeltasOrderedAlphanumericallyWithStateFirst(testCase)
% % %     % Make LFT where DeltaDelayZ is given to be last
% % %     delta = SequenceDelta(DeltaSlti('A'),...
% % %                           DeltaBounded('1'),...
% % %                           DeltaDelayZ(),...
% % %                           DeltaSltv('p'));
% % %     a = reshape([1:16], 4, 4);
% % %     b = [1:4]';
% % %     c = [1:4];
% % %     d = 0;
% % %     lft = Ulft(a, b, c, d, delta);
% % %     
% % %     % Check order
% % %     [~, reorder] = sort(delta.names);
% % %     reorder(strcmp('Z', delta.names)) = [];
% % %     reorder = [find(strcmp('Z', delta.names)), reorder];
% % %     
% % %     % Verify Delta order
% % %     for i = 1:length(delta.deltas)
% % %         verifyEqual(testCase, delta.deltas{reorder(i)}, lft.delta.deltas{i})
% % %     end
% % %     
% % %     % Verify correct rearrangement of matrices
% % %     for i = 1:length(delta.deltas)
% % %         % a matrix
% % %         for j = 1:length(delta.deltas)
% % %             verifyEqual(testCase, a(reorder(i), reorder(j)), lft.a{1}(i, j))
% % %         end
% % %         % b matrix
% % %         verifyEqual(testCase, b(reorder(i), 1), lft.b{1}(i, 1))
% % %         % c matrix
% % %         verifyEqual(testCase, c(1, reorder(i)), lft.c{1}(1, i))
% % %     end
% % % end    

% Previous test removed because the following requirement was changed

% % %  3. gatherLft shall ensure that all Delta objects in a Ulft are sorted by
% % %       their names alphanumerically, except that DeltaDelayZ or 
% % %       DeltaIntegrator shall be listed as the first Delta object.
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)