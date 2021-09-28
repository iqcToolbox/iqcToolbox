%  These tests are for a variety of utility functions which do not need a 
%  complete test suite dedicated to each function.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for utility functions.
classdef testUtilityFunctions < matlab.unittest.TestCase
methods (Test)

function testGetSubclasses(testCase)
% R1: getSubclasses shall accept the char array representing a superclass 
%     and return a cell array of char arrays representing the subclasses stored 
%     in the same directory of the given superclass.
    superclass = 'Delta';
    subclasses = {'DeltaBounded', 'DeltaDelayZ', 'DeltaDlti',...
                  'DeltaIntegrator', 'DeltaSlti', 'DeltaSltv',...
                  'DeltaSltvRateBnd', 'DeltaSltvRateBndImpl',...
                  'DeltaSltvRepeated'};
    subclasses_out = getSubclasses(superclass);
    for i = 1:length(subclasses)
        verifyTrue(testCase, ismember(subclasses{i}, subclasses_out))
    end
end 
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)