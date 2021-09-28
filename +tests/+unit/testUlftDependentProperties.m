%% Requirements:
% 1. Ulft.uncertain shall return "true" if the Ulft has any Deltas which are not
%     DeltaDelayZ or DeltaIntegrator
% 2. Ulft.uncertain shall return "false" if the Ulft has at most a DeltaDelayZ
%     or DeltaIntegrator
% 3. Ulft.timestep shall be empty if the Ulft lacks DeltaDelayZ or DeltaIntegrator
% 4. Ulft.timestep shall, if it has a DeltaDelayZ, return the timestep associated
%     with that DeltaDelayZ
% 5. Ulft.timestep shall, if it has a DeltaIntegrator, return 0

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for Ulft.plus
classdef testUlftDependentProperties < matlab.unittest.TestCase
methods (Test)
    function testUncertainTrue(testCase)
        lft = Ulft.random('num_delta', 2);
        verifyTrue(testCase, lft.uncertain);
        
        lft = Ulft.random('num_delta', 2, 'req_deltas', {'DeltaDelayZ'});
        verifyTrue(testCase, lft.uncertain);
        
        lft = Ulft.random('num_delta', 2, 'req_deltas', {'DeltaIntegrator'});
        verifyTrue(testCase, lft.uncertain)
    end
    
    function testUncertainFalse(testCase)
        lft = Ulft.random('num_delta', 0);
        verifyFalse(testCase, lft.uncertain);
        
        lft = Ulft.random('num_delta', 1, 'req_deltas', {'DeltaDelayZ'});
        verifyFalse(testCase, lft.uncertain);
        
        lft = Ulft.random('num_delta', 1, 'req_deltas', {'DeltaIntegrator'});
        verifyFalse(testCase, lft.uncertain);        
    end
    
    function testTimestep(testCase)
        timestep = -1;
        horizon_period = [2, 3];
        lft = toLft(DeltaDelayZ(1, timestep, horizon_period));
        verifyEqual(testCase, lft.timestep(1), timestep)
        
        timestep = 2.3;
        lft = toLft(DeltaDelayZ(2, timestep));
        verifyEqual(testCase, lft.timestep, timestep);
        
        lft = toLft(DeltaIntegrator());
        verifyEqual(testCase, lft.timestep, 0);
        
        lft = Ulft.random().removeUncertainty(1);
        verifyEmpty(testCase, lft.timestep)
    end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Added after v0.5.0 - Micah Fry (micah.fry@ll.mit.edu)