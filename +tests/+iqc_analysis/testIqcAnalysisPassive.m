%% Requirements:
%  1. IQC analysis shall produce an "infeasible problem" result when analyzing
%      passive uncertain systems that do not robustly satisfy the performance metric (passivity or stable)
%  2. IQC analysis shall produce a positive result for many passive uncertain systems
%      that robustly satisfy the performance metric. A positive result is not expected for every
%      robustly stable/passive system

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for IQC analysis with passive operators
classdef testIqcAnalysisPassive < matlab.unittest.TestCase
methods (Test)
    function testPassiveComponentsMakePassiveInterconnection(testCase)
    % Interconnection of two passive components makes a passive interconnection
        lower_lft = [-toLft(DeltaPassive('del_low')), 2; 2, 0];
        upper_lft = toLft(DeltaPassive('del_up'));
        lft = interconnect(upper_lft, lower_lft);
        lft = lft.addPerformance({PerformancePassive('pass')});
        options = AnalysisOptions('lmi_shift', 0, 'verbose', false);
        result = iqcAnalysis(lft, 'analysis_options', options);
        valid = check(result.debug.constraints('KYP LMI, 1')) > -1e-10;
        verifyTrue(testCase, valid)
    end
    
    function testContinuousTimeStable(testCase)
    % Interconnection of a strictly output passive linear system with a passive operator is robustly stable
        m = 1;
        k1 = 2;
        k2 = 3;
        s = tf('s');
        g_tf = s / (m * s^2 + k2 * s + k1) + 1;
        g = ss([-g_tf, g_tf; 1, 0]);
        g_lft = toLft(g);
        del = DeltaPassive('del_up');
        lft = interconnect(toLft(del), g_lft);
        options = AnalysisOptions('lmi_shift', 1e-9, 'verbose', false);
        result = iqcAnalysis(lft, 'analysis_options', options);
        verifyTrue(testCase, result.valid)
    end
    
    function testDiscreteTimeStable(testCase)
    % Interconnection of a strictly output passive linear system with a passive operator is robustly stable
        g_tf = ss(tf(1/(tf('z') + .5))) + 3;
        g = ss([-g_tf, g_tf; 1, 0]);
        g_lft = toLft(g);
        del = DeltaPassive('del_up');
        lft = interconnect(toLft(del), g_lft);
        options = AnalysisOptions('lmi_shift', 1e-9, 'verbose', false);
        result = iqcAnalysis(lft, 'analysis_options', options);
        verifyTrue(testCase, result.valid)
    end    
end
end

%%  CHANGELOG
% Oct. 7, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)