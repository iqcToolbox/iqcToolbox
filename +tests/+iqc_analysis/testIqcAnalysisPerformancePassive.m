%% Requirements:
%  1. IQC analysis shall produce an "infeasible problem" result when analyzing
%      uncertain systems that are not robustly passive
%  2. IQC analysis shall produce a positive result for many uncertain systems
%      that are strictly passive. A positive result is not expected for every
%      strictly passive system.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for IQC analysis with passive performance metrics
classdef testIqcAnalysisPerformancePassive < matlab.unittest.TestCase
methods (Test)
    function testContinuousTimePassive(testCase)
        a = [0 1; -2 -2];
        b = [0; 1];
        c = [-1, 2];
        d = 1.5;
        g = ss(a, b, c, d);
        lft = toLft(g);
        lft = addPerformance(lft, {PerformancePassive('test')});
        options = AnalysisOptions('lmi_shift', 1e-7, 'verbose', false);
        result = iqcAnalysis(lft, 'analysis_options', options);
        verifyTrue(testCase, result.valid)
        verifyEqual(testCase, result.performance, 0)
        
        lft = -lft;
        result = iqcAnalysis(lft, 'analysis_options', options);
        verifyFalse(testCase, result.valid)
    end
    
    function testDiscreteTimePassive(testCase)
        g = ss(tf(1/(tf('z') + .5))) + 3;
        lft = toLft(g);
        lft = addPerformance(lft, {PerformancePassive('test')});
        options = AnalysisOptions('lmi_shift', 1e-7, 'verbose', false);
        result = iqcAnalysis(lft, 'analysis_options', options);
        verifyTrue(testCase, result.valid)
        verifyEqual(testCase, result.performance, 0)
        
        lft = -lft;
        result = iqcAnalysis(lft, 'analysis_options', options);
        verifyFalse(testCase, result.valid)
    end
    
% The previous tests have systems with non-zero D matrices, and check for 
% strict input passivity.  The following systems are zero for some omega \in [-\infty, \infty]
% and they should pass the following tests (checking only passivity, not SIP),
% except double arithmetic in optimization leads to slight violation of constraints
%     function testContinuousTimePassive2(testCase)
%         g = ss(tf(1/(tf('s') + .5)));
%         lft = toLft(g);
%         lft = addPerformance(lft, {PerformancePassive('test')});
%         options = AnalysisOptions('lmi_shift', 0, 'verbose', false, 'solver', 'sdpt3');
%         result = iqcAnalysis(lft, 'analysis_options', options);
%         verifyTrue(testCase, result.valid)
%         verifyEqual(testCase, result.performance, 0)
%         
%         lft = -lft;
%         result = iqcAnalysis(lft, 'analysis_options', options);
%         verifyFalse(testCase, result.valid)
%     end
%     
%     function testDiscreteTimePassive(testCase)
%         g = ss(tf(1/(tf('z') + .5))) + 2;
%         lft = toLft(g);
%         lft = addPerformance(lft, {PerformancePassive('test')});
%         options = AnalysisOptions('lmi_shift', 0, 'verbose', false);
%         result = iqcAnalysis(lft, 'analysis_options', options);
%         verifyTrue(testCase, result.valid)
%         verifyEqual(testCase, result.performance, 0)
%         
%         lft = -lft;
%         result = iqcAnalysis(lft, 'analysis_options', options);
%         verifyFalse(testCase, result.valid)
%     end
end
end

%%  CHANGELOG
% Oct. 7, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)