%% Requirements:
%  1. IQC analysis shall produce an "infeasible problem" result when
%     analyzing sector-bounded uncertain systems that are not robustly stable
%  2. IQC analysis shall produce an upper-bound on worst-case performance
%     for many sector-bounded uncertain systems that are robustly stable. 
%     Producing an upper-bound for ALL SLTV-RB-uncertain-systems is not 
%     expected.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for IQC analysis with Sector-Bounded operators
classdef testIqcAnalysisSectorBounded < matlab.unittest.TestCase

methods (TestMethodSetup)
function seedAndReportRng(testCase)
    seed = floor(posixtime(datetime('now')));
    rng(seed);
    diagnose_str = ...
        sprintf(['Random inputs may be regenerated by calling: \n',...
                 '>> rng(%10d) \n',...
                 'before running the remainder of the test''s body'],...
                seed);
    testCase.onFailure(@() fprintf(diagnose_str));
end    
end

methods (Test)
function testSectorBoundedAnalysis(testCase)
    % Simple test on sole operator
    dim_outin = randi([1, 10]);
    lower_bound = -10 * rand;
    upper_bound = 10 * rand;
    horizon_period = [randi([0, 10]), randi([1, 10])];
    del = DeltaSectorBounded('test', dim_outin, lower_bound, upper_bound, horizon_period);
    lft = toLft(del);
    options = AnalysisOptions('verbose', false);
    result = iqcAnalysis(lft, 'analysis_options', options);
    worst_bound = max([-lower_bound, upper_bound]);
    diff_perf = abs(result.performance - worst_bound);
    testCase.verifyLessThan(diff_perf / worst_bound, 1e-3)
end

function testOpenLoopUnstableIsNotRobustlyStableWithSaturation(testCase)
    % Create an unstable system
    g = ss(1, [1, 1], [1; 1], 0, -1);
    % Stabilizing controller
    k = -dlqr(g.a, 1, 1, 1);
    dim_outin = 1;
    lower_bound = 0;
    upper_bound = 1;
    sat = DeltaSectorBounded('sat', dim_outin, lower_bound, upper_bound);
    sat_k = sat * k;
    lft = interconnect(sat_k, toLft(g));
    % IQC analysis should fail when checking if nominal system is stable
    options = AnalysisOptions('verbose', true);
    result = iqcAnalysis(lft, 'analysis_options', options);
    testCase.verifyFalse(result.valid)
    % Repeat analysis with different model
    sat_k = k - sat * k;
    lft = interconnect(sat_k, toLft(g));
    % IQC analysis should see that nominal system is now stable, but fail to conclude robust stability
    result = iqcAnalysis(lft, 'analysis_options', options);
    testCase.verifyFalse(result.valid)
end

function testRobustlyStableSystem(testCase)
    % System robustly stable to sector [-10, 50]
    s = tf('s');
    g = toLft(ss(0.1 * s / (2 + s) / (3 + s)));
    del = DeltaSectorBounded('sb', 1, -4, 50);
    lft = interconnect(-del, [1; 1] * g * [1, 1]);
    options = AnalysisOptions('verbose', true, 'lmi_shift', 1e-6, 'solver', 'mosek');
    result = iqcAnalysis(lft, 'analysis_options', options);
    testCase.verifyTrue(result.valid)
    
    % System robustly stable to saturation
    g = toLft(ss(4 / ((s + 1) * (s/2 + 1) * (s/3 + 1))));
    del = toLft(DeltaSectorBounded('sb', 1, 0, 1));
    lft = interconnect(-del, [1; 1] * g * [1, 1]);
    result = iqcAnalysis(lft, 'analysis_options', options);
    testCase.verifyTrue(result.valid)
    
    % But not a bigger sector
    del = toLft(DeltaSectorBounded('sb', 1, 0, 1.2));
    lft = interconnect(-del, [1; 1] * g * [1, 1]);
    result = iqcAnalysis(lft, 'analysis_options', options);
    testCase.verifyFalse(result.valid)
    
    % System robustly stable to sector [0, .45]
    g = toLft(ss((10 - s + s^2) / (25 + s + s^2)));
    del = DeltaSectorBounded('sb', 1, 0, .45);
    lft = interconnect(-del, [1; 1] * g * [1, 1]);
    result = iqcAnalysis(lft, 'analysis_options', options);
    testCase.verifyTrue(result.valid)
    
    % But not a bigger sector
    del = DeltaSectorBounded('sb', 1, 0, .6);
    lft = interconnect(-del, [1; 1] * g * [1, 1]);
    result = iqcAnalysis(lft, 'analysis_options', options);
    testCase.verifyFalse(result.valid)
end
end
end

%%  CHANGELOG
% Nov. 18, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)