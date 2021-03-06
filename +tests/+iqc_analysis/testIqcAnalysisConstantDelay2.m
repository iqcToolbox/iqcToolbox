%% Requirements:
%  1. IQC analysis shall produce an "infeasible problem" result when analyzing
%     constant time-delayed uncertain systems that are not robustly stable
%  2. IQC analysis shall produce an upper-bound on worst-case performance
%     for many constant delayed uncertain systems that are robustly stable. 
%     Producing an upper-bound for ALL delayed uncertain systems is not 
%     expected.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for IQC analysis with constant time delays (out = (delay - 1) (in))
classdef testIqcAnalysisConstantDelay2 < matlab.unittest.TestCase

methods (TestMethodSetup)
function seedAndReportRng(testCase)
    seed = floor(posixtime(datetime('now')));
    rng(seed, 'twister');
    diagnose_str = ...
        sprintf(['Random inputs may be regenerated by calling: \n',...
                 '>> rng(%10d) \n',...
                 'before running the remainder of the test''s body'],...
                seed);
    testCase.onFailure(@() fprintf(diagnose_str));
end    
end

methods (Test)
function testConstantDelayContinuousTime(testCase)
    % Check that robust stability is correctly concluded (no exponential rate bound)
    wn = 10; % natural freq rad/s
    zeta = 0.5; % damping coefficient
    s = tf('s');
    g = wn^2 / (s^2 + 2 * zeta * wn * s + wn^2);
    g = ss(g);
    delay_max = 0.1;% used to be 0.1251, with verification that result_no_kyp > result
    del = DeltaConstantDelay2('del', 1, delay_max); % 0.1 second delay
    g_del = (1 + del) * g;
    g_del_cl = interconnect(-1, [1; 1] * g_del * [1, 1]);
    options = AnalysisOptions('verbose', false, 'lmi_shift', 1e-6);
    result = iqcAnalysis(g_del_cl, 'analysis_options', options);
    testCase.verifyTrue(result.valid);
    m(2) = MultiplierConstantDelay2(del,...
                                    'discrete', false,...
                                    'constraint_q2_minus_q1_kyp', false,...
                                    'constraint_q2_plus_q1_kyp', false);
    result_no_kyp = iqcAnalysis(g_del_cl,...
                                'analysis_options', options,...
                                'multipliers_delta', m);
    testCase.verifyTrue(result_no_kyp.valid); 
%     testCase.verifyGreaterThan(result_no_kyp.performance, result.performance) 
% Prior line should come out true, but it is highly sensitive to each platform's solver, therefore, removing this test condition
    
    % Check non-KYP constraints with exponential convergence specification
    delay_max = 0.1;
    del = DeltaConstantDelay2('del', 1, delay_max); % 0.1 second delay
    g_del = (1 + del) * g;
    g_del_cl = interconnect(-1, [1; 1] * g_del * [1, 1]);
    expo = 0.1;
    options.exponential = expo;
    m(2) = MultiplierConstantDelay2(del,...
                                    'discrete', false,...
                                    'basis_poles', -expo * 1.1,...
                                    'exponential', expo,...
                                    'constraint_q2_minus_q1_kyp', false,...
                                    'constraint_q2_plus_q1_kyp', false);
    result = iqcAnalysis(g_del_cl,...
                         'analysis_options', options,...
                         'multipliers_delta', m);
    testCase.verifyTrue(result.valid);     
end

function testConstantDelayDiscreteTime(testCase)
    % Check that robust stability is correctly concluded (no exponential rate bound)
    wn = 10; % natural freq rad/s
    zeta = 0.5; % damping coefficient
    s = tf('s');
    g = ss(wn^2 / (s^2 + 2 * zeta * wn * s + wn^2));
    delay_max_c = 0.1; % 0.1 sec delay
    dt = 0.01; % 100 Hz frequency
    gd = c2d(g, dt);
    delay_max = delay_max_c / dt;
    del = DeltaConstantDelay2('del', 1, delay_max); % 10 step delay
    g_del = (1 + del) * gd;
    g_del_cl = interconnect(-1, [1; 1] * g_del * [1, 1]);
    options = AnalysisOptions('verbose', false, 'lmi_shift', 1e-6);
    result = iqcAnalysis(g_del_cl, 'analysis_options', options);
    testCase.verifyTrue(result.valid); 
    
    % Check that robust stability correctly fails (with big enough delay, no exponential specification)
    delay_max = 20;
    del = DeltaConstantDelay2('del', 1, delay_max); 
    g_del = (1 + del) * gd;
    g_del_cl = interconnect(-1, [1; 1] * g_del * [1, 1]);
    result = iqcAnalysis(g_del_cl, 'analysis_options', options);
    testCase.verifyFalse(result.valid); 
    
end
end
end

%%  CHANGELOG
% Apr. 18, 2022: Added after v0.9.0 - Micah Fry (micah.fry@ll.mit.edu)