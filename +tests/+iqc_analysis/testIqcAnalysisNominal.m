%% Requirements:
%  1. IQC analysis shall produce an upper-bound on worst-case performance
%     for any nominal system that is stable.
%  2. iqcAnalysis shall first check if the nominal system is unstable 
%     before analyzing the uncertain system.
%  3. If the nominal system is unstable, iqcAnalysis shall output that the
%     uncertain system is not robustly stable.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for IQC analysis for nominal LTI system
classdef testIqcAnalysisNominal < matlab.unittest.TestCase
methods (Test)

function testStableSystemPerformance(testCase)
    zero = [];
    pole = -.5;
    gain = 1;
    timestep = -1;
    g = ss(zpk(zero, pole, gain, timestep));
    lft = toLft(g);
    options = AnalysisOptions('verbose', false);
    [result, valid, ~, ~] = iqcAnalysis(lft, 'analysis_options', options);
    
    assertTrue(testCase, valid)
    verifyLessThan(testCase, abs(norm(g, 'inf') - result.performance), 1e-3);
end

function testUnstableNominalSystem(testCase)
    % Create unstable nominal system
    zero = [];
    pole = 1.5;
    gain = 1;
    timestep = -1;
    g = ss(zpk(zero, pole, gain, timestep));
    lft = toLft(g);
    
    % Check that it is asserted that the system is not stable
    options = AnalysisOptions('verbose', false);
    [result, valid, ~, ~] = iqcAnalysis(lft, 'analysis_options', options);
    verifyFalse(testCase, valid)
    
    % Create uncertain system
    lft = lft + toLft(DeltaSlti('slti'));
    
    % Check that it is asserted that the system is not robustly stable
    [result, valid, ~, ~] = iqcAnalysis(lft, 'analysis_options', options);
    verifyFalse(testCase, valid)
end

function testUnstableSpringMassDamper(testCase)
    % Spring-mass-damper parameters
    dt = 1;
    m = 1;
    c = 1;
    k = 1 + toLft(DeltaSlti('k', 1, -0.5, 0.5));
    
    % Discrete-time abcd matrices 
    % (euler discretization, unstable at slow sampling)
    a = [1,         dt; 
         -k * dt,   1 - c * dt];
    b = [0; dt / m];
    c = [1, 0];
    d = 0;
    
    % Uncertain system
    spring = toLft(a, b, c, d, dt);
    
    % Nominal system
    spring_nom = removeUncertainty(spring, 'k');
    spring_nominal_ss = ss(spring_nom.a{1}, b, c, d, dt);
    
    % Is IQC analysis asserting it is unstable?
    assertTrue(testCase, norm(eig(spring_nominal_ss.a)) >= 1)
    
    options = AnalysisOptions('verbose', false);
    [~, valid, ~, ~] = iqcAnalysis(spring_nom, 'analysis_options', options);
    verifyFalse(testCase, valid)
    
    [~, valid, ~, ~] = iqcAnalysis(spring, 'analysis_options', options);
    verifyFalse(testCase, valid)
end

function testValidFlagSetCorrectly(testCase)
    % This test is developed to check against hotfix-023.  It failed
    % previous to the fix.
    load(fullfile('+tests', '+iqc_analysis', 'dataIqcAnalysisNominal'),...
         'lft_test')
    options = AnalysisOptions('verbose', false);
    [~, valid, ~, ~] = iqcAnalysis(lft_test, 'analysis_options', options);
    verifyTrue(testCase, valid)
end

function testZeroDimensionZ(testCase)
    a_first = zeros(1, 0);
    c_first = zeros(1, 0);
    a = {a_first, 0}; 
    b = {0, 0}; 
    c = {c_first, 0}; 
    d = {1, 1};
    z = DeltaDelayZ([0, 1], -1, [1, 1])
    lft_feedthrough = Ulft(a, b, c, d, z, 'horizon_period', [1, 1]);
    options = AnalysisOptions('verbose', false);
    result = iqcAnalysis(lft_feedthrough, 'analysis_options', options);
    verifyTrue(testCase, result.valid)
    true_gain = 1;
    verifyLessThan(testCase, abs(true_gain - result.performance), 1e-2)
    
    a = {a_first, 0}; 
    b = {1, 1}; 
    c = {c_first, 1}; 
    d = {0, 0};
    lft_delay = Ulft(a, b, c, d, z, 'horizon_period', [1, 1]);
    result = iqcAnalysis(lft_delay, 'analysis_options', options);
    verifyTrue(testCase, result.valid)
    true_gain = 1;
    verifyLessThan(testCase, abs(true_gain - result.performance), 1e-2)
end

function testCorrectObjectiveScaling(testCase)
    % This test would have failed prior to the fix for hotfix-008
    lft = toLft(DeltaDelayZ());
    options = AnalysisOptions('verbose', false, 'lmi_shift', 1e-6);
    result = iqcAnalysis(lft, 'analysis_options', options);
    objective_scale = 1e3;
    mult_perf = MultiplierL2Induced(lft.performance.performances{1}, 1, 1,...
                                    'objective_scaling', objective_scale);
    result2 = iqcAnalysis(lft,...
                          'multipliers_performance', mult_perf,...
                          'analysis_options', options);
    diff_performance = abs(result.performance - result2.performance);
    verifyLessThan(testCase, diff_performance, 1e-3)
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)