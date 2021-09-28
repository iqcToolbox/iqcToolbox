%% Requirements:
%  1. IQC analysis shall produce an "infeasible problem" result when
%     analyzing SLTI-uncertain-systems that are not robustly stable
%  2. IQC analysis shall produce an upper-bound on worst-case performance
%     for many SLTI-uncertain-systems that are robustly stable. Producing
%     an upper-bound for ALL SLTI-uncertain-systems is not expected.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for IQC analysis with MultiplierSlti class.
classdef testIqcAnalysisSlti < matlab.unittest.TestCase
methods (Test)
function testNotRobustlyStableSystem(testCase)
    bnd = 2;
    dim_outin = 1;
    del = DeltaSlti('test', dim_outin, -bnd, bnd);
    bad_lft = toLft(del, 1, 1, 0, -1);
    options = AnalysisOptions('verbose', false,...
                              'solver', 'sdpt3',...
                              'lmi_shift', 1e-8);
    [~, valid, yalmip_report, ~] = iqcAnalysis(bad_lft,...
                                               'analysis_options', options);
    no_robustness = ~valid || yalmip_report.problem == 1;
    verifyTrue(testCase, no_robustness)
end

function testRobustlyStableSystem(testCase)
    bnd = 0.9;
    dim_outin = 1;
    del = DeltaSlti('test', dim_outin, -bnd, bnd);
    good_lft = toLft(del, 1, 1, 0, -1);
    options = AnalysisOptions('verbose', false,...
                          'solver', 'sdpt3',...
                          'lmi_shift', 1e-6);
    [result, valid, ~, ~] = iqcAnalysis(good_lft,...
                                        'analysis_options', options);
    assertTrue(testCase, valid)
    verifyLessThan(testCase, abs(result.performance - 10), 0.1)

    basis_length = 4;
    basis_poles = -0.35;
    mults(2) = MultiplierSlti(del,...
                              'basis_length', basis_length,...
                              'basis_poles', basis_poles);
    [result, valid, ~, ~] = iqcAnalysis(good_lft,...
                                        'multipliers_delta', mults,...
                                        'analysis_options', options);
    assertTrue(testCase, valid)
    verifyLessThan(testCase, result.performance, 10.1)

    basis_length = 5;
    basis_poles = linspace(.9, -.9, basis_length - 1)';
    mults(2) = MultiplierSlti(del,...
                              'basis_length', basis_length,...
                              'basis_poles', basis_poles);
    [result, valid, ~, ~] = iqcAnalysis(good_lft,...
                                        'multipliers_delta', mults,...
                                        'analysis_options', options);
    assertTrue(testCase, valid)
    verifyLessThan(testCase, result.performance, 10.1)     

    % Setup for feasibility check
    good_lft.performance.performances{1}.gain = 10.09;
    basis_length = 4;
    basis_poles = [.5 + .1i, .5 - .1i];
    mults(2) = MultiplierSlti(del,...
                              'basis_length', basis_length,...
                              'basis_poles', basis_poles,...
                              'constraint_q11_kyp', false);
    [result, valid, ~, ~] = iqcAnalysis(good_lft,...
                                        'multipliers_delta', mults,...
                                        'analysis_options', options);
    assertTrue(testCase, valid)
    verifyLessThan(testCase, result.performance, 10.1)
end

function testNotRobustlyStableContinuousSystem(testCase)
    a = -1 + DeltaSlti('test');
    b = 1;
    c = 1;
    d = 1;
    % A system which is not robustly stable (consider test = 1)
    bad_lft = toLft(a, b, c, d);
    
    % Check if robustly stable
    options = AnalysisOptions('verbose', false,...
                              'solver', 'sdpt3',...
                              'lmi_shift', 1e-6);
    [~, valid, ~, ~] = iqcAnalysis(bad_lft,...
                                   'analysis_options', options);
    no_robustness = ~valid;
    verifyTrue(testCase, no_robustness);
end

function testRobustlyStableContinuousSystem(testCase)
    dim_outin = 1;
    bnd = 1;
    del = DeltaSlti('test', dim_outin, -bnd, bnd);
    a = -2 + del;
    b = 1;
    c = 1;
    d = 0;
    % A system which is not robustly stable (conside test = 1)
    good_lft = toLft(a, b, c, d);
    
    % Check if robustly stable
    options = AnalysisOptions('verbose', false,...
                              'solver', 'sdpt3',...
                              'lmi_shift', 1e-7);
    [result, valid, ~, ~] = iqcAnalysis(good_lft,...
                                        'analysis_options', options);
    verifyTrue(testCase, valid);
    verifyLessThan(testCase, result.performance, 1.01)
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)