%% Requirements:
%  1. IQC analysis shall produce an "infeasible problem" result when
%     analyzing DLTI-uncertain-systems that are not robustly stable
%  2. IQC analysis shall produce an upper-bound on worst-case performance
%     for many DLTI-uncertain-systems that are robustly stable. Producing
%     an upper-bound for ALL DLTI-uncertain-systems is not expected.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for IQC analysis with MultiplierDlti class.
classdef testIqcAnalysisDlti < matlab.unittest.TestCase
methods (Test)
function testSmallGainNotRobustlyStableSystem(testCase)
    % Define nominal system (hinf_norm = 1)
    zero = [];
    pole = -.5;
    gain = .5;
    timestep = -1;
    g = [1; 1] * ss(zpk(zero, pole, gain, timestep)) * [1, 1];
    lft_g = toLft(g);
    
    % Define uncertainty (upper_bound = 1)
    dim_outin = 1;
    upper_bound = 1;
    lft_dlti = toLft(DeltaDlti('test', dim_outin, dim_outin, upper_bound));
    
    % Not robustly stable interconnection
    lft = interconnect(lft_dlti, lft_g);
    
    opts = AnalysisOptions('verbose', false,...
                           'solver', 'sdpt3',...
                           'lmi_shift', 1e-7);
    [~, valid] = iqcAnalysis(lft, 'analysis_options', opts);
    assertFalse(testCase, valid)
end

function testSmallGainRobustlyStableSystem(testCase)
    % Redefine nominal system (hinf_norm = .98)
    zero = [];
    pole = -.5;
    gain = .49;
    timestep = -1;
    lft_g = toLft([1; 1] * ss(zpk(zero, pole, gain, timestep)) * [1, 1]);
    
    % Define uncertainty (upper_bound = 1)
    dim_outin = 1;
    upper_bound = 1;
    lft_dlti = toLft(DeltaDlti('test', dim_outin, dim_outin, upper_bound));
    
    % Robustly stable interconnection
    lft = interconnect(lft_dlti, lft_g);
    
    opts = AnalysisOptions('verbose', false,...
                           'solver', 'sdpt3',...
                           'lmi_shift', 1e-7);
    [result, valid] = iqcAnalysis(lft, 'analysis_options', opts);
    assertTrue(testCase, valid)
    verifyTrue(testCase, isfinite(result.performance));

    
    % Redefine uncertainty (upper_bound = .9)
    upper_bound = .9;
    lft_dlti = toLft(DeltaDlti('test', dim_outin, dim_outin, upper_bound));
    
    % Robustly stable interconnection
    lft = interconnect(lft_dlti, lft_g);
    [~, valid] = iqcAnalysis(lft, 'analysis_options', opts);
    assertTrue(testCase, valid)
    verifyTrue(testCase, isfinite(result.performance));
    
    % Redefine nominal system (hinf_norm = 1)
    gain = 0.5;
    lft_g = toLft([1; 1] * ss(zpk(zero, pole, gain, timestep)) * [1, 1]);
    
    % Robustly stable interconnection
    lft = interconnect(lft_dlti, lft_g);
    [~, valid] = iqcAnalysis(lft, 'analysis_options', opts);
    assertTrue(testCase, valid)
    verifyTrue(testCase, isfinite(result.performance));
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)