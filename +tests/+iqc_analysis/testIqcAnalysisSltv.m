%% Requirements:
%  1. IQC analysis shall produce an "infeasible problem" result when
%     analyzing SLTV-uncertain-systems that are not robustly stable
%  2. IQC analysis shall produce an upper-bound on worst-case performance
%     for many SLTV-uncertain-systems that are robustly stable. 
%     Producing an upper-bound for ALL norm-bounded-uncertin-systems is not 
%     expected.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for IQC analysis with norm-bounded operators
classdef testIqcAnalysisSltv < matlab.unittest.TestCase
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
    lft_bnd = toLft(DeltaSltv('bnd', dim_outin, -upper_bound, upper_bound));
    
    % Not robustly stable interconnection
    lft = interconnect(lft_bnd, lft_g);
    
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
    lft_bnd = toLft(DeltaSltv('bnd', dim_outin, -upper_bound, upper_bound));
    
    % Robustly stable interconnection
    lft = interconnect(lft_bnd, lft_g);
    
    opts = AnalysisOptions('verbose', false,...
                           'solver', 'sdpt3',...
                           'lmi_shift', 1e-6);
    [result, valid] = iqcAnalysis(lft, 'analysis_options', opts);
    assertTrue(testCase, valid)
    verifyTrue(testCase, isfinite(result.performance));

    
    % Redefine uncertainty (upper_bound = .9)
    upper_bound = .9;
    lft_bnd = toLft(DeltaSltv('bnd', dim_outin, -upper_bound, upper_bound));
    
    % Robustly stable interconnection
    lft = interconnect(lft_bnd, lft_g);
    [~, valid] = iqcAnalysis(lft, 'analysis_options', opts);
    assertTrue(testCase, valid)
    verifyTrue(testCase, isfinite(result.performance));
    
    % Redefine nominal system (hinf_norm = 1)
    gain = 0.5;
    lft_g = toLft([1; 1] * ss(zpk(zero, pole, gain, timestep)) * [1, 1]);
    
    % Robustly stable interconnection
    lft = interconnect(lft_bnd, lft_g);
    [~, valid] = iqcAnalysis(lft, 'analysis_options', opts);
    assertTrue(testCase, valid)
    verifyTrue(testCase, isfinite(result.performance));
end

function testRobustlyStableTimeInvariantSystem(testCase)
    zero = [];
    pole = -.5;
    gain = 1;
    timestep = -1;
    g = ss(zpk(zero, pole, gain, timestep));
    norm_g = norm(g, 'inf');
    lft_g = toLft(g);
    
    dim_outin = 1;
    upper_bound = 2;
    lft_bnd = toLft(DeltaSltv('bnd', dim_outin, -upper_bound, upper_bound));
    
    opts = AnalysisOptions('verbose', false,...
                           'solver', 'sdpt3',...
                           'lmi_shift', 1e-7);
    [result, valid] = iqcAnalysis(lft_g + lft_bnd, 'analysis_options', opts);
    assertTrue(testCase, valid)
    verifyLessThan(testCase,...
                   abs(result.performance - (norm_g + upper_bound)),...
                   1e-3)
   
    [result, valid] = iqcAnalysis(lft_g * lft_bnd, 'analysis_options', opts);
    assertTrue(testCase, valid)
    verifyLessThan(testCase,...
                   abs(result.performance - (norm_g * upper_bound)),...
                   1e-3)
end

function testRobustlyStableTimeVaryingSystem(testCase)
    % Create time-varying Z
    z_dims = [1, 2, 4, 3];
    z_hp = [3, 1];
    z_timestep = -1;
    z_tv = DeltaDelayZ(z_dims, z_timestep, z_hp);
    
    len = length(z_dims);
    a_tv = cell(1, len);
    b_tv = cell(1, len);
    c_tv = cell(1, len);
    d_tv = cell(1, len);
    for i=1:length(z_dims) - 1
        a_tv{i} = i * ones(z_dims(i + 1), z_dims(i));
        b_tv{i} = i * ones(z_dims(i + 1), 1);
        c_tv{i} = i * ones(1, z_dims(i));
        d_tv{i} = i * ones(1);
    end
    a_tv{end} = 0.5 * eye(z_dims(end - z_hp(2) + 1));
    b_tv{end} = len * ones(z_dims(end - z_hp(2) + 1), 1);
    c_tv{end} = len * ones(1, z_dims(end));
    d_tv{end} = len * ones(1);

    lft_ztv = Ulft(a_tv, b_tv, c_tv, d_tv, z_tv, 'horizon_period', z_hp);
    scale_gain = 1 / 100;
    lft_ztv = lft_ztv * matchHorizonPeriod(toLft(scale_gain),...
                                           lft_ztv.horizon_period);
    
    % Create time-varying bounded operator
    dim_outin = ones(1, 4);
    upper_bound = 2;
    bnd_hp = [3, 1];
    lft_sltv = toLft(DeltaSltv('sltv',...
                               dim_outin,...
                               -upper_bound,...
                               upper_bound,...
                               bnd_hp));
    
    opts = AnalysisOptions('verbose', false,...
                           'solver', 'sdpt3',...
                           'lmi_shift', 1e-7);
    [result, valid] = iqcAnalysis(lft_ztv, 'analysis_options', opts);
    assertTrue(testCase, valid);
    nominal_norm = result.performance;

    [result, valid] = iqcAnalysis(lft_ztv + lft_sltv, 'analysis_options', opts);
    assertTrue(testCase, valid);
    expected_perf = nominal_norm + upper_bound;
    verifyLessThan(testCase, result.performance, expected_perf)
    % verifyLessThan(testCase,...
    %                abs(nominal_norm + upper_bound - result.performance),...
    %                1e-3)
    % Note: commented test discarded because measured performance need not
    % equal upper bound on performance (due to triangle inequality)
    
    [result, valid] = iqcAnalysis(lft_ztv * lft_sltv, 'analysis_options', opts);
    assertTrue(testCase, valid);
    expected_perf = nominal_norm * upper_bound;
    percent_error = abs(result.performance - expected_perf) / expected_perf;
    verifyLessThan(testCase, percent_error, 0.01)
end

function testTimeVaryingGain(testCase)
    hp = [3, 1];
    lft_delay = toLft(DeltaDelayZ(1, -1, hp));
    
    bnd = 10;
    lft_sltv  = toLft(DeltaSltv('sltv', 1, -bnd, bnd, hp));
    
    opts = AnalysisOptions('verbose', false,...
                           'solver', 'sdpt3',...
                           'lmi_shift', 1e-7);
    [result, valid] = iqcAnalysis(lft_delay + lft_sltv, 'analysis_options', opts);
    assertTrue(testCase, valid)
    assertLessThan(testCase, abs(bnd + 1 - result.performance), 1e-3)
    [result, valid] = iqcAnalysis(lft_delay * lft_sltv, 'analysis_options', opts);
    assertTrue(testCase, valid)
    assertLessThan(testCase, abs(bnd * 1 - result.performance), 1e-3)
    
    bnd = 10  * [ones(1, hp(1)), zeros(1, hp(2))];
    lft_sltv  = toLft(DeltaSltv('sltv', 1, -bnd, bnd, hp));
    
    [result, valid] = iqcAnalysis(lft_delay * lft_sltv, 'analysis_options', opts);
    assertTrue(testCase, valid)
    verifyLessThan(testCase, abs(bnd(1) * 1 - result.performance), 1e-3)
end

function testTimeVaryingOrConstantDecisionVariables(testCase)
    hp = [5, 3];
    lft_delay = toLft(DeltaDelayZ(1, -1, hp));
    
    bnd = linspace(0.1, 1, sum(hp));
    lft_sltv  = toLft(DeltaSltv('sltv', 1, -bnd, bnd, hp));
    
    opts = AnalysisOptions('verbose', false,...
                           'solver', 'sdpt3',...
                           'lmi_shift', 1e-7);
    lft = lft_delay + lft_sltv;
    [result, valid] = iqcAnalysis(lft, 'analysis_options', opts);
    assertTrue(testCase, valid)
    performance_tv = result.performance;
    mults(2) = MultiplierSltv(lft.delta.deltas{2}, 'quad_time_varying', false);
    [result, valid] = iqcAnalysis(lft,...
                                  'analysis_options', opts,...
                                  'multipliers_delta', mults);
    assertTrue(testCase, valid)
    verifyLessThanOrEqual(testCase, performance_tv, result.performance)
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)