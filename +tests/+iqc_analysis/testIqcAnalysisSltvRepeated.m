%% Requirements:
%  1. IQC analysis shall produce an "infeasible problem" result when
%     analyzing repeated SLTV-uncertain-systems that are not robustly stable
%  2. IQC analysis shall produce an upper-bound on worst-case performance
%     for many SLTV-uncertain-systems that are robustly stable. 
%     Producing an upper-bound for ALL norm-bounded-uncertin-systems is not 
%     expected.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for IQC analysis with norm-bounded operators
classdef testIqcAnalysisSltvRepeated < matlab.unittest.TestCase
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
    region_type = 'box';
    region_data = {[-upper_bound, upper_bound]};
    lft_bnd = toLft(DeltaSltvRepeated('test',...
                                      dim_outin,...
                                      region_type,...
                                      region_data));
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
    region_type = 'box';
    region_data = {[-upper_bound, upper_bound]};
    lft_bnd = toLft(DeltaSltvRepeated('test',...
                                      dim_outin,...
                                      region_type,...
                                      region_data));
    
    % Robustly stable interconnection
    lft = interconnect(lft_bnd, lft_g);
    
    opts = AnalysisOptions('verbose', false,...
                           'solver', 'sdpt3',...
                           'lmi_shift', 1e-7);
    [result, valid] = iqcAnalysis(lft, 'analysis_options', opts);
    assertTrue(testCase, valid)
    verifyTrue(testCase, isfinite(result.performance));

    
    % Redefine uncertainty (upper_bound = .9)
    upper_bound = .9;
    region_type = 'box';
    region_data = {[-upper_bound, upper_bound]};
    lft_bnd = toLft(DeltaSltvRepeated('test',...
                                      dim_outin,...
                                      region_type,...
                                      region_data));
    
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

function testPolytopeSmallGainNotRobustlyStableSystem(testCase)
    % Define nominal system (hinf_norm = 1)
    zero = [];
    pole = -.5;
    gain = .5;
    timestep = -1;
    g = [1; 1; 1] * ss(zpk(zero, pole, gain, timestep)) * [1, 1, 1];
    lft_g = toLft(g);
    
    % Define uncertainty (upper_bound = 1)
    names = {'a', 'b'};
    dim_outin = [1; 1];
    bnd = 1;
    v1 = [ bnd;  bnd];
    v2 = [ bnd; -bnd];
    v3 = [-bnd;  bnd];
    v4 = [-bnd; -bnd];
    region_type = 'polytope';
    region_data = {[v1, v2, v3, v4]};
    lft_bnd = toLft(DeltaSltvRepeated(names,...
                                      dim_outin,...
                                      region_type,...
                                      region_data));
    % Not robustly stable interconnection
    lft = interconnect(lft_bnd, lft_g);
    
    opts = AnalysisOptions('verbose', false,...
                           'solver', 'sdpt3',...
                           'lmi_shift', 1e-7);
    [~, valid] = iqcAnalysis(lft, 'analysis_options', opts);
    assertFalse(testCase, valid)
end

function testPolytopeSmallGainRobustlyStableSystem(testCase)
    % Redefine nominal system (hinf_norm = .98)
    zero = [];
    pole = -.5;
    gain = .49 / 2;
    timestep = -1;
    g = [1; 1; 1] * ss(zpk(zero, pole, gain, timestep)) * [1, 1, 1];
    lft_g = toLft(g);
    
    % Define uncertainty (upper_bound = 1)
    names = {'a', 'b'};
    dim_outin = [1; 1];
    bnd = 1;
    v1 = [ bnd;  bnd];
    v2 = [ bnd; -bnd];
    v3 = [-bnd;  bnd];
    v4 = [-bnd; -bnd];
    region_type = 'polytope';
    region_data = {[v1, v2, v3, v4]};
    lft_bnd = toLft(DeltaSltvRepeated(names,...
                                      dim_outin,...
                                      region_type,...
                                      region_data));

    % Robustly stable interconnection
    lft = interconnect(lft_bnd, lft_g);
    
    opts = AnalysisOptions('verbose', false,...
                           'solver', 'sdpt3',...
                           'lmi_shift', 1e-7);
    [result, valid] = iqcAnalysis(lft, 'analysis_options', opts);
    assertTrue(testCase, valid)
    verifyTrue(testCase, isfinite(result.performance));

    
    % Redefine uncertainty (upper_bound = .9)
    bnd = 0.9;
    v1 = [ bnd;  bnd];
    v2 = [ bnd; -bnd];
    v3 = [-bnd;  bnd];
    v4 = [-bnd; -bnd];
    region_data = {[v1, v2, v3, v4]};
    lft_bnd = toLft(DeltaSltvRepeated(names,...
                                      dim_outin,...
                                      region_type,...
                                      region_data));
    
    % Robustly stable interconnection
    lft = interconnect(lft_bnd, lft_g);
    [~, valid] = iqcAnalysis(lft, 'analysis_options', opts);
    assertTrue(testCase, valid)
    verifyTrue(testCase, isfinite(result.performance));
    
    % Redefine nominal system (hinf_norm = 1)
    gain = 0.5 / 2;
    g = [1; 1; 1] * ss(zpk(zero, pole, gain, timestep)) * [1, 1, 1];
    lft_g = toLft(g);
    
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
    g = [1; 1] * ss(zpk(zero, pole, gain, timestep)) * [1, 1];
    norm_g = norm(g, 'inf');
    lft_g = toLft(g);
    
    dim_outin = 2;
    upper_bound = 2;
    box = {[-upper_bound, upper_bound;
            -upper_bound, upper_bound]};
    lft_bnd = toLft(DeltaSltvRepeated({'a', 'b'},...
                                      [dim_outin / 2; dim_outin / 2],...
                                      'box',...
                                      box));

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

function testRobustlyStableTimeInvariantSystemPolytope(testCase)
    zero = [];
    pole = -.5;
    gain = 1;
    timestep = -1;
    g = [1; 1] * ss(zpk(zero, pole, gain, timestep)) * [1, 1];
    norm_g = norm(g, 'inf');
    lft_g = toLft(g);
    
    dim_outin = 2;
    upper_bound = 2;
    vertices = {[upper_bound, upper_bound, -upper_bound, -upper_bound;
                 upper_bound,-upper_bound, upper_bound, -upper_bound]};
    lft_bnd = toLft(DeltaSltvRepeated({'a', 'b'},...
                                      [dim_outin / 2; dim_outin / 2],...
                                      'polytope',...
                                      vertices));

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
    dim_outin = 2 * ones(1, len);
    for i=1:length(z_dims) - 1
        a_tv{i} = i * ones(z_dims(i + 1), z_dims(i));
        b_tv{i} = i * ones(z_dims(i + 1), dim_outin(i));
        c_tv{i} = i * ones(dim_outin(i), z_dims(i));
        d_tv{i} = i * ones(dim_outin(i));
    end
    a_tv{end} = 0.5 * eye(z_dims(end - z_hp(2) + 1));
    b_tv{end} = len * ones(z_dims(end - z_hp(2) + 1), dim_outin(i));
    c_tv{end} = len * ones(dim_outin(i), z_dims(end));
    d_tv{end} = len * ones(dim_outin(i));

    lft_ztv = Ulft(a_tv, b_tv, c_tv, d_tv, z_tv, 'horizon_period', z_hp);
    scale_gain = 1 / 100 * eye(2);
    lft_ztv = lft_ztv * matchHorizonPeriod(toLft(scale_gain),...
                                           lft_ztv.horizon_period);
    % Create time-varying bounded operator
    upper_bound = 2;
    bnd_hp = [3, 1];
    total_time = sum(bnd_hp);
    lft_sltv1 = toLft(DeltaSltv('sltv1',...
                               dim_outin / 2,...
                               -upper_bound,...
                               upper_bound,...
                               bnd_hp));
    lft_sltv2 = toLft(DeltaSltv('sltv2',...
                               dim_outin / 2,...
                               -upper_bound,...
                               upper_bound,...
                               bnd_hp));
    lft_sltv = [lft_sltv1, 0; 0, lft_sltv2];
    
    box = [-upper_bound, upper_bound;
           -upper_bound, upper_bound];
    lft_sltvr = toLft(DeltaSltvRepeated({'sltv1', 'sltv2'},...
                                        [dim_outin / 2 ; dim_outin / 2],...
                                        'box',...
                                        repmat({box}, 1, total_time),...
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
    
    [result, valid] = iqcAnalysis(lft_ztv + lft_sltvr, 'analysis_options', opts);
    assertTrue(testCase, valid);
    expected_perf = nominal_norm + upper_bound;
    verifyLessThan(testCase, result.performance, expected_perf)
    % verifyLessThan(testCase,...
    %                abs(nominal_norm + upper_bound - result.performance),...
    %                1e-3)
    % Note: commented test discarded because measured performance need not
    % equal upper bound on performance (due to triangle inequality)
    
    [result, valid] = iqcAnalysis(lft_ztv * lft_sltvr, 'analysis_options', opts);
    assertTrue(testCase, valid);
    expected_perf = nominal_norm * upper_bound;
    percent_error = abs(result.performance - expected_perf) / expected_perf;
    verifyLessThan(testCase, percent_error, 0.01)
end

function testRobustlyStableTimeVaryingSystemPolytope(testCase)
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
    dim_outin = 2 * ones(1, len);
    for i=1:length(z_dims) - 1
        a_tv{i} = i * ones(z_dims(i + 1), z_dims(i));
        b_tv{i} = i * ones(z_dims(i + 1), dim_outin(i));
        c_tv{i} = i * ones(dim_outin(i), z_dims(i));
        d_tv{i} = i * ones(dim_outin(i));
    end
    a_tv{end} = 0.5 * eye(z_dims(end - z_hp(2) + 1));
    b_tv{end} = len * ones(z_dims(end - z_hp(2) + 1), dim_outin(i));
    c_tv{end} = len * ones(dim_outin(i), z_dims(end));
    d_tv{end} = len * ones(dim_outin(i));

    lft_ztv = Ulft(a_tv, b_tv, c_tv, d_tv, z_tv, 'horizon_period', z_hp);
    scale_gain = 1 / 100 * eye(2);
    lft_ztv = lft_ztv * matchHorizonPeriod(toLft(scale_gain),...
                                           lft_ztv.horizon_period);
    
    % Create time-varying bounded operator (not-repeated multiplier)
    upper_bound = 2;
    bnd_hp = [3, 1];
    total_time = sum(bnd_hp);
    lft_sltv1 = toLft(DeltaSltv('sltv1',...
                               dim_outin / 2,...
                               -upper_bound,...
                               upper_bound,...
                               bnd_hp));
    lft_sltv2 = toLft(DeltaSltv('sltv2',...
                               dim_outin / 2,...
                               -upper_bound,...
                               upper_bound,...
                               bnd_hp));
    lft_sltv = [lft_sltv1, 0; 0, lft_sltv2];
    
    % Create time-varying bounded operator (repeated multiplier, polytope)
    vertices = repmat({[upper_bound, upper_bound,-upper_bound,-upper_bound;
                        upper_bound,-upper_bound, upper_bound,-upper_bound]},...
                      1, total_time);
    lft_sltvr = toLft(DeltaSltvRepeated({'sltv1', 'sltv2'},...
                                        [dim_outin / 2 ; dim_outin / 2],...
                                        'polytope',...
                                        vertices,...
                                        bnd_hp));
    
    opts = AnalysisOptions('verbose', false,...
                           'solver', 'sdpt3',...
                           'lmi_shift', 1e-7);
    % Nominal analysis
    [result, valid] = iqcAnalysis(lft_ztv, 'analysis_options', opts);
    assertTrue(testCase, valid);
    nominal_norm = result.performance;

    % Non-repeated multiplier analysis
    [result, valid] = iqcAnalysis(lft_ztv + lft_sltv, 'analysis_options', opts);
    assertTrue(testCase, valid);
    expected_perf = nominal_norm + upper_bound;
    verifyLessThan(testCase, result.performance, expected_perf)
    % verifyLessThan(testCase,...
    %                abs(nominal_norm + upper_bound - result.performance),...
    %                1e-3)
    % Note: commented test discarded because measured performance need not
    % equal upper bound on performance (due to triangle inequality)
    
    % Repeated multiplier (polytope) analysis
    [result, valid] = iqcAnalysis(lft_ztv + lft_sltvr, 'analysis_options', opts);
    assertTrue(testCase, valid);
    expected_perf = nominal_norm + upper_bound;
    verifyLessThan(testCase, result.performance, expected_perf)
    % verifyLessThan(testCase,...
    %                abs(nominal_norm + upper_bound - result.performance),...
    %                1e-3)
    % Note: commented test discarded because measured performance need not
    % equal upper bound on performance (due to triangle inequality)
    
    [result, valid] = iqcAnalysis(lft_ztv * lft_sltvr, 'analysis_options', opts);
    assertTrue(testCase, valid);
    expected_perf = nominal_norm * upper_bound;
    percent_error = abs(result.performance - expected_perf) / expected_perf;
    verifyLessThan(testCase, percent_error, 0.01)
end

function testTimeVaryingOrConstantDecisionVariables(testCase)
    hp = [5, 3];
    dim_outin = 4;
    lft_delay = toLft(DeltaDelayZ(dim_outin, -1, hp));
    
    % Time-varying bounds for a box
    bnd = num2cell(linspace(0.1, 1, sum(hp)));
    scale = 2;
    box = cellfun(@(bnd) [-bnd, bnd; -scale * bnd, scale * bnd],...
                  bnd,...
                  'UniformOutput', false);
    lft_sltv  = toLft(DeltaSltvRepeated({'a', 'b'},...
                                        [dim_outin / 2; dim_outin / 2],...
                                        'box', box,...
                                        hp));
    
    opts = AnalysisOptions('verbose', false,...
                           'solver', 'sdpt3',...
                           'lmi_shift', 1e-7);
    lft = lft_delay + lft_sltv;
    [result, valid] = iqcAnalysis(lft, 'analysis_options', opts);
    assertTrue(testCase, valid)
    performance_tv = result.performance;

    % Compare results when using vertices to characterize parameter space
    vertices = boxToVertices(box);
    lft_sltv  = toLft(DeltaSltvRepeated({'a', 'b'},...
                                        [dim_outin / 2; dim_outin / 2],...
                                        'polytope', vertices,...
                                        hp));
    lft = lft_delay + lft_sltv;
    [result, valid] = iqcAnalysis(lft, 'analysis_options', opts);
    assertTrue(testCase, valid)
    verifyLessThan(testCase, abs(result.performance - performance_tv), 1e-4)

    % Repeat box results for future comparison
    lft_sltv  = toLft(DeltaSltvRepeated({'a', 'b'},...
                                        [dim_outin / 2; dim_outin / 2],...
                                        'box', box,...
                                        hp));
    lft = lft_delay + lft_sltv;
    [result, valid] = iqcAnalysis(lft, 'analysis_options', opts);
    assertTrue(testCase, valid)
    performance_tv = result.performance;

    % Does box result become worse when using time-invariant quad? (Yes)
    mults(2) = MultiplierSltvRepeated(lft.delta.deltas{2},...
                                      'quad_time_varying', false);
    [result, valid] = iqcAnalysis(lft,...
                                  'analysis_options', opts,...
                                  'multipliers_delta', mults);
    assertTrue(testCase, valid)
    verifyLessThanOrEqual(testCase, performance_tv, result.performance)
end

function testTimeVaryingOrConstantDecisionVariablesPolytope(testCase)
    hp = [5, 3];
    dim_outin = 4;
    lft_delay = toLft(DeltaDelayZ(dim_outin, -1, hp));
    
    bnd = num2cell(linspace(0.1, 1, sum(hp)));
    scale = 2;
    box = cellfun(@(bnd) [-bnd, bnd; -scale * bnd, scale * bnd],...
                  bnd,...
                  'UniformOutput', false);
    vertices = boxToVertices(box);
    lft_sltv  = toLft(DeltaSltvRepeated({'a', 'b'},...
                                        [dim_outin / 2; dim_outin / 2],...
                                        'polytope', vertices,...
                                        hp));
    
    opts = AnalysisOptions('verbose', false,...
                           'solver', 'sdpt3',...
                           'lmi_shift', 1e-7);
    lft = lft_delay + lft_sltv;
    [result, valid] = iqcAnalysis(lft, 'analysis_options', opts);
    assertTrue(testCase, valid)
    performance_tv = result.performance;

    % Does polytope result become worse when using time-invariant quad?
    % (Yes)
    mults(2) = MultiplierSltvRepeated(lft.delta.deltas{2},...
                                      'quad_time_varying', false);
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