%% Requirements:
%  1. IQC analysis shall produce an "infeasible problem" result when
%     analyzing repeated SLTV-RB-uncertain-systems that are not robustly stable
%  2. IQC analysis shall produce an upper-bound on worst-case performance
%     for many SLTV-RB-uncertain-systems that are robustly stable. 
%     Producing an upper-bound for ALL SLTV-RB-uncertain-systems is not 
%     expected.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for IQC analysis with norm-bounded operators
classdef testIqcAnalysisSltvRateBnd < matlab.unittest.TestCase
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
    bnd = 1;
    lower_bound = -bnd;
    upper_bound = bnd;
    lower_rate = -bnd;
    upper_rate = bnd;
    lft_bnd = toLft(DeltaSltvRateBnd('test',...
                                      dim_outin,...
                                      lower_bound,...
                                      upper_bound,...
                                      lower_rate,...
                                      upper_rate));
                                        
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
    bnd = 1;
    lower_bound = -bnd;
    upper_bound = bnd;
    lower_rate = -bnd;
    upper_rate = bnd;
    horizon_period = [0, 1];
    basis_length = 2;
    lft_bnd = toLft(DeltaSltvRateBnd('test',...
                                      dim_outin,...
                                      lower_bound,...
                                      upper_bound,...
                                      lower_rate,...
                                      upper_rate,...
                                      horizon_period,...
                                      basis_length));
    
    % Robustly stable interconnection
    lft = interconnect(lft_bnd, lft_g);
    
    opts = AnalysisOptions('verbose', false,...
                           'solver', 'sdpt3',...
                           'lmi_shift', 1e-7);
    [result, valid] = iqcAnalysis(lft, 'analysis_options', opts);
    assertTrue(testCase, valid)
    verifyTrue(testCase, isfinite(result.performance));

    
    % Redefine uncertainty (upper_bound = .9)
    bnd = .9;
    lower_bound = -bnd;
    upper_bound = bnd;
    lower_rate = -2 * bnd;
    upper_rate = 2 * bnd;
    lft_bnd = toLft(DeltaSltvRateBnd('test',...
                                      dim_outin,...
                                      lower_bound,...
                                      upper_bound,...
                                      lower_rate,...
                                      upper_rate));
    
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
    bnd = 1;
    dim_outin = 2;
    lower_bound = -bnd;
    upper_bound = bnd;
    lower_rate = -bnd;
    upper_rate = bnd;
    basis_length = 3;
    horizon_period = [0, 1];
    lft_bnd = toLft(DeltaSltvRateBnd('test',...
                                      dim_outin,...
                                      lower_bound,...
                                      upper_bound,...
                                      lower_rate,...
                                      upper_rate,...
                                      horizon_period,...
                                      basis_length));
                                  
    %  Not robustly stable interconnection
    lft = interconnect(lft_bnd, lft_g);
    
    opts = AnalysisOptions('verbose', false,...
                           'solver', 'sdpt3',...
                           'lmi_shift', 1e-7);
    % Make polytope sltv-rb lft from box sltv-rb
    [~, valid, ~, new_lft] = iqcAnalysis(lft, 'analysis_options', opts);
    vertices = boxToVertices(new_lft.delta.deltas{2}.region_data);
    del_poly = DeltaSltvRateBndImpl('test',...
                                              dim_outin,...
                                              'polytope',...
                                              vertices,...
                                              horizon_period,...
                                              basis_length);
    a = new_lft.a; 
    b = new_lft.b; 
    c = new_lft.c;
    d = new_lft.d;
    delta = SequenceDelta({new_lft.delta.deltas{1}, del_poly});
    new_lft_poly = Ulft(a, b, c, d, delta,...
                        'horizon_period', new_lft.horizon_period,...
                        'performance', new_lft.performance,...
                        'disturbance', new_lft.disturbance);
    [~, valid] = iqcAnalysis(new_lft_poly, 'analysis_options', opts);
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
    
    % Setup for Polytopic rate-bounded uncertainty (upper_bound = 1)
    bnd = 1;
    dim_outin = 2;
    lower_bound = -bnd;
    upper_bound = bnd;
    lower_rate = -bnd;
    upper_rate = bnd;
    basis_length = 3;
    horizon_period = [0, 1];
    lft_bnd = toLft(DeltaSltvRateBnd('test',...
                                      dim_outin,...
                                      lower_bound,...
                                      upper_bound,...
                                      lower_rate,...
                                      upper_rate,...
                                      horizon_period,...
                                      basis_length));
                                  
    %  Intermediate lft representation
    lft = interconnect(lft_bnd, lft_g);
    
    opts = AnalysisOptions('verbose', false,...
                           'solver', 'sdpt3',...
                           'lmi_shift', 1e-7);
    % Make polytope sltv-rb lft from box sltv-rb
    [~, valid, ~, new_lft] = iqcAnalysis(lft, 'analysis_options', opts);
    vertices = boxToVertices(new_lft.delta.deltas{2}.region_data);
    del_poly = DeltaSltvRateBndImpl('test',...
                                              dim_outin,...
                                              'polytope',...
                                              vertices,...
                                              horizon_period,...
                                              basis_length);
    a = new_lft.a; 
    b = new_lft.b; 
    c = new_lft.c;
    d = new_lft.d;
    delta = SequenceDelta({new_lft.delta.deltas{1}, del_poly});
    % Final lft representation
    new_lft_poly = Ulft(a, b, c, d, delta,...
                        'horizon_period', new_lft.horizon_period,...
                        'performance', new_lft.performance,...
                        'disturbance', new_lft.disturbance);
    [result, ~] = iqcAnalysis(new_lft_poly, 'analysis_options', opts);
    is_valid = all(check(result.debug.constraints) >= -1e-8, 'all');
    assertTrue(testCase, is_valid)
    verifyTrue(testCase, isfinite(result.performance));

    
%     % Redefine uncertainty (upper_bound = .9)
%     bnd = 0.9;
%     v1 = [ bnd;  bnd];
%     v2 = [ bnd; -bnd];
%     v3 = [-bnd;  bnd];
%     v4 = [-bnd; -bnd];
%     region_data = {[v1, v2, v3, v4]};
%     lft_bnd = toLft(DeltaSltvRepeated(names,...
%                                       dim_outin,...
%                                       region_type,...
%                                       region_data));
%     
%     % Robustly stable interconnection
%     lft = interconnect(lft_bnd, lft_g);
%     [~, valid] = iqcAnalysis(lft, 'analysis_options', opts);
%     assertTrue(testCase, valid)
%     verifyTrue(testCase, isfinite(result.performance));
%     
%     % Redefine nominal system (hinf_norm = 1)
%     gain = 0.5 / 2;
%     g = [1; 1; 1] * ss(zpk(zero, pole, gain, timestep)) * [1, 1, 1];
%     lft_g = toLft(g);
%     
%     % Robustly stable interconnection
%     lft = interconnect(lft_bnd, lft_g);
%     [~, valid] = iqcAnalysis(lft, 'analysis_options', opts);
%     assertTrue(testCase, valid)
%     verifyTrue(testCase, isfinite(result.performance));
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)