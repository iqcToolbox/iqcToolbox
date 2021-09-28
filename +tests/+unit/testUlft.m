%% Requirements:
%  1. Ulft shall express the upper-lft representation of an
%  uncertain/nonlinear system (M,Δ), with additional IQC-based elements, as
%  shown below:
%
%                              ┌───────┐
%                              │       │
%                     ┌────────► delta ├────────┐
%                     │        │       │        │
%                     │        └───────┘        │
%                     │                         │
%                     │   ┌─────────────────┐   │
%                     │   │                 │   │
%                     └───┤   a    │   b    ◄───┘
%                         │        │        │
%                         │ ───────┼─────── │
%                         │        │        │  disturbance
%                   ◄──┬──┤   c    │   d    ◄──┬───
%                      │  │                 │  │
%                      │  └─────────────────┘  │
%                      │      ___________      │
%                      │     /           \     │
%                      └───►X performance X◄───┘
%                            \___________/
%
%  where
%     - "delta" expresses the operators/uncertainties/nonlinearities
%               perturbing the nominal system (for dynamical systems, 
%               "delta" contains the delay operator (Z) or integration 
%               operator (1/s).
%     - "disturance" expresses the sets of signals disturbing the LFT
%     - "a", "b", "c", "d" (cell arrays of matrices) describe how "delta" 
%                          and "disturbance" enter the system dynamics
%     - "horizon_period" describes the horizon and period of discrete-time
%                        eventually periodic systems (LTI systems have
%                        horizon_period = [0, 1])
%
%  1.1 The nominal system for Ulft objects may be continuous-time LTI,
%      discrete-time LTI, or discrete-time LTV (eventually periodic).
%  1.2 Time-varying Ulft objects may also have time-varying dimensions, in
%      both the system dynamics (a, b, c, d), uncertainties (delta),
%      disturbances (disturbance), and performances (performance).
%  2. Ulft shall be constructable with input arguments "a", "b", "c", "d", 
%     and "delta".  The constructor shall allow users to specify any of the
%     other properties, and shall provide default arguments otherwise.
%  3. Ulft objects shall have the following binary operators which result
%     in a Ulft object:
%         - multiplication (*)
%         - subtraction (-)
%         - addition (+)
%         - power (^), with natural exponents
%  4. Ulft objects shall have the following unary operators which result
%     in a Ulft object:
%         - inversion (inv)
%         - unary minus (-)
%  5. Ulft objects shall have the following n-ary operators which result
%     in a Ulft object:
%         - horizontal concatenation (horzcat, or [a, b])
%         - veritcal concatenation (vertcat, or [a; b])
%  6. Ulft objects shall automatically combine the same uncertainty
%     appearing in multiple instances as a single instance with repeated
%     copies
%  
%     

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for Ulft.
classdef testUlft < matlab.unittest.TestCase
methods (Test)
    
function testSimpleConstructor(testCase)
    % Create simple lft with standard Z
    z = DeltaDelayZ();
    a = {0};
    b = {1};
    c = {1};
    d = {0};
    lft = Ulft(a, b, c, d, z);
    
    % Check correctness
    verifyEqual(testCase, lft.a, a)
    verifyEqual(testCase, lft.b, b)
    verifyEqual(testCase, lft.c, c)
    verifyEqual(testCase, lft.d, d)
    
    verifyEqual(testCase, lft.delta.deltas{1}, z);
    
    verifySize(testCase, lft.disturbance.disturbances, [1, 1])
    verifyEqual(testCase, lft.disturbance.disturbances{1}.name, 'default_l2')
    verifyClass(testCase, lft.disturbance.disturbances{1}, 'DisturbanceL2')
    
    verifySize(testCase, lft.performance.performances, [1, 1])
    verifyEqual(testCase, lft.performance.performances{1}.name, 'default_l2')
    verifyClass(testCase,...
                lft.performance.performances{1},...
                'PerformanceL2Induced')
end

function testTimeVaryingConstruction(testCase)
    % Create lft with time-varying Z
    z_dims = [1,2,4,3];
    z_horizon_period = [3,1];
    z_timestep = -1;
    z = DeltaDelayZ(z_dims, z_timestep, z_horizon_period);
    
    len = length(z_dims);
    a = cell(1, len);
    b = cell(1, len);
    c = cell(1, len);
    d = cell(1, len);
    for i=1:length(z_dims) - 1
        a{i} = i * ones(z_dims(i + 1), z_dims(i));
        b{i} = i * ones(z_dims(i + 1), 1);
        c{i} = i * ones(1, z_dims(i));
        d{i} = i * ones(1);
    end
    a{end} = len * eye(z_dims(end - z_horizon_period(2) + 1), z_dims(end));
    b{end} = len * ones(z_dims(end - z_horizon_period(2) + 1), 1);
    c{end} = len * ones(1, z_dims(end));
    d{end} = len * ones(1);

    lft = Ulft(a, b, c, d, z, 'horizon_period', z_horizon_period);
    
    % Check correctness
    verifyEqual(testCase, lft.a, a)
    verifyEqual(testCase, lft.b, b)
    verifyEqual(testCase, lft.c, c)
    verifyEqual(testCase, lft.d, d)
    
    verifyEqual(testCase, lft.delta.deltas{1}, z);
end

function testNaryOperations(testCase)
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
    a_tv{end} = len * eye(z_dims(end - z_hp(2) + 1), z_dims(end));
    b_tv{end} = len * ones(z_dims(end - z_hp(2) + 1), 1);
    c_tv{end} = len * ones(1, z_dims(end));
    d_tv{end} = len * ones(1);

    lft_ztv = Ulft(a_tv, b_tv, c_tv, d_tv, z_tv, 'horizon_period', z_hp);
    
    % Create time-varying bounded nonlinearity (w/ time-varying dimensions)
    bnd_hp = [3, 1];
    bnd_dim_out = [3, 1, 2, 2];
    bnd_dim_in = [2, 1, 3, 4];
    bnd_bound = 1;

    bnd = DeltaBounded('bnd_tv', bnd_dim_out, bnd_dim_in, bnd_bound, bnd_hp);

    len = length(bnd_dim_out);
    a_bnd = cell(1, len);
    b_bnd = cell(1, len);
    c_bnd = cell(1, len);
    d_bnd = cell(1, len);
    for i=1:length(bnd_dim_out)
        a_bnd{i} = ones(bnd_dim_in(i), bnd_dim_out(i));
        b_bnd{i} = ones(bnd_dim_in(i), 1);
        c_bnd{i} = ones(1, bnd_dim_out(i));
        d_bnd{i} = ones(1);
    end

    lft_bnd = Ulft(a_bnd, b_bnd, c_bnd, d_bnd, bnd, 'horizon_period', bnd_hp);
    
    % Combine two time-varying LFTs
    lft_tv_add = lft_ztv + lft_bnd;
    lft_tv_sub = lft_ztv - lft_bnd;
    lft_tv_mult = lft_ztv * lft_bnd;
    lft_tv_horz = [lft_ztv, lft_bnd];
    lft_tv_vert = [lft_ztv; lft_bnd];
    
    % Create time-invariant Z
    z = DeltaDelayZ();
    a_ti = {0};
    b_ti = {1};
    c_ti = {1};
    d_ti = {0};
    lft_zti = Ulft(a_ti, b_ti, c_ti, d_ti, z);
    
    % Create SLTI uncertainty
    a_slti = DeltaSlti('slti');
    b_slti = 1;
    c_slti = 1;
    d_slti = 0;
    timestep = -1;
    lft_slti = toLft(a_slti, b_slti, c_slti, d_slti, timestep);
    
    % Combine two time-invariant LFTs
    lft_ti_add = lft_zti + lft_slti;
    lft_ti_sub = lft_zti - lft_slti;
    lft_ti_mult = lft_zti * lft_slti;
    lft_ti_horz = [lft_zti, lft_slti];
    lft_ti_vert = [lft_zti; lft_slti];
end

function testPlusMultiplicationOperator(testCase)
    z = toLft(DeltaDelayZ(5));
    dlti = toLft(DeltaDlti('dlti', 3, 5));
    slti = toLft(DeltaSlti('slti', 3));
    lft1 = (ones(3, 5) * z) + dlti + (slti * ones(3, 5));
    
    
    verifyEqual(testCase,...
                lft1.performance,...
                SequencePerformance(PerformanceL2Induced('default_l2')))
    verifyEqual(testCase,...
                lft1.disturbance,...
                SequenceDisturbance(DisturbanceL2('default_l2', {[]})))    
end

function testConsistentDeltaAndA(testCase)
% This test was driven by the hotfix-005 issue, wherein Ulft would not
% prevent constructing ill-formed objects. Those ill-formed objects
% occurred when delta was set to empty, and the a matrix was not.
verifyError(testCase,...
            @() Ulft(ones(3,2), ones(3,1), ones(1,2), ones(1,1), {}),...
            ?MException)
end

end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)