%% Requirements:
%  1. MultiplierSltv shall represent the multiplier for static, linear,
%        time-varying uncertainties.
%        See Section VI.C in "System Analysis via Integral Quadratic
%        Constraints" (Megretski and Rantzer, 1997)
%        See Section 5.2.1 in "Robust Stability and Performance Analysis
%        Based on Integral Quadratic Constraints" (Veenman, et al, 2016)
%  1.1 MultiplierSltv shall contain the filter, quad, and constraints such
%        that multiplier = filter' * quad * filter, and the constraints on
%        quad guarantee that the DeltaDlti uncertainty delta is within the 
%        set described by the IQC defined by multiplier.
%  1.2 MultiplierSltv shall define a multiplier w/ the following structure:
%                         *                                  
%             /eye   0  \    /q11 q12\   /eye   0  \
%     Π(jw) = |         |  . |       | . |         |
%             \0     eye/    \q21 q22/   \0     eye/
%               filter'        quad         filter
%  1.2 MultiplierSltv shall have a default configuration, such that it may
%        be constructed with a DeltaSltv uncertainty as the only argument
%  1.3 MultiplierSltv shall be able to flexibly construct it's properties
%        under user input, allowing the following properties to vary 
%        according to user inputs:
%             quad_time_varying: by allowing the quadratic to be
%                                time-varying instead of constant
%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for MultiplierSltv.
classdef testMultiplierSltv < matlab.unittest.TestCase
methods (Test)
function testDefaultConstructor(testCase)
    d = DeltaSltv('test');
    m = MultiplierSltv(d);

    % Standard property check
    verifyEqual(testCase, m.name,           d.name)
    verifyEqual(testCase, m.horizon_period, d.horizon_period)
    verifyEqual(testCase, m.upper_bound,    d.upper_bound)
    verifyEqual(testCase, m.lower_bound,    d.lower_bound)
    verifyEqual(testCase, m.dim_outin,      d.dim_out)
    verifyEqual(testCase, m.dim_outin,      d.dim_in)

    % Check defaults
    verifyTrue(testCase, m.quad_time_varying)
    
    % Check filter
    verifyEqual(testCase, lftToSs(m.filter_lft), ss(eye(2)))
end

function testUnequalBoundsError(testCase)
    hp = [2, 3];
    ub = 1;
    lb = -2;
    verifyError(testCase, ...
                @() MultiplierSltv(DeltaSltv('test', 1, lb, ub, hp)),...
                ?MException)
end

function testConstraintParameters(testCase)
    hp = [2, 6];
    ub = linspace(1, 6, sum(hp));
    lb = -ub;
    dim_outin = 1;
    del = DeltaSltv('test', dim_outin, lb, ub, hp);
    quad_time_varying = true;
    mult = MultiplierSltv(del, 'quad_time_varying', quad_time_varying);
    verifyTrue(testCase, mult.quad_time_varying);
    verifyEqual(testCase, length(mult.decision_vars), 2 * sum(hp))
    verifyEqual(testCase, length(mult.constraints), sum(hp))
    
    quad_time_varying = false;
    mult = MultiplierSltv(del, 'quad_time_varying', quad_time_varying);
    verifyFalse(testCase, mult.quad_time_varying);
    verifyEqual(testCase, length(mult.decision_vars), 2)
    verifyEqual(testCase, length(mult.constraints), 1)    
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)