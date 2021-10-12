function mult = combineAllMultipliers(mult_del, mult_dis, mult_perf, dim_out)
%% COMBINEALLMULTIPLIERS function to create a single, diagonally-augmented
%  multiplier composed of all uncertainty, disturbance, and performance multipliers.
%
%  mult = combineAllMultipliers(mult_del, mult_dis, mult_perf, dim_out)
%
%  Variables:
%  ---------
%    Input:
%      mult_del : MultiplierDeltaCombined :: MultiplierDelta subclass representing all combined Delta multipliers
%      mult_dis : MultiplierDisturbanceCombined :: MultiplierDisturbance subclass representing all combined Disturbance multipliers
%      mult_perf : MultiplierPerformanceCombined :: MultiplierPerformance subclass representing all combined Performance multipliers
%      dim_out : arry of naturals :: output dimensions of Ulft these multipliers pertain to (size(lft, 1))
%    Output:
%      mult : MultiplierDeltaCombined :: Representation of combined multiplier. The fact that it is a MultiplierDelta object
%                                        has no intrinsic meaning, it is simply an available class 
%
%  See also iqcAnalysis, MultiplierDeltaCombined, MultiplierDisturbanceCombined, MultiplierPerformanceCombined

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

validateattributes(mult_del,...
                   'MultiplierDeltaCombined',...
                   {'nonempty'},...
                   mfilename)
validateattributes(mult_dis,...
                   'MultiplierDisturbanceCombined',...
                   {'nonempty'},...
                   mfilename)
validateattributes(mult_perf,...
                   'MultiplierPerformanceCombined',...
                   {'nonempty'},...
                   mfilename)
validateattributes(dim_out, 'numeric', {'integer', 'positive'}, mfilename)

% Create combined filter, create combined quad
a_del   = mult_del.filter.a;
b1_del  = mult_del.filter.b1;
b2_del  = mult_del.filter.b2;
c1_del  = mult_del.filter.c1;
c2_del  = mult_del.filter.c2;
d11_del = mult_del.filter.d11;
d12_del = mult_del.filter.d12;
d21_del = mult_del.filter.d21;
d22_del = mult_del.filter.d22;

a_perf   = mult_perf.filter.a;
b1_perf  = mult_perf.filter.b1;
b2_perf  = mult_perf.filter.b2;
c1_perf  = mult_perf.filter.c1;
c2_perf  = mult_perf.filter.c2;
d11_perf = mult_perf.filter.d11;
d12_perf = mult_perf.filter.d12;
d21_perf = mult_perf.filter.d21;
d22_perf = mult_perf.filter.d22;

a_dis = mult_dis.filter.a;
b_dis = mult_dis.filter.b;
c_dis = mult_dis.filter.c;
d_dis = mult_dis.filter.d;

q11_del = mult_del.quad.q11;
q12_del = mult_del.quad.q12;
q21_del = mult_del.quad.q21;
q22_del = mult_del.quad.q22;

q11_perf = mult_perf.quad.q11;
q12_perf = mult_perf.quad.q12;
q21_perf = mult_perf.quad.q21;
q22_perf = mult_perf.quad.q22;

q   = mult_dis.quad.q;

total_time = length(a_del);
a   = cell(1, total_time);
b1  = cell(1, total_time);
b2  = cell(1, total_time);
c1  = cell(1, total_time);
c2  = cell(1, total_time);
d11 = cell(1, total_time);
d12 = cell(1, total_time);
d21 = cell(1, total_time);
d22 = cell(1, total_time);
q11 = cell(1, total_time);
q12 = cell(1, total_time);
q21 = cell(1, total_time);
q22 = cell(1, total_time);
for i = 1:total_time
    dim_out_dis = size(c_dis{i}, 1);
    
    % Constructing filter matrices for combined, augmented filter
    a{i}   = blkdiag(a_del{i}, a_perf{i}, a_dis{i});
    b1{i}  = blkdiag(b1_del{i}, b1_perf{i});
    b2{i}  = blkdiag(b2_del{i}, [b2_perf{i}; b_dis{i}]);
    c1{i}  = blkdiag(c1_del{i}, c1_perf{i});
    c2{i}  = blkdiag(c2_del{i}, c2_perf{i}, c_dis{i});
    d11{i} = blkdiag(d11_del{i}, d11_perf{i});
    d12{i} = blkdiag(d12_del{i}, d12_perf{i});
    d21{i} = blkdiag(d21_del{i}, [d21_perf{i}; zeros(dim_out_dis, dim_out(i))]);
    d22{i} = blkdiag(d22_del{i}, [d22_perf{i}; d_dis{i}]);
    
    % Constructing combined, performance-augmented quad
    q11{i} = blkdiag(q11_del{i}, q11_perf{i});
    q12{i} = blkdiag(q12_del{i}, [q12_perf{i}, zeros(dim_out(i), dim_out_dis)]);
    q21{i} = blkdiag(q21_del{i}, [q21_perf{i}; zeros(dim_out_dis, dim_out(i))]);
    q22{i} = blkdiag(q22_del{i}, q22_perf{i}, q{i});
end
filter.a   = a;
filter.b1  = b1;
filter.b2  = b2;
filter.c1  = c1;
filter.c2  = c2;
filter.d11 = d11;
filter.d12 = d12;
filter.d21 = d21;
filter.d22 = d22;

quad.q11 = q11;
quad.q12 = q12;
quad.q21 = q21;
quad.q22 = q22;

discrete = horzcat(mult_perf.discrete, mult_dis.discrete, mult_del.discrete);
if ~isempty(discrete)
    assert(all(discrete(1) == discrete),...
           'combineAllMultipliers:combineAllMultipliers',...
           'Cannot combine Multipliers that are discrete and continuous')
    discrete = discrete(1);
end

mult                = MultiplierPerformanceCombined(MultiplierPerformanceDefault());
mult.filter         = filter;
mult.quad           = quad;
mult.name           = [mult_del.name, mult_perf.name, mult_dis.name];
mult.decision_vars  = [mult_del.decision_vars,...
                       mult_perf.decision_vars,...
                       mult_dis.decision_vars];
mult.constraints    = [mult_del.constraints,...
                       mult_perf.constraints,...
                       mult_dis.constraints];
mult.objective      = mult_perf.objective;
mult.horizon_period = mult_perf.horizon_period;
mult.discrete       = discrete;
end