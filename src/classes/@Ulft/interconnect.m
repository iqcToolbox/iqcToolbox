function lft_out = interconnect(upper_lft, lower_lft)
%% INTERCONNECT function to create a single Ulft object from an upper 
%     interconnection of two Ulft objects.
%
%     lft_out = interconnect(upper_lft, lower_lft)
%
%     Variables:
%     ---------
%       Input:
%         upper_lft : Ulft object :: upper system existing in the desired
%                                    interconnection
%         lower_lft : Ulft object :: lower system existing in the desired
%                                    interconnection
%       Output:
%         lft_out : Ulft object :: resultant lft
%
%           ┌─────────┐
%           │         │
%     ┌─────►upper_lft├─────┐
%     │     │         │     │
%     │     └─────────┘     │
%     │                     │       ┌────\       ┌────────────┐
%     │     ┌─────────┐     │       │ to  }   ◄──┤  lft_out   ◄───
%     └─────┤         ◄─────┘       └────/       └────────────┘
%           │lower_lft│
%     ◄─────┤         ◄──────
%           └─────────┘
%
%     See also Ulft.   

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Check correctness and consistency of inputs
if ~isa(upper_lft, 'Ulft')
    ul = toLft(upper_lft); 
else
    ul = upper_lft;
end
if ~isa(lower_lft, 'Ulft')
    ll = toLft(lower_lft); 
else
    ll = lower_lft;
end

validateattributes(ul, {'Ulft'}, {'nonempty'}, mfilename)
validateattributes(ll, {'Ulft'}, {'nonempty'}, mfilename)

if ~isequal(ul.horizon_period, ll.horizon_period)
    horizon_period = commonHorizonPeriod([ul.horizon_period;
                                          ll.horizon_period]);
    ul = matchHorizonPeriod(ul, horizon_period);
    ll = matchHorizonPeriod(ll, horizon_period);
end

assert(all(ul.horizon_period == ll.horizon_period),...
       'Ulft:interconnect',...
       'Both Ulft objects must have the same horizon_period')
assert(all(size(ul, 1) < fliplr(size(ll, 2))),...
       'Ulft:interconnect',...
       ['upper_lft output channels must be smaller than',...
        ' the input channels of lower_lft'])
assert(all(size(ul, 2) < fliplr(size(ll, 1))),...
       'Ulft:interconnect',...
       ['upper_lft input channels must be smaller than',...
        ' the output channels of lower_lft'])
defaultPerformance = @(lft) length(lft.performance.names) == 1 ...
                            && strcmp(lft.performance.names{1}, 'default_l2');
assert(defaultPerformance(ul) && defaultPerformance(ll),...
       'Ulft:interconnect',...
       'Both Ulft objects must have a default performance measure')
defaultDisturbance = @(lft) length(lft.disturbance.names) == 1 ...
                            && strcmp(lft.disturbance.names{1}, 'default_l2');
assert(defaultDisturbance(ul) && defaultPerformance(ll),...
       'Ulft:interconnect',...
       'Both Ulft objects must have a default disturbance spec')

%% Make intermediate matrices
dim_uin  = size(ul, 2);
dim_uout = size(ul, 1);

total_time = sum(ul.horizon_period);
b1  = cell(1, total_time);
b2  = cell(1, total_time);
c1  = cell(1, total_time);
c2  = cell(1, total_time);
d11 = cell(1, total_time);
d12 = cell(1, total_time);
d21 = cell(1, total_time);
d22 = cell(1, total_time);
q1  = cell(1, total_time);
q2  = cell(1, total_time);
for i = 1:total_time
    b1{i}  = ll.b{i}(:, 1:dim_uout);
    b2{i}  = ll.b{i}(:, dim_uout + 1 : end);
    c1{i}  = ll.c{i}(1:dim_uin, :);
    c2{i}  = ll.c{i}(dim_uin + 1 : end, :);
    d11{i} = ll.d{i}(1:dim_uin, 1:dim_uout);
    d12{i} = ll.d{i}(1:dim_uin, dim_uout + 1 : end);
    d21{i} = ll.d{i}(dim_uin + 1 : end, 1:dim_uout);
    d22{i} = ll.d{i}(dim_uin + 1 : end, dim_uout + 1 : end);
    q1{i}  = inv(eye(dim_uin(i)) - d11{i} * ul.d{i});
    q2{i}  = inv(eye(dim_uout(i)) - ul.d{i} * d11{i});
end

%% Make combined_lft matrices
a = cell(1, total_time);
b = cell(1, total_time);
c = cell(1, total_time);
d = cell(1, total_time);
for i = 1:total_time
    a11 = ll.a{i} + b1{i} * ul.d{i} * q1{i} * c1{i};
    a12 = b1{i} * q2{i} * ul.c{i};
    a21 = ul.b{i} * q1{i} * c1{i};
    a22 = ul.a{i} + ul.b{i} * q1{i} * d11{i} * ul.c{i};
    a{i} = [a11, a12;
            a21, a22];
    b{i} = [b2{i} + b1{i} * ul.d{i} * q1{i} * d12{i};
            ul.b{i} * q1{i} * d12{i}];
    c{i} = [c2{i} + d21{i} * q2{i} * ul.d{i} * c1{i}, d21{i} * q2{i} * ul.c{i}];
    d{i} = d22{i} + d21{i} * ul.d{i} * q1{i} * d12{i};
end
delta = SequenceDelta(ll.delta.deltas, ul.delta.deltas);

lft_out = Ulft(a, b, c, d, delta, 'horizon_period', ul.horizon_period);
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)