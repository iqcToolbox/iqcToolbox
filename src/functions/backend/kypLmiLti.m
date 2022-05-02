function lmi_mat = kypLmiLti(filter, quad, kyp_var, exponential)
%% KYPLMILTI function for generating the KYP-based LMI of LTI systems
%  Typically "kyp_var" and "quad" will be an sdpvar object, and the returned
%  LMI matrix will be in terms of those sdpvars.
%
%  lmi_mat = kypLmiLti(filter, quad, kyp_var, exponential)
%
%  The frequency domain inequality:
%     filter(freq)' * quad * filter(freq) < 0 for all freq in domain(freq),
%  holds iff there exists a symmetric matrix kyp_var such that:
%     lmi_mat < 0, where lmi_mat is defined in terms of filter, quad, and kyp_var.
%  This function returns the lmi_mat, given a filter, quad, and kyp_var.  Note,
%  the lmi_mat is different if the filter is continuous-time or discrete-time.
%
%  Variables:
%  ---------
%     Input:
%       filter : ss object :: filter in the KYP frequency domain inequality
%       quad   : sdpvar or double array :: quad in the KYP frequency domain inequality
%       kyp_var: sdpvar or double array :: kyp_var which must exist in equivalent KYP LMI condition
%       exponential : double or empty:: scalar value which determines exponential rate of convergence being
%                               verified (if empty, will set exponential = 1 for discrete-time
%                                                        and exponential = 0 for continuous-time)
%     Output:
%       lmi_mat : sdpvar or double array :: lmi_mat which must be positive or negative definite 
%                                           whenever KYP FDI is positive or negative

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

validateattributes(filter, {'ss'}, {'nonempty'}, mfilename)
assert(isa(quad, 'sdpvar') || isa(quad, 'numeric'),...
       'kypLmiLti:kypLmiLti',...
       'quad must be either an sdpvar or numeric')
assert(isa(kyp_var, 'sdpvar') || isa(kyp_var, 'numeric'),...
       'kypLmiLti:kypLmiLti',...
       'kyp_var must be either an sdpvar or numeric')
if nargin < 4
    exponential = [];
end
validateattributes(exponential, {'numeric'}, {'nonnegative'})

% Set default exponential rate for continuous- or discrete-time systems
if isempty(exponential)
    if filter.Ts == 0
        exponential = 0;
    else
        exponential = 1;
    end
end

abcd = [filter.a, filter.b;
        filter.c, filter.d];
[dim_state, dim_in] = size(filter.b);
if filter.Ts == 0
    kyp_mat = [2 * exponential * kyp_var, kyp_var;
               kyp_var,                   zeros(dim_state)];
    lmi_mat = [eye(dim_state), zeros(dim_state, dim_in); abcd]' * ...
              blkdiag(kyp_mat, quad) * ...
              [eye(dim_state), zeros(dim_state, dim_in); abcd];
else
    lmi_mat = abcd' * blkdiag(kyp_var, quad) * abcd - ...
              blkdiag(exponential ^ 2 * kyp_var, zeros(dim_in));
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)