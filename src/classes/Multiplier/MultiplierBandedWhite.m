classdef MultiplierBandedWhite < MultiplierDisturbance
%% MULTIPLIERBANDEDWHTIE class for disturbance signals which are white within
%  a given frequency band (DisturbanceBandedWhite).
%  Extends the base class MultiplierDisturbance.
%
%  extended methods:
%    MultiplierBandedWhite(disturbance, varargin) :: Constructor
%
%  extended properties:
%    chan_in : cell array of column vectors of naturals :: The input channels of the LFT pertaining 
%                                                          to this disturbance
%    dim_in : row of naturals :: The input dimension of the LFT pertaining to this disturbance
%    omega : positive double :: The frequency which defines the band [-omega, omega] in which
%                                the disturbance is white
%    poles : array of poles :: A row or column vector of real, stable poles that parameterize
%                               the multiplier
%
%  See also MultiplierBandedWhite.MultiplierBandedWhite

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    chan_in (1, :) cell
    dim_in (1, :) double
    omega (1, 1) double {mustBeNonnegative}
    poles double
end

methods
function this_mult = MultiplierBandedWhite(disturbance, dim_in_lft, discrete, varargin)
%% MULTIPLIERBANDEDWHITE constructor
%
%  multiplier = MultiplierBandedWhite(disturbance, dim_in_lft, discrete, 'poles', -0.5)
%  multiplier = MultiplierBandedWhite(disturbance, dim_in_lft, discrete) assumes the input above
%
%  Variables:
%  ---------
%     Input:
%       disturbance : DisturbanceBandedWhite object
%       dim_in_lft : row of naturals :: input dimension of LFT associated with disturbance
%       discrete : logical :: Whether or not this multiplier is in discrete-time
%     Output:
%       multiplier : MultiplierBandedWhite object
%
%  See also DisturbanceBandedWhite
input_parser = inputParser;
addRequired(input_parser,...
            'disturbance',...
            @(dis) validateattributes(dis,...
                                      {'DisturbanceBandedWhite'},...
                                      {'nonempty'},...
                                      mfilename));
addRequired(input_parser,...
            'dim_in_lft',...
            @(dim) validateattributes(dim,...
                                      {'numeric'},...
                                      {'integer', 'row', 'positive'},...
                                      mfilename));
addRequired(input_parser,...
            'discrete',...
            @(disc) validateattributes(disc, {'logical'}, {'nonempty'}))
addParameter(input_parser,...
             'poles',...
             -0.5,...
             @(poles) validateattributes(poles,...
                                         {'numeric'},...
                                         {'real', 'vector'}));

parse(input_parser, disturbance, dim_in_lft, discrete, varargin{:})
disturbance = input_parser.Results.disturbance;
dim_in_lft  = input_parser.Results.dim_in_lft;
discrete    = input_parser.Results.discrete;
poles       = input_parser.Results.poles;

assert(sum(disturbance.horizon_period) == length(dim_in_lft),...
       'MultiplierBandedWhite:MultiplierBandedWhite',...
       'dim_in_lft must have the same length as total_time of the disturbance')
consistent_dims = cellfun(@(chan, dim) maxEmpty(chan) <= dim,...
                          disturbance.chan_in,...
                          num2cell(dim_in_lft));
assert(all(consistent_dims),...
       'MultiplierBandedWhite:MultiplierBandedWhite',...
       'Disturbance channels exceeds the input dimensions of lft')
assert(discrete,...
       'MultiplierBandedWhite:MultiplierBandedWhite',...
       'Currently only implements this multiplier in discrete-time')  % This can be relaxed in the future
if discrete
    assert(all(abs(poles) < 1),...
           'MultiplierBandedWhite:MultiplierBandedWhite',...
           'For discrete-time, all poles must be within unit circle')
    assert(disturbance.omega <= pi,...
           'MultiplierBandedWhite:MultiplierBandedWhite',...
           ['For discrete-time, frequency band [-omega, omega] is constrained',...
            ' such that omega <= pi'])
    timestep = -1;
end
% Will need to insert conditions for continuous-time

   
% Define multiplier
this_mult.name              = disturbance.name;
this_mult.horizon_period    = disturbance.horizon_period;
this_mult.chan_in           = disturbance.chan_in;
this_mult.omega             = disturbance.omega;
this_mult.dim_in            = dim_in_lft;
this_mult.discrete          = discrete;
this_mult.poles             = poles;

% Define basis function and filter
poles_length = length(this_mult.poles);
poles_cell = reshape(num2cell(this_mult.poles), poles_length, 1);
poles_cell{end + 1, 1} = [];
basis = zpk(cell(poles_length + 1, 1),...
            poles_cell,...
            ones(poles_length + 1, 1),...
            timestep);
basis_lft = toLft(ss(basis));
rows = eye(this_mult.dim_in(1));
chan_select = toLft(rows(this_mult.chan_in{1}, :));
chan_select = chan_select.matchHorizonPeriod(this_mult.horizon_period);
filter_lft = basis_lft * chan_select;
filter.a = filter_lft.a;
filter.b = filter_lft.b;
filter.c = filter_lft.c;
filter.d = filter_lft.d;

% Define decision vars and quad
decision_vars = cell(1, poles_length + 1);
decision_constraint = 0;
for i = 1:(poles_length + 1)
    x = sdpvar(1);
    omega = this_mult.omega;
    if i == 1
        term = omega;
    else
        pole = this_mult.poles(i - 1);
        if pole == 0
            term = sin(omega);
        else
            term = -atan(-pole * sin(omega) / (1 - pole * cos(omega))) / pole;
        end
    end
    decision_constraint = decision_constraint + term * x;
    decision_vars{i} = x;
end
constraint = (decision_constraint >= 0):['Banded White Multiplier, ',...
                                         this_mult.name,...
                                         ', integral(xi, omega) >= 0'];         %#ok<BDSCA>
q = [zeros(poles_length),           vertcat(decision_vars{2:end});
     horzcat(decision_vars{2:end}), 2 * decision_vars{1}];
quad.q = repmat({q}, 1, sum(this_mult.horizon_period));

% Collect into multiplier
this_mult.filter        = filter;
this_mult.quad          = quad;
this_mult.decision_vars = decision_vars;
this_mult.constraints   = constraint;
end
end
end

%%  CHANGELOG
% Nov. 23, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)