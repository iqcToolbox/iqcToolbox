classdef MultiplierDeltaDefault < MultiplierDelta
%% MUTLIPLIERDELTADEFAULT class. Used as a placeholder for empty values in
%  a MultiplierDelta array. Extends the base class MultiplierDelta
%
%  extended methods:
%    MultiplierDeltaDefault(delta) :: Constructor
%
%  See also MultiplierDeltaDefault.MultiplierDeltaDefault, MultiplierPerformanceDefault, MultiplierDisturbanceDefault

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

methods
function this_mult = MultiplierDeltaDefault(delta)
%% MULTIPLIERDELTADEFAULT constructor
%
%  this_mult = MultiplierDeltaDefault() creates the empty default multiplier to fill arrays
%  this_mult = MultiplierDeltaDefault(delta) creates a partially-filled multiplier if
%                                            delta is DeltaDelayZ or DeltaIntegrator
%                                            This is to use the same multiplier when creating
%                                            Multipliers in iqcAnalysis for DeltaDelayZ or DeltaIntegrator
%
%  Variables:
%  ---------
%    Input:
%       delta : DeltaDelayZ or DeltaIntegrator object 
%    Output:
%       this_mult : MultiplierDeltaDefault object
%
%  See also MultiplierDeltaDefault, MultiplierPerformanceDefault, MultiplierDisturbanceDefault

    if nargin == 0
    % Not constructed from a Delta
        [filter.a{1}, filter.b1{1},  filter.b2{1},...
        filter.c1{1}, filter.d11{1}, filter.d12{1},...
        filter.c2{1}, filter.d21{1}, filter.d22{1}] = deal([]);
        [quad.q11{1}, quad.q12{1},...
         quad.q21{1}, quad.q22{1}] = deal([]);

    else
        assert(isa(delta, 'DeltaDelayZ') || isa(delta, 'DeltaIntegrator'),...
               'MultiplierDeltaDefault:MultiplierDeltaDefault',...
               ['If constructing an default Delta multiplier, you must',...
                'only provide DeltaDelayZ or DeltaIntegrator objects'])
        
        this_mult.horizon_period = delta.horizon_period;
        total_time = sum(delta.horizon_period);
        
        [filter.a, filter.b1,  filter.b2,...
        filter.c1, filter.d11, filter.d12,...
        filter.c2, filter.d21, filter.d22] = deal(cell(1,total_time));
        [quad.q11, quad.q12,...
         quad.q21, quad.q22] = deal(cell(1, total_time));             
    end
    this_mult.filter = filter;
    this_mult.quad   = quad;
end    
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)