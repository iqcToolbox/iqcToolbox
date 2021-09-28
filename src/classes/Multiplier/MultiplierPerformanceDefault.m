classdef MultiplierPerformanceDefault < MultiplierPerformance
%% MUTLIPLIERPERFORMANCEDEFAULT class. Used as a placeholder for empty values in
%  a MultiplierPerformance array. Extends the base class MultiplierPerformance
%
%  extended methods:
%    MultiplierPerformanceDefault() :: Constructor
%
%  See also MultiplierPerformanceDefault.MultiplierPerformanceDefault, MultiplierDisturbanceDefault, MultiplierDeltaDefault

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

methods
    
function this_mult = MultiplierPerformanceDefault()
%% MULTIPLIERPERFORMANCEDEFAULT constructor
%
%  this_mult = MultiplierPerformanceDefault() creates the empty default multiplier to fill arrays
%
%  Variables:
%  ---------
%    Output:
%       this_mult : MultiplierPerformanceDefault object
%
%  See also MultiplierPerformanceDefault, MultiplierDisturbanceDefault, MultiplierDeltaDefault
    
[filter.a{1}, filter.b1{1},  filter.b2{1},...
 filter.c1{1}, filter.d11{1}, filter.d12{1},...
 filter.c2{1}, filter.d21{1}, filter.d22{1}] = deal([]);
[quad.q11{1}, quad.q12{1},...
 quad.q21{1}, quad.q22{1}] = deal([]);

this_mult.filter = filter;
this_mult.quad   = quad;
end    
end

end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)