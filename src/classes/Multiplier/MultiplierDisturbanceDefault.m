classdef MultiplierDisturbanceDefault < MultiplierDisturbance
%% MUTLIPLIERDISTURBANCEDEFAULT class. Used as a placeholder for empty values in
%  a MultiplierDisturbance array. Extends the base class MultiplierDisturbance
%
%  extended methods:
%    MultiplierDisturbanceDefault() :: Constructor
%
%  See also MultiplierDisturbanceDefault.MultiplierDisturbanceDefault, MultiplierPerformanceDefault, MultiplierDeltaDefault

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

methods
    function this_mult = MultiplierDisturbanceDefault()
    %% MULTIPLIERDISTURBANCEDEFAULT constructor
    %
    %  this_mult = MultiplierDisturbanceDefault() creates the empty default multiplier to fill arrays
    %
    %  Variables:
    %  ---------
    %    Output:
    %       this_mult : MultiplierDisturbanceDefault object
    %
    %  See also MultiplierDisturbanceDefault, MultiplierPerformanceDefault, MultiplierDeltaDefault
    
        [filter.a{1}, filter.b{1},...
         filter.c{1}, filter.d{1}] = deal([]);
        quad.q{1} = [];
        
        this_mult.filter = filter;
        this_mult.quad   = quad;        
    end
end

end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)