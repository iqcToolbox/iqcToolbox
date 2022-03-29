function multipliers = initMultiplierDisturbance(array_length)
%% INITMULTIPLIERDISTURBANCE function for initializing an array of 
%  MultiplierDisturbanceDefault objects (i.e., an empty array of MultiplierDisturbance objects)
%
%  multipliers = initMultiplierDisturbance(array_length)
%
%  Variables:
%  ---------
%     Input:
%       array_length : natural :: the length of the multiplier array
%     Output:
%       multipliers : array of MultiplierDisturbance objects :: each element is a MultiplierDeltaDefault
%
%  See also MultiplierDisturbance, MultiplierDisturbanceDefault

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

validateattributes(array_length, {'numeric'}, {'nonempty','integer','nonnegative'})
if array_length > 0
    multipliers(array_length) = MultiplierDisturbanceDefault;
else
    multipliers  = MultiplierDisturbance.empty();
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)