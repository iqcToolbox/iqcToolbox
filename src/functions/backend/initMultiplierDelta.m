function multipliers = initMultiplierDelta(array_length)
%% INITMULTIPLIERDELTA function for initializing an array of 
%  MultiplierDeltaDefault objects (i.e., an empty array of MultiplierDelta
%  objects)
%
%  multipliers = initMultiplierDelta(array_length)
%
%  Variables:
%  ---------
%     Input:
%       array_length : natural :: the length of the multiplier array
%     Output:
%       multipliers : array of MultiplierDelta objects :: each element is a MultiplierDeltaDefault
%
%  See also MultiplierDelta, MultiplierDeltaDefault

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

validateattributes(array_length, {'numeric'}, {'nonempty','integer','positive'})

multipliers(array_length) = MultiplierDeltaDefault;
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)