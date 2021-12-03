function power_lft = power(this_lft, exponent)
%% .^ POWER overloaded function for Ulft objects.
%
%     power_lft = this_lft .^ exponent
%     power_lft = power(this_lft, exponent)
% 
%     Variables:
%     ---------
%       Input:
%         this_lft : Ulft object :: lft to be raised
%         exponent : natural number :: power to raise the lft
%       Output:
%         power_lft : Ulft object :: lft raised to the given power
%
%     This method is no different than Ulft.mpower
%
%     See also Ulft, mtimes, mpower.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

power_lft = this_lft ^ exponent;
%%  CHANGELOG
% Dec. 02, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)