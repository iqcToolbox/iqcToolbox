function power_lft = mpower(this_lft, exponent)
%% ^ MPOWER overloaded function for Ulft objects.
%
%     power_lft = this_lft ^ exponent
%     power_lft = mpower(this_lft, exponent)
% 
%     Variables:
%     ---------
%       Input:
%         this_lft : Ulft object :: lft to be raised
%         exponent : natural number :: power to raise the lft
%       Output:
%         power_lft : Ulft object :: lft raised to the given power
%
%     See also Ulft, mtimes, power.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

validateattributes(exponent, 'numeric', {'integer', 'positive'})
power_lft = this_lft;
for i = 1:exponent - 1
    power_lft = power_lft * this_lft;    
end

end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)