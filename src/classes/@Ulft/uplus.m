function plus_lft = uplus(this_lft)
%% + UPLUS overloaded function for Ulft objects.
%
%     this_lft = +this_lft
%
%     Variables:
%     ---------
%       Input:
%         this_lft : Ulft object :: the lft on the right side of the "+" symbol
%       Output:
%         plus_lft : Ulft object :: the resultant lft
%
%     See also Ulft, plus, minus, uminus, mtimes.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

plus_lft = this_lft;
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)