function inv_lft = inv(this_lft)
%% INV overloaded function for Ulft objects.
%
%     inv_lft = inv(this_lft)
%
%     Variables:
%     ---------
%       Input:
%         this_lft : Ulft object :: the lft to be inverted
%       Output:
%         inv_lft : Ulft :: the inverse lft
%
%     See also Ulft.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Check validity
total_time = sum(this_lft.horizon_period);
for i = 1:total_time
    assert(cond(this_lft.d{i}) < 1e11,...
           'Ulft:inv',...
           'lft is not invertible (a d-matrix is not invertible)')
end
defaultPerformance = @(lft) length(lft.performance.names) == 1 ...
                            && strcmp(lft.performance.names{1}, 'default_l2');
assert(defaultPerformance(this_lft),...
       'Ulft:inv',...
       'Ulft object must have a default performance measure')
defaultDisturbance = @(lft) length(lft.disturbance.names) == 1 ...
                            && strcmp(lft.disturbance.names{1}, 'default_l2');
assert(defaultDisturbance(this_lft),...
       'Ulft:inv',...
       'Ulft object must have a default disturbance spec')
    
%% Invert
inv_lft = this_lft;
for i = 1:sum(this_lft.horizon_period)
    inv_lft.d{i} = inv(this_lft.d{i});
    inv_lft.c{i} = -inv_lft.d{i} * this_lft.c{i};
    inv_lft.b{i} = this_lft.b{i} * inv_lft.d{i};
    inv_lft.a{i} = this_lft.a{i} - this_lft.b{i} * inv_lft.d{i} * this_lft.c{i};
end

end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)