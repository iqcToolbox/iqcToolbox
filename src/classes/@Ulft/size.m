function varargout = size(this_lft, dimension)
%% SIZE overloaded function for calculating size of output/input channnels of Ulft object
%
%     [dim_out, dim_in] = size(this_lft)
%     dim_out = size(this_lft)
%     dim_out_or_in = size(this_lft, dimension)
%
%     Variables:
%     ---------
%       Input:
%         this_lft : Ulft objection :: the lft under consideration
%         dimension: natural number :: desired dimension to measure
%                                      (1 - output dimension, 
%                                       2 - input dimension)
%       Output:
%         varargout : array(s) of naturals :: number of lft input/output channels
%
%     See also Ulft.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Calculate size
dim_out = cellfun(@(d_cell) size(d_cell, 1), this_lft.d);
dim_in  = cellfun(@(d_cell) size(d_cell, 2), this_lft.d);

if nargin == 1
    varargout = {dim_out, dim_in};
else
    validateattributes(dimension,...
                       {'numeric'},...
                       {'nonempty', 'integer', 'positive'})
    if dimension == 1
        varargout = {dim_out};
    elseif dimension == 2
        varargout = {dim_in};
    else
        error('Ulft:size', '2nd argument must be 1 (output) or 2 (input)')
    end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)