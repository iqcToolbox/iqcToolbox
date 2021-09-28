function this_lft = removeUncertainty(this_lft, delta)
%% REMOVEUNCERTAINTY method for removing specified uncertainties by setting
%     them to zero in the lft.
%
%     lft_out = removeUncertainty(this_lft, 2)
%     lft_out = removeUncertainty(this_lft, 'mass')
%     lft_out = removeUncertainty(this_lft, [2, 6, 3])
%     lft_out = removeUncertainty(this_lft, {'mass', 'Cd', 'Cl'})
%
%     Variables:
%     ---------
%       Input:
%         this_lft : Ulft object :: the lft whose specified uncertainties are to be removed
%         delta : string/cell of strings/array of natural numbers :: the names of the delta objects to be removed 
%       Output:
%         this_lft : Ulft object :: the lft whose specified uncertainties are removed
%
%     See also Ulft, removeDisturbance, removePerformance.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Check and process inputs
if isempty(delta); return; end
validateattributes(delta, {'numeric', 'char', 'cell'}, {'nonempty'})

isChar = @(del) isa(del, 'char');
% Convert all input types to be indices
if isChar(delta)
    delta_ind = find(strcmp(delta, this_lft.delta.names), 1);
    if isempty(delta_ind)
        warning('Ulft:removeUncertainty',...
                ['The named uncertainty, ',...
                 delta,...
                 ', does not appear in the lft'])
    end
elseif isa(delta, 'cell') && all(cellfun(isChar, delta), 'all')
    delta_ind = [];
    for i = 1:length(delta)
        delta_i = delta{i};
        delta_ind = [delta_ind, find(strcmp(delta_i, this_lft.delta.names))];
        if isempty(find(strcmp(delta_i, this_lft.delta.names), 1))
            warning('Ulft:removeUncertainty',...
                    ['The named uncertainty, ',...
                     delta_i,...
                     ', does not appear in the lft'])
        end
    end
elseif isnumeric(delta)
    validateattributes(delta, 'numeric', {'positive',...
                                          'integer',...
                                          'finite',...
                                          'nonnan'})
        bad_inds = delta(delta > length(this_lft.delta.names));
        if ~isempty(bad_inds)
            warning('Ulft:removeUncertainty',...
                    ['The numbered uncertainty, ',...
                     num2str(bad_inds),...
                     ', does not appear in the lft'])
        end
    delta_ind  = delta;
end

%% Remove pertinent sections of a, b, c matrices, reconstruct new lft
% Define the indices of deltas to keep
keep_ind = setdiff(1:length(this_lft.delta.names), delta_ind);

if length(keep_ind) < length(this_lft.delta.names)
    total_time = sum(this_lft.horizon_period);
    a_block = cell(1, total_time);
    b_block = cell(1, total_time);
    c_block = cell(1, total_time);
    [dim_out, dim_in] = size(this_lft);
    % Make separate blocks of matrices pertinent to each delta
    for i = 1:total_time
        a_block{i} = mat2cell(this_lft.a{i},...
                              this_lft.delta.dim_ins(:,i),...
                              this_lft.delta.dim_outs(:,i));
        b_block{i} = mat2cell(this_lft.b{i},...
                              this_lft.delta.dim_ins(:,i),...
                              dim_in(i));
        c_block{i} = mat2cell(this_lft.c{i},...
                              dim_out(i),...
                              this_lft.delta.dim_outs(:,i));
    end

    % Make a, b, c matrices with desired uncertainties cut out  
    a = cell(1, total_time);
    b = cell(1, total_time);
    c = cell(1, total_time);
    for i = 1:total_time
        if ~isempty(keep_ind)
            a{i} = cell2mat(a_block{i}(keep_ind, keep_ind));
            b{i} = cell2mat(b_block{i}(keep_ind, :));
            c{i} = cell2mat(c_block{i}(:, keep_ind));
        else
            a{i} = zeros(0);
            b{i} = zeros(0, dim_in(i));
            c{i} = zeros(dim_out(i), 0);
        end
    end
    
    % Define reduced list of deltas
    new_delta = this_lft.delta.deltas(keep_ind);
    
    % Make reduced lft
    this_lft = Ulft(a, b, c, this_lft.d, new_delta,...
                    'horizon_period', this_lft.horizon_period,...
                    'disturbance', this_lft.disturbance,...
                    'performance', this_lft.performance);    
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)