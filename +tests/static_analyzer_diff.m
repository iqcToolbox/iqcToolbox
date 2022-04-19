function [diff_info, diff_files, new_error] = static_analyzer_diff(old_info, old_files)
%% STATIC_ANALYZER_DIFF for analyzing all source code
%
%  [diff_info, diff_file, new_error] =  tests.static_analyzer_diff(old_info, old_files)
%
%  This determines if there are any different errors arising from the src directory 
%  as compared to the old errors expressed in the old_info cell array. The input variables
%  match the format of the output from checkcode() 
%  (e.g., [old_info, old_files] = checkcode('Ulft'))
%
%  Variables:
%  ---------
%    Input:
%      old_info : cell array of structs :: list of results from checkcode()
%      old_files : cell array of char arrays :: list of filenames pertaining to old_info
%    Output:
%      diff_info : cell array of structs :: results from checkcode that differ from old_info
%      diff_files : cell array of char arrays :: list of filenames pertaining to diff_info
%      new_error : boolean :: flag if analyzed files have new errors

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Get current files
this_pth = mfilename('fullpath');
this_pth(end - length(mfilename):end) =  [];
src_path = fullfile(this_pth, '..', 'src');
src_files = dir(fullfile(src_path, '**', '*.m'));
src_files = cellfun(@(dir, name) fullfile(dir, name),...
                    {src_files.folder}, {src_files.name},...
                    'UniformOutput', false);

%% Determine if new errors from matching files
diff_info = cell(1, length(old_info));
diff_files = cell(1, length(old_info));
new_error = false;
for i = 1:length(old_files)
    ind_new = find(strcmp(old_files{i}, src_files));
    if ~isempty(ind_new)
        info = checkcode(src_files{ind_new},...
                        ['-config=', fullfile(this_pth, 'chkcode_rules.txt')]);
        if ~isequal(info, old_info{i})
            info_no_line = rmfield(info, 'line');
            old_info_no_line = rmfield(old_info{i}, 'line');
            for j = 1:length(info)
                in_old = find(arrayfun(@(old) isequal(info_no_line(j), old),...
                                       old_info_no_line), 1);
                if isempty(in_old)
                    new_error = true;
                    diff_info{i} = info;
                    diff_files{i} = src_files{ind_new};
                    break
                end
            end
        end
        src_files(ind_new) = [];
    end
end
diff_info(cellfun(@(c) isempty(c), diff_info)) = [];
diff_files(cellfun(@(c) isempty(c), diff_files)) = [];

%% Determine if new errors from new files
new_info = checkcode(src_files,...
                    ['-config=', fullfile(this_pth, 'chkcode_rules.txt')]);
bad_inds = find(cellfun(@(info) ~isempty(info), new_info));
if ~isempty(bad_inds)
    new_error = true;
    new_info = reshape(new_info(bad_inds), 1, length(new_info));
    new_file = src_files(bad_inds);
    diff_info = [diff_info, new_info];
    diff_files = [diff_files, new_file];
end

if new_error
    if length(diff_info) == 1
        disp(diff_files{1})
    end    
    checkcode(diff_files, ['-config=',fullfile(this_pth, 'chkcode_rules.txt')]);
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Added after v0.5.0 - Micah Fry (micah.fry@ll.mit.edu)