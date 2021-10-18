function [info, files, msg] = static_analyzer_base()
%% STATIC_ANALYZER_BASE for analyzing all source code
%
%  [info, file] =  tests.static_analyzer_base
%
%  Runs checkcode on all files in src directory
%
%  Variables:
%  ---------
%    Output:
%      info : cell array of structs :: results from checkcode
%      files : cell array of char arrays :: list of filenames pertaining to info

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
%% Check code

[info, files] = checkcode(src_files,...
                          ['-config=', fullfile(this_pth, 'chkcode_rules.txt')]);
msg = checkcode(src_files,...
                          ['-config=', fullfile(this_pth, 'chkcode_rules.txt')], '-string');
%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Added after v0.5.0 - Micah Fry (micah.fry@ll.mit.edu)