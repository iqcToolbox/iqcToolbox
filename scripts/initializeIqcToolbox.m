function initializeIqcToolbox()
%% INITIALIZEIQCTOOLBOX function for adding paths of iqcToolbox and deps
%  initializeIqcToolbox()

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Load configuration file
tmp_dir = fullfile(mfilename('fullpath'), '..', '..');
dir_res = dir(tmp_dir);
top_dir = dir_res(1).folder;
cfg_fname = fullfile(top_dir, 'configuration', 'install_configuration.mat');
assert(isfile(cfg_fname),...
       'initializeIqcToolbox:initializeIqcToolbox',...
       ['iqcToolbox has not been correctly installed,',...
        ' run installIqcToolbox first']);

load(fullfile(top_dir, 'configuration', 'install_configuration.mat'),...
     'install_complete',...
     'install_dir_sdpt3',...
     'install_dir_yalmip',...
     'install_dir_lpsolve',...
     'save_added_paths')

assert(install_complete,...
       'initializeIqcToolbox:initializeIqcToolbox',...
       ['iqcToolbox has not been correctly installed,',...
        ' run installIqcToolbox first']);

%% Add iqcToolbox paths
addpath(fullfile(top_dir));
subdirectories = {'configuration',...
                   'sandbox',...
                   'scripts',...
                   'src'};
for i = 1:length(subdirectories)
    addpath(genpath(fullfile(top_dir, subdirectories{i})))
end
addpath(fullfile(top_dir, 'dependencies')); 

%% Add sdpt3 paths
addpath(genpath(install_dir_sdpt3))
rmpath(genpath(fullfile(install_dir_sdpt3, '.git')))

%% Add yalmip paths
addpath(genpath(install_dir_yalmip))
rmpath(genpath(fullfile(install_dir_yalmip, '.git')))

%% Add lpsolve paths, put library files on path
addpath(genpath(install_dir_lpsolve))
if ispc
    try mxlpsolve('make_lp', 1, 1);
    catch
        setenv('PATH', [getenv('PATH'), install_dir_lpsolve, ';'])
        if save_added_paths
            if ~exist(fullfile('C:','Windows','lpsolve55.dll'), 'file')    
                [~, ~] = system(['copy ',...
                                 fullfile(install_dir_lpsolve, 'lpsolve55.dll'),...
                                 ' ',...
                                 fullfile('C:','Windows')]);
            end
            ls(fullfile('C:', 'Windows', 'lp*'))
        end
        curr_dir = pwd;
        cd(install_dir_lpsolve)
        pwd
        cd(curr_dir)
        mxlpsolve('make_lp', 1, 1);
    end
else
    if save_added_paths
        if ~exist(fullfile('/','usr','lib','liblpsolve55.a'), 'file')    
            [~, ~] = system(['sudo cp ',...
                             fullfile(install_dir_lpsolve, 'liblpsolve55.a'),...
                             ' ',...
                             fullfile('/','usr','lib')]);
        end
        if ~exist(fullfile('/','usr','lib','liblpsolve55.so'), 'file')    
            [~, ~] = system(['sudo cp ',...
                             fullfile(install_dir_lpsolve, 'liblpsolve55.so'),...
                             ' ',...
                             fullfile('/','usr','lib')]);
        end
    else
        curr_dir = pwd;
        cd(install_dir_lpsolve)
        mxlpsolve('make_lp', 1, 1);
        cd(curr_dir)
    end
end
mxlpsolve('make_lp', 1, 1);

%% Add iqcToolbox startup path (to avoid using other startup.m in directory
% This should be the last addpath
addpath(genpath(fullfile(top_dir, 'scripts')))

%% Save paths (if specified)
try
    if save_added_paths
        savepath
    else
        throw(MException('initializeIqcToolbox', 'Paths not saved'))
    end
    [msg, id] = lastwarn;
    if strcmp(id, 'MATLAB:SavePath:PathNotSaved')
        throw(MException(id, msg))
    end
catch
    s1 = sprintf([...
        'Current initialization has not saved the added iqcToolbox paths\n',...
        'and this script will need to be rerun at each startup of MATLAB.']);
    warning('on', 'initializeIqcToolbox:initializeIqcToolbox');
    warning('initializeIqcToolbox:initializeIqcToolbox', s1);
    s2 = sprintf([...
        'To avoid rerunning this script each time, call the following: \n',...
        ' >>> \n',...
        ' clear all \n',...
        ' load(''%s'') \n',...
        ' save_added_paths = true; \n',...
        ' save(''%s'') \n',...
        ' cd(''%s'') \n',...
        ' initializeIqcToolbox \n',...
        ' >>> \n',...
        'Issuing the above commands may require admin privileges \n\n'],...
        cfg_fname,...
        cfg_fname,...
        fullfile(top_dir, 'scripts'));
    disp(s2)
end

disp('iqcToolbox initialized, try running tests.run_tests_script to test tbx.')

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Added incorporation of lpsolve - Micah Fry
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)