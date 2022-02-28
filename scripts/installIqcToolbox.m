function installIqcToolbox(save_added_paths, interactive)
%% INSTALLIQCTOOLBOX function for installing the iqcToolbox
%
%  installIqcToolbox(save_added_paths, interactive)
%  installIqcToolbox() assumes save_added_paths == false, interactive == true
%
%  Variables:
%  ---------
%    Input:
%       save_added_paths : boolean :: whether or not the added paths should be permanently saved
%       interactive : boolean :: whether or not the user must interactively affirm installation of dependencies


%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Process inputs
switch nargin
    case 0
        save_added_paths = false;
        interactive = true;
    case 2
        save_added_paths = logical(save_added_paths);
        interactive = logical(interactive);
    otherwise
        error('installIqcToolbox:installIqcToolbox',...
              ['User must provide two boolean inputs or ',...
               'no inputs (thereby permitting default values) ',...
               'to installIqcToolbox(...)'])
end

%% Define iqcToolbox directory
tmp_dir = fullfile(mfilename('fullpath'), '..', '..');
dir_res = dir(tmp_dir);
top_dir = dir_res(1).folder;
%% Check if iqcToolbox is already installed
cfg_fname = fullfile(top_dir, 'configuration', 'install_configuration.mat');
if isfile(cfg_fname)
    load(cfg_fname, 'install_complete')
    if install_complete
        initializeIqcToolbox;
        return
    end
end

%% Check if required dependencies are already installed
% yalmip
try
    sdpvar;
    dependency_installed_yalmip = true;
    tmp_dir = fullfile(which('sdpvar'), '..', '..');
    dir_res = dir(tmp_dir);
    install_dir_yalmip = dir_res(1).folder;
catch
    dependency_installed_yalmip = false;
end

% SDPT3
try
    B = graph2(10);
    feas = 1;
    [blk,At,C,b] = thetaproblem(B, feas);
    ops.printlevel = 0;
    [obj, X, y, Z, info] = sqlp(blk, At, C, b, ops);
    if info.termcode == 0
        dependency_installed_sdpt3 = true;
    end
    tmp_dir = fullfile(which('sqlp'), '..', '..');
    dir_res = dir(tmp_dir);
    install_dir_sdpt3 = dir_res(1).folder;
catch
    dependency_installed_sdpt3 = false;
end

% C/C++ Compiler
try
    mex -setup C
    dependency_installed_ccompiler = true;
catch
    dependency_installed_ccompiler = false;
end

% LPSOLVE
try
    lp=mxlpsolve('make_lp',0,4);
    mxlpsolve('add_constraint',lp,[3, 2, 2, 1],3,4);
    mxlpsolve('add_constraint',lp,[0, 4, 3, 1],2,3);
    mxlpsolve('set_obj_fn',lp,[2, 3, -2, 3]);
    mxlpsolve('set_verbose',lp,3)
    result=mxlpsolve('solve',lp);
    obj=mxlpsolve('get_objective', lp);
    x=mxlpsolve('get_variables', lp);
    dependency_installed_lpsolve = true;
    tmp_dir = fullfile(which('lp_maker'), '..');
    dir_res = dir(tmp_dir);
    install_dir_lpsolve = dir_res(1).folder;
catch
    dependency_installed_lpsolve = false;
end

%% Install essential packages
if ~dependency_installed_yalmip
    try
        if interactive
            affirmInstallation('yalmip')
        end
        disp(['Installing from the fork ',...
              '<a href="https://github.com/iqcToolbox/yalmip">YALMIP</a>, ',...
              'Copyright (c) 2012-2021 by Johan Löfberg'])
        install_dir_yalmip = fullfile(top_dir, 'dependencies', 'yalmip');
        tmp_dir = tempname;
        mkdir(tmp_dir);
        websave(fullfile(tmp_dir, 'yalmip.zip'),...
                'https://github.com/iqcToolbox/yalmip/archive/develop.zip');
        unzip(fullfile(tmp_dir, 'yalmip.zip'), tmp_dir);
        movefile(fullfile(tmp_dir, 'YALMIP-develop'), install_dir_yalmip)
        dependency_installed_yalmip = true;
    catch
        install_dir_yalmip = [];
    end
end

comp = computer;
comp = comp(1 : end - 2);
switch comp
case 'PCWIN'
    if ~dependency_installed_ccompiler
        error('installIqcToolbox:installIqcToolbox',...
              ['In order to compile the SDPT3 dependency, MATLAB ',...
              'must have access to a C/C++ Compiler. \n',...
              'Various free compilers are available for different',...
              ' platforms (Windows, Linux, Mac).',...
              '\n\nFor Windows machines, MATLAB provides a MinGW C/C++ ',...
              ' compiler. Visit and follow the instructions on ',...
              '<a href="https://www.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler?s_tid=prof_contriblnk">this site</a> \n',...
              'Or, alternatively, in the main MATLAB',...
              ' window, go to Environment -> Add-Ons -> Get Add-ons, ',...
              ' & search for the MinGW C/C++ Compiler Add On.',...
              '\n\n More information on other compilers and other platforms',...
              ' may be found ',...
              '<a href="https://www.mathworks.com/support/requirements/supported-compilers.html">here</a>.',...
              ])
    end
    
    if ~dependency_installed_lpsolve
        if interactive
            affirmInstallation('lpsolve')
        end
        disp(['Installing ',...
              '<a href="http://lpsolve.sourceforge.net/">LPSOLVE 5.5</a>, ',...
              'Copyright (c) 2004 - 2021 by ',...
              'M. Berkelaar, K. Eikland, and P. Notebaert'])
        install_dir_lpsolve = fullfile(top_dir, 'dependencies', 'lpsolve');
        tmp_dir = tempname;
        mkdir(tmp_dir);
        lpsolve_ex = 'https://github.com/iqcToolbox/lpsolve552/raw/master/lp_solve_5.5.2.0_exe_win64.zip';
        lpsolve_dev = 'https://github.com/iqcToolbox/lpsolve552/raw/master/lp_solve_5.5.2.0_dev_win64.zip';
        lpsolve_mat = 'https://github.com/iqcToolbox/lpsolve552/raw/master/lp_solve_5.5.2.0_MATLAB_exe_win64.zip';
        websave(fullfile(tmp_dir, 'lpsolve_ex.zip'), lpsolve_ex);
        websave(fullfile(tmp_dir, 'lpsolve_dev.zip'), lpsolve_dev);
        websave(fullfile(tmp_dir, 'lpsolve_mat.zip'), lpsolve_mat);
        unzip(fullfile(tmp_dir, 'lpsolve_ex.zip'), tmp_dir);
        unzip(fullfile(tmp_dir, 'lpsolve_dev.zip'), tmp_dir);
        unzip(fullfile(tmp_dir, 'lpsolve_mat.zip'), tmp_dir);
        % Need to rename mxlpsolve file to interface correctly with YALMIP
        movefile(fullfile(tmp_dir, 'mxlpsolve.m'),...
                 fullfile(tmp_dir, 'mxlpsolve_old.m'))
        movefile(fullfile(tmp_dir), install_dir_lpsolve)
        dependency_installed_lpsolve = true;
    end
    
case 'GLNXA'
    if ~dependency_installed_lpsolve
        if interactive
            affirmInstallation('lpsolve')
        end
        disp(['Installing ',...
              '<a href="http://lpsolve.sourceforge.net/">LPSOLVE 5.5</a>, ',...
              'Copyright (c) 2004 - 2021 by ',...
              'M. Berkelaar, K. Eikland, and P. Notebaert'])
        install_dir_lpsolve = fullfile(top_dir, 'dependencies', 'lpsolve');
        tmp_dir = tempname;
        mkdir(tmp_dir);
        lpsolve_ex = 'https://github.com/iqcToolbox/lpsolve552/raw/master/lp_solve_5.5.2.0_exe_ux64.tar.gz';
        lpsolve_dev = 'https://github.com/iqcToolbox/lpsolve552/raw/master/lp_solve_5.5.2.0_dev_ux64.tar.gz';
        lpsolve_mat = 'https://github.com/iqcToolbox/lpsolve552/raw/master/lp_solve_5.5.2.0_MATLAB_exe_ux64.tar.gz';
        fname_ex = fullfile(tmp_dir, 'lpsolve_ex.tar.gz');
        fname_dev = fullfile(tmp_dir, 'lpsolve_dev.tar.gz');
        fname_mat = fullfile(tmp_dir, 'lpsolve_mat.tar.gz');
        opts = weboptions('Timeout', 10);
        while ~exist(fname_ex, 'file')
            websave(fname_ex, lpsolve_ex, opts);
        end
        while ~exist(fname_dev, 'file')
            websave(fname_dev, lpsolve_dev, opts);
        end
        while ~exist(fname_mat, 'file')
            websave(fname_mat, lpsolve_mat, opts);
        end
        gunzip(fullfile(tmp_dir, 'lpsolve_ex.tar.gz'), tmp_dir);
        gunzip(fullfile(tmp_dir, 'lpsolve_dev.tar.gz'), tmp_dir);
        gunzip(fullfile(tmp_dir, 'lpsolve_mat.tar.gz'), tmp_dir);
        untar(fullfile(tmp_dir, 'lpsolve_ex.tar'), tmp_dir);
        untar(fullfile(tmp_dir, 'lpsolve_dev.tar'), tmp_dir);
        untar(fullfile(tmp_dir, 'lpsolve_mat.tar'), tmp_dir);
        % Need to rename mxlpsolve file to interface correctly with YALMIP
%         movefile(fullfile(tmp_dir, 'mxlpsolve.mexa64'),...
%              fullfile(tmp_dir, 'lp_solve.mexa64'))
        movefile(fullfile(tmp_dir), install_dir_lpsolve)
        dependency_installed_lpsolve = true;
    end
end

if ~dependency_installed_sdpt3
    try
        if interactive
            affirmInstallation('SDPT3')
        end
        disp(['Installing ',...
              '<a href="https://github.com/Kim-ChuanToh/SDPT3">SDPT3</a>, ',...
              'Copyright (c) 1997 by Kim-Chuan Toh, '...
              'Michael J. Todd, and Reha H. Tutuncu'])
        install_dir_sdpt3 = fullfile(top_dir, 'dependencies', 'sdpt3');
        tmp_dir = tempname;
        mkdir(tmp_dir)
        websave(fullfile(tmp_dir, 'SDPT3.zip'),...
                'https://github.com/Kim-ChuanToh/SDPT3/zipball/master');
        unzip(fullfile(tmp_dir, 'SDPT3.zip'), tmp_dir);
        dir_res = dir(tmp_dir);
        index = cell2mat({dir_res.isdir}) & ...
                ~strcmp({dir_res.name}, '.') & ... 
                ~strcmp({dir_res.name}, '..');
        tmp_move = fullfile(dir_res(index).folder, dir_res(index).name);
        movefile(tmp_move, install_dir_sdpt3)
        curr_dir = pwd;
        cd(install_dir_sdpt3)
        Installmex(1)
        cd(curr_dir)
        dependency_installed_sdpt3 = true;
    catch
        install_dir_sdpt3 = [];
    end
end
   
%% Save configuration
install_complete = dependency_installed_yalmip && ...
                   dependency_installed_sdpt3 && ...
                   dependency_installed_lpsolve;
mkdir(fullfile(top_dir, 'configuration'));
save(cfg_fname, 'install_complete',...
                'install_dir_sdpt3',...
                'install_dir_yalmip',...
                'install_dir_lpsolve',...
                'save_added_paths');
            
%% Run initialization
initializeIqcToolbox

function affirmInstallation(dependency)
% Variables for user interaction ([y/n])
y = 'y';
Y = y;
n = 'n';
N = n;
affirm = false;
while ~affirm
    yesno = input([dependency, ' is not detected on your system and ',...
                  'must be installed. Do you wish to install? [''y''/''n'']',...
                  '\n ']);
    affirm = strcmp(yesno, 'y');
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Added interactive behavior and installation of lpsolve - Micah Fry
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)