classdef SequencePerformance
%% SEQUENCEPERFORMANCE class for collecting different Performance objects.
% This class assists in conducting the same operation to all the 
% Performance objects in its sequence
%
%   methods:
%     SequencePerformance(varargin) :: Constructor method
%     disp(this_performances) :: Display method
%     matchHorizonPeriod(this_performances, horizon_period) :: Matches performances properties to new horizon_period
%
%   properties:
%     performances : a cell array of Performance objects
%
%   dependent properties:
%     chan_ins : (num_del x 1) cell array of each performances chan_in property
%     chan_outs : (num_del x 1) cell array of each performances chan_out property
%     horizon_periods : (num_del x 2) array of performance horizon_periods
%     names : (1 x num_del) cell array of performance names
%     types : (1 x num_del) cell array of each performance's class
%
%   See also SequencePerformance.SequencePerformance, SequenceDelta, SequenceDisturbance, Ulft.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
    performances (1, :) cell
end

properties (Dependent)
    chan_outs
    chan_ins
    horizon_periods
    names
    types
end

methods
    function this_performances = SequencePerformance(varargin)
    %% SEQUENCEPERFORMANCE constructor
    %
    %  this_performances = SequencePerformance(dis1, dis2, dis3)
    %  this_performances = SequencePerformance(cell_of_diss)
    %  this_performances = SequencePerformance(cell_of_diss, dis_4, cell_of_diss_2)
    %
    %  Variables:
    %  ---------
    %     Input:
    %       any combination of cells of Performance objects and Performance objects    
    %     Output:
    %       this_performances : SequencePerformance object :: collection of given Performances
    %
    %  See also SequencePerformance.
    
        % Allocate empty cell of appropriate length
        len_all = sum(cellfun(@(arg) length(arg), varargin));
        performances = cell(1, len_all);
        
        % Check input arguments and fill
        isPerformance = @(perf) isa(perf, 'Performance');
        isCellOfPerformances = @(perf) iscell(perf) ...
                                       && ...
                                       all(cellfun(isPerformance, perf));
        ind_perf = 1;
        for i = 1:nargin
            perf = varargin{i};
            assert(isPerformance(perf) || isCellOfPerformances(perf),...
                   'SequencePerformance:SequencePerformance',...
                   ['Input objects must be Performances ',...
                    'and/or Cells of Performances'])
            len_perf = length(perf);
            if isPerformance(perf)
                performances{ind_perf} = perf;
            else
                performances(ind_perf : ind_perf + len_perf - 1) = horzcat(perf);
            end
            ind_perf = ind_perf + len_perf;
        end
        
        % Check if like-named performances are indeed equivalent
        names = cellfun(@(perf) perf.name, performances,'UniformOutput',false);
        [~, unique_indices] = unique(names);
        for i = unique_indices'
            perf = performances{i};
            perf_group = performances(strcmp(names{i}, names));
            for j = 1:length(perf_group)
                % Remove checks on chan_in, chan_out
                warn_state = warning('query', 'MATLAB:structOnObject');
                warning('off', 'MATLAB:structOnObject');
                perf_no_chan = rmfield(struct(perf), {'chan_in', 'chan_out'});
                perf_grp_no_dim = rmfield(struct(perf_group{j}),...
                                         {'chan_in', 'chan_out'});
                warning(warn_state);
                % Check all other properties
                assert(isequaln(perf_no_chan, perf_grp_no_dim),...                        
                       'SequencePerformance:SequencePerformance',...
                       ['There cannot be multiple Performances that share',...
                        'the same name but have different properties'])
            end
        end
        
        horizon_periods = cell2mat(cellfun(@(perf) perf.horizon_period,...
                                           performances',...
                                           'UniformOutput', false));
        horizon_period = commonHorizonPeriod(horizon_periods);
        performances = cellfun(@(perf) matchHorizonPeriod(perf, horizon_period),...
                                 performances,...
                                 'UniformOutput', false);
        
        % Assign to class property
        this_performances.performances = performances;
    end     
    
    function disp(this_performances)
    %% DISP function for SequencePerformance object
    %
    %  disp(this_performances) (e.g., disp(SequencePerformance(PerformanceL2('l2')) )
    %
    %  Variables:
    %  ---------
    %     Input:
    %       this_perf : SequencePerformance object  
    %
    %  See also Ulft.disp, Performance.disp
    
        fprintf([repmat('*', [1,80]), '\n',...
                 repmat('*', [1,31]),'   PERFORMANCES   ', repmat('*', [1,31]), '\n']);
        cellfun(@(perf) disp(perf), this_performances.performances);
    end
    
    function this_perf = matchHorizonPeriod(this_perf, horizon_period)
    %% MATCHHORIZONPERIOD function to apply matchHorizonPeriod to all Performances
    %  in SequencePerformance object
    %
    %  this_performances = matchHorizonPeriod(this_performances, total_time)
    %  Variables:
    %  ---------
    %     Input:
    %       this_performances : SequencePerformance object
    %       new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence   
    %     Output:
    %       this_performances : SequencePerformance object
    %
    %  See also Ulft.matchHorizonPeriod
    
        this_perf.performances = ...
            cellfun(@(perf, hp) matchHorizonPeriod(perf, hp),...
                    this_perf.performances,...
                    repmat({horizon_period}, size(this_perf.performances)),...
                    'UniformOutput', false);
    end
       
    %% Getter functions
    function chan_outs = get.chan_outs(this_performances)
        channel_cell = cellfun(@(perf) perf.chan_out,...
                               this_performances.performances',...
                               'UniformOutput', false);
        chan_outs = vertcat(channel_cell{:});
    end
    
    function chan_ins = get.chan_ins(this_performances)
        channel_cell = cellfun(@(perf) perf.chan_in,...
                               this_performances.performances',...
                               'UniformOutput', false);
        chan_ins = vertcat(channel_cell{:});
    end
       
    function names = get.names(this_performances)
        names = cellfun(@(perf) perf.name,...
                        this_performances.performances,...
                        'UniformOutput', false);
    end
    
    function types = get.types(this_performances)
        types = cellfun(@(perf) class(perf),...
                        this_performances.performances,...
                        'UniformOutput', false);
    end
    
    function horizon_period = get.horizon_periods(this_performances)
        horizon_period = cell2mat(cellfun(@(perf) perf.horizon_period,...
                                          this_performances.performances',...
                                          'UniformOutput', false));
    end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)