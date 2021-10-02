classdef DeltaIntegrator < Delta
%% DELTAINTEGRATOR class for the integrator operator Z, extends
% the base class Delta.
%
%   extended methods:
%     DeltaIntegrator(dim_state, timestep, horizon_period) :: Constructor
%     disp(this_delta) :: Display method
%     matchHorizonPeriod(this_delta, new_horizon_period) 
%                        :: Matches delta properties to new horizon_period
%     deltaToMultiplier(this_delta, varargin)
%                        :: Method for constructing a multiplier from a delta
%     normalizeDelta(this_delta) :: Changes unnormalized LFT to normalized one
%     toRct(this_delta) :: Converts Delta to object needed in lftToRct()
%
%   See also Delta, DeltaDelayZ

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

properties
end

methods
    function this_delta = DeltaIntegrator(dim_outin)
    %% DELTAINTEGRATOR constructor
    %
    %  d = DeltaIntegrator(dim_outin)
    %  d = DeltaIntegrator() assumes dim_outin == 1
    %
    %  Variables:
    %  ---------
    %    Input:
    %       dim_outin : natural :: output/input dimensions of integrator
    %    Output:
    %       this_delta : DeltaIntegrator object
    %
    %  See also Delta, Delta.Delta, DeltaIntegrator

        % Defining defaults for missing arguments
        switch nargin
            case 0
                dim_outin = 1;                    
        end
        horizon_period = [0, 1];
        validateattributes(dim_outin, 'numeric', {'numel', 1})
        % Calling Delta constructor
        this_delta@Delta('Integrator', dim_outin, dim_outin, horizon_period);                    
    end

    function disp(this_delta)
    %% DISP function for DeltaIntegrator object
    %
    %  disp(delta_integrator_obj) (e.g., disp(DeltaIntegrator())
    %
    %  Variables:
    %  ---------
    %     Input:
    %        this_delta : DeltaIntegrator object
    %
    %  See also DeltaIntegrator, Delta.disp, Ulft.disp

        disp@Delta(this_delta, 'integration operator')              
    end

    function this_delta = matchHorizonPeriod(this_delta, new_horizon_period)
    %% MATCHHORIZONPERIOD function to ensure properties of DeltaIntegrator 
    %  have the same length (which is no greater than 1)
    %
    %  matchHorizonPeriod(this_delta, new_horizon_period) will simply check that the desired horizon_period is admissible ([0, 1])
    %  matchHorizonPeriod(this_delta) checks that this_delta.horizon_period is admissible
    %
    %  Variables:
    %  ---------
    %     Input:
    %       this_del : DeltaIntegrator object
    %       new_horizon_period : 1 x 2 array of naturals :: [horizon, period] (in timesteps) of new sequence
    %    Output:
    %       this_delta : DeltaIntegrator object
    %
    %  Note, because the constructor forces the time length of
    %  dim_outin to be 1, this function merely returns an error if
    %  requested to extend property lengths beyond 1.
    %
    %  See also DeltaIntegrator
    
    if nargin == 1
    % Check correctness of this_delta.horizon_period
        assert(all(this_delta.horizon_period == [0, 1]),...
               ['The horizon_period for DeltaIntegrator is not [0,1]',...
                ' and is incompatible with the DeltaIntegrator class']);
    else
    % Check correctness of new_horizon_period
        assert(all(new_horizon_period == [0, 1]),...
               ['The new_horizon_period for DeltaIntegrator is not [0,1]',...
                ' and is incompatible with the DeltaIntegrator class']);
    end
    end

    function multiplier = deltaToMultiplier(this_del, varargin)                 %#ok<VANUS>
    %% DELTATOMULTIPLIER function to generate a multiplier from this object
    %  This function only exists to concretize the abstract class Delta.
    %  There is no multiplier for DeltaIntegrator, therefore this returns the
    %  default (or empty) Multiplier.
    %
    %  multiplier = deltaToMultiplier(this_del)
    %
    %  Variables:
    %  ---------
    %     Input:
    %       this_del : DeltaIntegrator object
    %     Output:
    %       multiplier : MultiplierDeltaDefault object
    %
    %  See also DeltaIntegrator
        multiplier = MultiplierDeltaDefault(this_del);
    end

    function [del_diff, del_ave, del_scale, del_norm] = normalizeDelta(this_delta)
    %% NORMALIZEDELTA function for DeltaIntegrator object. This method extendeds
    %  the default method, but because DeltaIntegrator objects are not
    %  normalized, the extended behavior essentially mimics that of the default
    %  behavior.
    %
    %    [del_diff, del_ave, del_scale, del_norm] = normalizeDelta(this_delta)
    %
    %    Variables:
    %    ---------
    %      Input:
    %         this_delta : Delta object
    %      Output:
    %         del_diff : double matrix :: Matrix for modifying portions of LFT a, b, c, d
    %                                       matrices pertaining to this Delta
    %         del_ave : double matrix :: Matrix for modifying portions of LFT a, b, c, d
    %                                       matrices pertaining to this Delta
    %         del_scale : double matrix :: Matrix for modifying portions of LFT a, b, c, d
    %                                       matrices pertaining to this Delta
    %         del_norm : Delta object (or empty) :: Normalized Delta object. If empty, indicates
    %                                      that no normalization took place and a warning is thrown.
    %
    %    See also Ulft.normalizeLft, Delta.normalizeDelta
    
    [del_diff, del_ave, del_scale] = normalizeDelta@Delta(this_delta);
    del_norm = this_delta;
    end
    
    function output = toRct(this_delta)                                         %#ok<MANU>
    %% TORCT function for DeltaIntegrator object. This extends the base method by 
    %  not returning any RCT object, because the RCT does not have separate
    %  uncertainties to express the integration operation.
        %
    %    rct_obj = toRct(this_delta)
    % 
    %    Variables:
    %    ---------
    %      Input:
    %        this_delta : Delta object to be converted
    %      Output:
    %        rct_obj : empty array
    %
    %    See also Delta.toRct, lftToRct, rctToLft
    output = [];
    end

    function value = sample(~, ~)                                               %#ok<STOUT>
        %% SAMPLE function for DeltaIntegrator, which does not represent an uncertainty.
        error('DeltaIntegrator:sample', 'DeltaIntegrator cannot be sampled.');
    end

    function validateSample(~, ~, ~)
        %% VALIDATESAMPLE function for DeltaIntegrator, which does not represent an uncertainty.
        error('DeltaIntegrator:validateSample',...
              'DeltaIntegrator cannot be sampled.');
    end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Added sample and validateSample methods - Jason Nezvadovitz (jason.nezvadovitz@ll.mit.edu)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)