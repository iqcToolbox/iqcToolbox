function delta_out = toDelta(rct_obj, varargin)
%% TODELTA - Function for converting uncertain RCT objects to iqcToolbox Deltas.
%
%     delta_out = toDelta(rct_obj, varagin)
%
%     If input is ultidyn, output is a DeltaDlti.
%     If input is ureal, output is a DeltaSlti. 
%     Note: rct_obj.NominalValue is lost when converting to Delta objects.
% 
%     Variables:
%     ---------
%       Input:
%         rct_obj : ureal or ultidyn (from Robust Control Toolbox)
%         varagin : map of delta names and objects from lftToRct
%       output:
%         delta_out : Delta object :: the resultant delta
%
%     See also rctToLft, lftToRct

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Calling appropriate functions for different input arguments

if nargin == 1
   
    if isa(rct_obj, 'ureal')
        delta_out = urealToDelta(rct_obj);
    elseif isa(rct_obj, 'ultidyn')
        delta_out = ultidynToDelta(rct_obj);
    else
        error('toDelta:toDelta', ['When providing a single argument for toDelta, ',...
              'that input must be a ureal or ultidyn']);
    end
elseif nargin == 2
    delta_out = normalizeLft(toLft(varargin{1}(erase(rct_obj.Name, 'Normalized')))).delta.deltas{1};
end

end

function delta_out = ultidynToDelta(ulti_in)
    goodDelta = @(delta) isa(delta, 'ultidyn') &&...
                         isequal(delta.Type, 'GainBounded');
    assert(goodDelta(ulti_in),...
           'toDelta:toDelta',...
           strcat('Cannot make a DeltaDlti obect where uncertainties are not',...
                  ' GainBounded ultidyns'))
    delta_out = DeltaDlti(ulti_in.Name,...
                          size(ulti_in.OutputUnit,1),...
                          size(ulti_in.InputUnit,1),...
                          ulti_in.Bound);
end

function delta_out = urealToDelta(ureal_in)
    nom_check = @(delta) isequal(delta.NominalValue, 0);
    
    if ~nom_check(ureal_in)
       warning('toDelta:toDelta',...
                strcat('Nominal value information will not be retained',...
                       ' in the DeltaSlti'));
    end
    
    delta_out = DeltaSlti(ureal_in.Name,...
                          1,...
                          ureal_in.Range(1),...
                          ureal_in.Range(2));
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Philip Mulford (philip.mulford@ll.mit.edu)