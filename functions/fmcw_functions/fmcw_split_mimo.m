function [varargout] = fmcw_split_mimo(vdat)
%FMCW_SPLIT_MIMO seperates multiple antenna contributions in vdat.
%   Files that are recorded in MIMO mode (multiple input, multiple output),
%   contain information from various Tx-Rx combinations. This function
%   disentagles them. The output is either a vdat-structure containing all
%   of the combinations as seperate elements or, if several outputs are
%   assigned, the individual vdats are given as seperate structures. In any
%   case, the vdat elements/structures will look like they were recorded as
%   a single non-mimo file.
%
%   input: original vdat
%   output: vdat-structure containing multiple elements for each antenna
%           combination OR several vdats for every antenna combination%
%
%   Falk Oraschewski,
%   2022 December 24, Neumayer Station III

nout = max(nargout,1); % Number of output arguments

% Get number of Tx-Rx combinations
n_Tx = numel(vdat.TxAnt);
n_Rx = numel(vdat.RxAnt);
n_comb = n_Tx * n_Rx; % Number of antenna combinations.

% 
comb_Tx = repelem((1:n_Tx),1,n_Rx);
comb_Rx = repmat((1:n_Rx),1,n_Tx);

vdat_all = repmat(vdat,1,n_comb);

% Disentangle vdat
for comb = 1:n_comb
    n_chirp = vdat.ChirpsInBurst/n_comb;
    vdat_all(comb).ChirpsInBurst = n_chirp;
    vdat_all(comb).TxAnt = vdat.TxAnt * 0;
    vdat_all(comb).TxAnt(comb_Tx(comb)) = 1;
    vdat_all(comb).RxAnt = vdat.RxAnt * 0;
    vdat_all(comb).RxAnt(comb_Rx(comb)) = 1;
    n_v = length(vdat.v)/n_comb;
    vdat_all(comb).v = vdat.v(n_v*(comb-1)+1:n_v*comb);
    vdat_all(comb).Startind = vdat.Startind(1:n_chirp);
    vdat_all(comb).Endind = vdat.Endind(1:n_chirp);
    vdat_all(comb).vif = vdat.vif(n_chirp*(comb-1)+1:n_chirp*comb,:);
    vdat_all(comb).chirpNum = vdat.chirpNum(1:n_chirp);
    vdat_all(comb).chirpAtt = vdat.chirpAtt(n_chirp*(comb-1)+1:n_chirp*comb);
    vdat_all(comb).chirpTime = vdat.chirpTime(n_chirp*(comb-1)+1:n_chirp*comb);
end

% Process output depending on number of output arguments.
if nout > n_comb
    error(['FMCW_SPLIT_MIMO: Too many output variables assigned. For the selected input-vdat a maximum of ',num2str(n_comb), ' output arguments is possible.'])
elseif nout > 1
    if nout < n_comb
        disp(['The file contains data from ',num2str(n_comb),' Tx-Rx-combinations, but only the first ',num2str(nout),' combinantions are loaded.'])
    end

    for k = 1:nout
        varargout{k} = vdat_all(k);
    end
else
    if n_comb > 1
        disp(['The vdat contains data from ',num2str(n_comb),' Tx-Rx-combinations which are split into different structure elements. To load them as individual vdats assign ',num2str(n_comb),' output variables.'])
    end
    varargout{1} = vdat_all;
end