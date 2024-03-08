function [Rcoarse,varargout] = fmcw_range_crop(maxrange,Rcoarse,varargin)

% Rcoarse,Rfine,spec_cor,spec] = fmcw_range_crop(maxrange,Rcoarse,Rfine,spec_cor,spec)
%
% Falk Oraschewski
% 17 March 2022

n = find(Rcoarse<=maxrange,1,'last');
Rcoarse = Rcoarse(1:n);

for k = 1:numel(varargin)
    varargout{k} = varargin{k}(:,1:n);
end
