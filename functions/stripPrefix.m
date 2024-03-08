function [folder, prefix] = stripPrefix(dirNames)
% STRIPPREFIX Find and remove the common prefix from a set of directory names
%   [folders, prefix] = STRIPPREFIX(folderNames) returns the
%   common prefix of a set of directory names, and removes it from each
%   folder name. 
% 
%   Input:
%       dirNames: Cell array of strings containing the names of
%       directories
%   Output:
%       dir: Cell array of subdirectory names with the common prefix being
%       removed.
%       prefix: String of the common prefix for all dirs.
%
%   Falk Oraschewski
%   04.04.2023

% Split folder names into parts
parts = cellfun(@(x) strsplit(x, filesep), dirNames, 'UniformOutput', false);
parts = vertcat(parts{:});

nParts = size(parts,2); % Number of different parts

% Find the common prefix of all the parts
commonPrefix = '';
for i = 1:nParts
    if all(strcmp(parts(1,i), parts(2:end,i)))
        commonPrefix = [commonPrefix parts{1,i} filesep];
    else
        break
    end
end

% Remove common prefix from each folder name
folder = cellfun(@(x) strrep(x, commonPrefix, ''), dirNames, 'UniformOutput', false);

% Remove any trailing file separator characters from the prefix
prefix = strip(commonPrefix, 'right', filesep);
end