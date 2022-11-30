function str = strcrop(str,delim,pos)
% STRCROP leverages strsplit to isolate just one part of a string, but all
% in one line of code. Also works on cells of strings.
% 
% Inputs:   str - The string or cell of strings. Must be a 'char' or a cell
%                 containing all 'char' data.
% 
%         delim - The delimiter (or a cell of multiple delimiters).
%                 Default value '.'
% 
%           pos - The position of the split string you want to keep.
%                 Default value [1].

% Check inputs
if ~(ischar(str) || (iscell(str) && all(cellfun(@ischar,str))))
    error('''str'' must be a char or cell of chars');
end
if ~exist('delim','var'), delim = '.'; end
if ~exist('pos'  ,'var'), pos   =  1 ; end

% Define in-line functions
if iscell(str)
      fa = @(str) strsplit(str,delim); fb = @(str) str{pos};
      f1 = @(str) cellfun(fa,str,'UniformOutput',false);
      f2 = @(str) cellfun(fb,str,'UniformOutput',false);
else, f1 = @(str) strsplit(str,delim); f2 = @(str) str{pos};
end

% Perform the crop
str = f2(f1(str));