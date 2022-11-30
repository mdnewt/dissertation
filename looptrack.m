function looptrack(ii,L,t,name)
% LOOPTRACK.M outputs an indicator of the progress of a loop in the
% command line.
% 
% -Required inputs:
%     ii - Number of the current iteration
%      L - Total number of iterations
% -Optional inputs:
%      t - Iteration run time (execute t = tic at the beginning of the loop)
%   name - Name of the current sample
% 
% -Possible calls:
%   1) looptrack(ii,L)
%   2) looptrack(ii,L,t)
%   3) looptrack(ii,L,[],name)

if (exist('t','var') && ~isempty(t)) && (exist('name','var') && ~isempty(name))
    fprintf([num2str(ii),'/',num2str(L),' - ',name,' - ',num2str(toc(t)),' sec\n']);
elseif (exist('t','var') && ~isempty(t))
    fprintf([num2str(ii),'/',num2str(L),' - ',num2str(toc(t)),' sec\n']);
elseif (exist('name','var') && ~isempty(name))
    fprintf([num2str(ii),'/',num2str(L),' - ',name,'\n']);
else
    fprintf([num2str(ii),'/',num2str(L),'\n']);
end

