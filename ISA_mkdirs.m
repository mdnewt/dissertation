function ISA_mkdirs(p,z,v)
% ISA_MKDIRS creates the required directory structure to run the iterative
% shape averaging (ISA) algorithm.
%  - 'p' is the base directory where the new directories will be placed.
%  - 'z' is a logical input; if z = [true]. The optional '0 - Pre-training'
%    directory is created. Default value is [false].
%  - 'v' is a logical input; if v = [true]. The optional validation sub-
%    directories are created. Default value is [false].
% 
% This function is part of the iterative shape averaging (ISA) toolbox.
% Version history: V1 - 2017 May 02.
%                  V2 - 2018 Oct 24. - Removed 'Training' and 'Unknown
%                       Samples' folders - these are no longer necessary.

if ~exist('p','var') || isempty(p), p = uigetdir; end

% Create directories
d = [p,'\1 - Average Atlas'           ]; if ~exist(d,'dir'), mkdir(d); end
d = [p,'\1 - Average Atlas\Iterations']; if ~exist(d,'dir'), mkdir(d); end
d = [p,'\2 - Registrations'           ]; if ~exist(d,'dir'), mkdir(d); end

% Create optional pre-training directories
if exist('z','var') && z
    d = [p,'\0 - Pre-training'       ]; if ~exist(d,'dir'), mkdir(d); end
    d = [p,'\0 - Pre-training\Stacks']; if ~exist(d,'dir'), mkdir(d); end
    d = [p,'\0 - Pre-training\Masks' ]; if ~exist(d,'dir'), mkdir(d); end
end

% Create optional LOO cross-validation directories
if exist('v','var') && v
    d = [p,'\1 - Average Atlas\Validation\Iterations\']; if ~exist(d,'dir'), mkdir(d); end
end