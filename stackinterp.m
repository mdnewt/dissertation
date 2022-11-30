function stackOut = stackinterp(varargin)
% STACKINTERP accepts the following inputs:
% 1) stackinterp(stackIn,[sizeOut]); where [sizeOut] is a vector describing
%    the size of the output stack.
% 2) stackinterp(stackIn,[sizeOut],'method'); where 'method' is either
%    'interp3' or 'tform'. Default is 'interp3'.

stackIn = varargin{1}   ;
sizeOut = varargin{2}   ;
sizeIn  = size( stackIn);
classIn = class(stackIn);
if ~any(strcmp(classIn,{'single','double'})), stackIn = single(stackIn); end

if nargin == 2
   [x,y,z] = meshgrid(linspace(1,sizeIn(2),sizeOut(2)),...
                      linspace(1,sizeIn(1),sizeOut(1)),...
                      linspace(1,sizeIn(3),sizeOut(3)));
    stackOut = interp3(stackIn,x,y,z);
elseif nargin ==3
    method = varargin{3};
    switch method
        case 'interp3'
           [x,y,z] = meshgrid(linspace(1,sizeIn(2),sizeOut(2)),...
                              linspace(1,sizeIn(1),sizeOut(1)),...
                              linspace(1,sizeIn(3),sizeOut(3)));
            stackOut = interp3(stackIn,x,y,z);
        case 'tform'
            dr = sizeOut/sizeIn;
            T  = maketform('affine',[dr 0 0;0 dr 0;0 0 dr;0 0 0]);
            R  = makeresampler({'linear','linear','linear'},'fill');
            stackOut = round(tformarray(stackIn,T,R,[1,2,3],[1,2,3],sizeOut,[],0));
    end
elseif nargin ==4
    method = varargin{3};
    interp = varargin{4};
    switch method
        case 'interp3'
           [x,y,z] = meshgrid(linspace(1,sizeIn(2),sizeOut(2)),...
                              linspace(1,sizeIn(1),sizeOut(1)),...
                              linspace(1,sizeIn(3),sizeOut(3)));
            stackOut = interp3(stackIn,x,y,z,interp);
        case 'tform'
            dr = sizeOut/sizeIn;
            T  = maketform('affine',[dr 0 0;0 dr 0;0 0 dr;0 0 0]);
            R  = makeresampler({interp,interp,interp},'fill');
            stackOut = round(tformarray(stackIn,T,R,[1,2,3],[1,2,3],sizeOut,[],0));
    end
end

if ~any(strcmp(classIn,{'single','double'}))
    eval(['stackOut = ',classIn,'(stackOut);']);
end

end