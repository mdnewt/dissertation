function [stack_re,TF,R_re,R,angles] = imstack_reorient(stack,order,offset,TF,skip,vox)
% IMSTACK_REORIENT performs reorientation on a 3D image stack. Orientations
% are determined by drawing alignment lines in all 3 planes. The rotation
% transform and reference frame are preserved so that the inverse transform
% can be applied later to map volumes-of-interest back to the original
% image. This function is intended for use with medical/biological images
% and thus the function is described with reference to the traditional
% anatomic planes (i.e. axial, coronal, and sagittal).
% 
% Inputs:  stack - 3D image stack (required). It is assumed that the image 
%                  stack is oriented as [coronal, sagittal, axial], and all
%                  function documentation refers to the planes in these
%                  terms.
% 
%          order - Optional input, given as a vector, describing the order
%                  in which sequential rotation is to be performed. Default
%                  value is [3,2,1] (i.e. [axial, sagittal, coronal]). If
%                  input 'TF' is supplied, input 'order' is unused and
%                  therefore can be given as an empty matrix.
% 
%         offset - Optional input, given as an integer vector of length 3,
%                  corresponding to the dimensionimstas of the image. Each
%                  value, if set to 1, adds a +90° rotation to the
%                  specified dimension. If set to -1, a -90° rotation is
%                  applied. If set to 0 [default], no rotation is applied.
% 
%             TF - Optional input, given as a valid geometric transform or
%                  displacement field compatible with the imwarp function.
%                  If specified, automatically applies this transform in 
%                  lieu of manual determination.
% 
%           skip - Optional input, given as a logical. If true, the actual
%                  reorientation is skipped and the TF is output.
% 
%            vox - Optional input. If given, used to create R or R_re.
% 
% Acceptable input formats: 1) imstack_reorient(stack);
%                           2) imstack_reorient(stack,order);
%                           3) imstack_reorient(stack,order,offset);
%                           4) imstack_reorient(stack,[],offset);
%                           5) imstack_reorient(stack,[],[],TF);
% 
% Input formats 1 and 2 will result in manual reorientation, while input
% format 3 will result in automated reorientation using the supplied
% transform or displacement field.
% 
% Outputs: stack_re - The reoriented output image stack.
% 
%                TF - The rotation matrix used to perform reorientation,
%                     stored as a MATLAB affine transformation object. The
%                     inverse transform will be used to backproject the
%                     reoriented reference frame onto the original
%                     reference frame.
% 
%              R_re - The spatial reference frame of the reoriented output
%                     stack.
% 
%                 R - The spatial reference frame of the input stack.
% 
%            angles - The angles of rotation in each plane. This output is
%                     optional, as it is not used in the imstack_derotate
%                     function.
% 
% Note: This function uses the JavaFrame reference to maximize the figure
%       window during assisted reorientation, and this produces a warning
%       message due to the eventual obsolescence of the JavaFrame property.
%       This warning message can be disabled, however we have not chosen to
%       do so automatically within the function.
% 
% Copyright 2017 - Michael Newton, Samantha Hartner, and Tristan Maerz.
% Contact: Michael Newton (E: mnewton@ltu.edu)

% Check to see if order is defined
if ~exist('order','var') || isempty(order)
    order = [3,2,1];
elseif length(order) ~=3 || ~all(ismember([1,2,3],order))
    order = [3,2,1];
    warning(['''order'' has been incorrectly specified. Proceeding ',...
        'with default value.']);
end

% Check to see if offset is defined
if ~exist('offset','var') || isempty(offset)
    offset = [0,0,0];
elseif length(offset) ~= 3 || max(abs(offset)) > 1
    offset = [0,0,0];
    warning(['''order'' has been incorrectly specified. Proceeding ',...
        'with default value.']);
end

% Check to see if skip is defined
if ~exist('skip','var') || isempty(skip)
    skip = false;
elseif length(skip) ~= 1 || max(abs(offset)) > 1
    skip = false;
    warning(['''skip'' has been incorrectly specified. Proceeding ',...
        'with default value (false).']);
end

% Check to see if vox is defined
if ~exist('vox','var') || isempty(skip)
    vox = 1;
elseif length(skip) ~= 1 || max(abs(offset)) > 1
    vox = 1;
    warning(['''vox'' has been incorrectly specified. Proceeding ',...
        'with default value 1.']);
end

%% Only run manual reorientation if TF is not specified as an input
if ~exist('TF','var') || isempty(TF)

% Create base identity matrix
base = diag(ones(4,1));

% Define function for determining reorientation angle while accounting for
% rotation of the displayed image
getang = @(pos) 90+(atan2(pos(1,2)-pos(2,2),pos(1,1)-pos(2,1))*180/pi);

%% Outer 'while' loop provides the ability to fix incorrect reorientations
w = 1;
while w

% Reset temporary stack and rotation matrix at beginning of each attempt
stack_temp = stack;
rotmat = base;

%% Inner 'for' loop through 3 sequential rotations in the order specified
for ii = 1:3
switch order(ii)
    
    case 1 % Coronal reorientation (dimension 1)
        
        % Permute
        % ** Note that the coronal permutation results in mirroring across
        %    the sagittal midline, necessitating the use of 'flip'
        stack_temp = permute(stack_temp,[2,3,1]);
        stack_temp = flip(   stack_temp,1      );
        
        % Determine reference slice
        C = round(size(stack_temp)./2);
        f = figure;
        imshow(imrotate([squeeze(stack_temp(C(1),:,:));...
                         squeeze(stack_temp(:,C(2),:))],-90),[]);
        title('Select reference slice (horizontal) for reorientation');
        g = ginput(1); g = round(g(2));
        
        % Draw alignment line
        imshow(imrotate(stack_temp(:,:,g),-90),[]);
        set(gcf,'WindowState','maximize');
        titleswitch(offset(ii));
        h = imline();
        pos = h.getPosition();
        close(f);
        
        % Calculate angle and perform reorientation
        corAngle = getang(pos) + (90 .* offset(ii));
%         corAngle = (atan2(pos(1,2)-pos(2,2),pos(1,1)-pos(2,1))*180/pi)+90;
        stack_temp = imrotate(stack_temp,corAngle);
        
        % Apply reorientation angle to overall rotation matrix
%         rx = base; rx(2:3,2:3) = [cosd(corAngle),-sind(corAngle) ;...
%                                   sind(corAngle), cosd(corAngle)];
%         rotmat = rotmat * rx;
        ry = base; ry(1:3,1:3) = [cosd(-corAngle) ,0,sind(-corAngle);...
                                  0              ,1,0             ;...
                                 -sind(-corAngle) ,0,cosd(-corAngle)];
        rotmat = rotmat * ry;
        
        % Display
        f = figure; imshow(imrotate(stack_temp(:,:,g),-90),[]);
        title('Any button to close'); waitforbuttonpress; close(f);
        
        % De-permute (and flip)
        stack_temp = flip(   stack_temp,1      );
        stack_temp = permute(stack_temp,[3,1,2]);
    
    case 2 % Sagittal reorientation (dimension 2)
        
        % Permute
        stack_temp = permute(stack_temp,[3,1,2]);
        
        % Determine reference slice
        C = round(size(stack_temp)./2);
        f = figure;
        imshow([squeeze(stack_temp(C(1),:,:));...
                squeeze(stack_temp(:,C(2),:))],[]);
        title('Select reference slice (vertical) for reorientation');
        g = ginput(1); g = round(g(1));
        
        % Draw alignment line
        imshow(stack_temp(:,:,g),[]);
        set(gcf,'WindowState','maximize');
        titleswitch(offset(ii));
        h = imline();
        pos = h.getPosition();
        close(f);
        
        % Calculate angle and perform reorientation
        sagAngle = getang(pos) + (90 .* offset(ii));
%         sagAngle = 90+(atan2(pos(1,2)-pos(2,2),pos(1,1)-pos(2,1))*180/pi);
        stack_temp = imrotate(stack_temp,sagAngle);
        
        % Apply reorientation angle to overall rotation matrix
        rx = base; rx(2:3,2:3) = [cosd(sagAngle),-sind(sagAngle) ;...
                                  sind(sagAngle), cosd(sagAngle)];
        rotmat = rotmat * rx;
%         ry = base; ry(1:3,1:3) = [cosd(sagAngle) ,0,sind(sagAngle);...
%                                   0              ,1,0             ;...
%                                  -sind(sagAngle) ,0,cosd(sagAngle)];
%         rotmat = rotmat * ry;
        
        % Display
        f = figure; imshow(stack_temp(:,:,g),[]);
        title('Any button to close'); waitforbuttonpress; close(f);
        
        % De-permute
        stack_temp = permute(stack_temp,[2,3,1]);
    
    case 3 % Axial reorientation (dimension 3)
        
        % Determine reference slice
        C = round(size(stack_temp)./2);
        f = figure;
        imshow([imrotate(squeeze(stack_temp(C(1),:,:)),-90),...
                imrotate(squeeze(stack_temp(:,C(2),:)),-90)],[]);
        title('Select horizontal reference slice for reorientation');
        g = ginput(1); g = round(g(2));
        
        % Draw alignment line
        imshow(stack_temp(:,:,g),[]);
        set(gcf,'WindowState','maximize');
        titleswitch(offset(ii));
        h = imline();
        pos = h.getPosition();
        close(f);
        
        % Calculate angle and perform reorientation
        axAngle = getang(pos) + (90 .* offset(ii));
%         axAngle = 90+(atan2(pos(1,2)-pos(2,2),pos(1,1)-pos(2,1))*180/pi);
        stack_temp = imrotate(stack_temp,axAngle);
        
        % Apply reorientation angle to overall rotation matrix
        rz = base; rz(1:2,1:2) = [cosd(axAngle),-sind(axAngle) ;...
                                  sind(axAngle), cosd(axAngle)];
        rotmat = rotmat * rz;
                            
        % Display
        f = figure; imshow(stack_temp(:,:,g),[]);
        title('Any button to close'); waitforbuttonpress; close(f);
        
end
end

% Show final orientation and confirm this is acceptable - if not, process
% will reset
C = round(size(stack_temp)./2);
f = figure; imshow([stack_temp(:,:,C(3)), squeeze(stack_temp(:,C(2),:));...
            imrotate(squeeze(stack_temp(C(1),:,:)),-90),...
            zeros(size(stack_temp,3),size(stack_temp,3))],[]);
title(['Is the current orientation acceptable? Mouse to proceed, key ',...
       'to try again.']);
w = waitforbuttonpress;
close(f);
        
end

clear stack_temp % Clear stack_temp to free memory

% Compute affine transformation
angles = [corAngle,sagAngle,axAngle];
TF = affine3d(rotmat);

else
    
angles = nan; % Placeholder value for output 'angles'

end

% Perform reorientation
if ~skip && exist('vox','var')
    R = imref3d(size(stack),vox,vox,vox);
   [stack_re,R_re] = imwarp(stack,R,TF);
elseif ~skip
    R = imref3d(size(stack),vox,vox,vox);
   [stack_re,R_re] = imwarp(stack,R,TF);
else stack_re = []; R_re = [];
end

%% titleswitch subfunction
function titleswitch(toggle)
    switch(toggle)
    case 0
    title(['Draw alignment line. In the plane shown, the image ',...
           'will be reoriented so that the line travels from the ',...
           'top to bottom of the image']);
    case 1
    title(['Draw alignment line. In the plane shown, the image ',...
           'will be reoriented so that the line travels from the ',...
           'left to right of the image']);
    case -1
    title(['Draw alignment line. In the plane shown, the image ',...
           'will be reoriented so that the line travels from the ',...
           'right to left of the image']);
    end
end

end