function [tx,ty,tz] = align_centers(mov,rM,fix,rF,mode)
% 'mode' can have two values: 'com' - center of mass [default]
%                             'geo' - geometric

if ~exist('mode','var'), mode = 'com'; end
switch mode
case 'com'
    % Align centers
    sumFixedIntensity  = sum(fix(:));
    sumMovingIntensity = sum(mov(:));
    szF = size(fix); szM = size(mov);
    if     isa(rM,'imref2d') && isa(rF,'imref2d')
        dim = 2;
       [xFixed ,yFixed ] = meshgrid(1:szF(2),1:szF(1));
       [xMoving,yMoving] = meshgrid(1:szM(2),1:szM(1));
    elseif isa(rM,'imref3d') && isa(rF,'imref3d')
        dim = 3;
       [xFixed ,yFixed ,zFixed ] = meshgrid(1:szF(2),1:szF(1),1:szF(3));
       [xMoving,yMoving,zMoving] = meshgrid(1:szM(2),1:szM(1),1:szM(3));
    else,  error('Inconsistency in image and/or reference frame dimensionality.');
    end

    fixedXCOM  = (rF.PixelExtentInWorldX .* (sum(xFixed( :).*double(fix(:)))...
               ./ sumFixedIntensity)) + rF.XWorldLimits(1);
    fixedYCOM  = (rF.PixelExtentInWorldY .* (sum(yFixed( :).*double(fix(:)))...
               ./ sumFixedIntensity)) + rF.YWorldLimits(1);
    movingXCOM = (rM.PixelExtentInWorldX .* (sum(xMoving(:).*double(mov(:)))...
               ./ sumMovingIntensity)) + rM.XWorldLimits(1);
    movingYCOM = (rM.PixelExtentInWorldY .* (sum(yMoving(:).*double(mov(:)))...
               ./ sumMovingIntensity)) + rM.YWorldLimits(1);
    tx = fixedXCOM - movingXCOM;
    ty = fixedYCOM - movingYCOM;
    if dim == 3
        fixedZCOM  = (rF.PixelExtentInWorldZ .* (sum(xFixed( :).*double(fix(:)))...
                   ./ sumFixedIntensity)) + rF.ZWorldLimits(1);
        movingZCOM = (rM.PixelExtentInWorldZ .* (sum(yMoving(:).*double(mov(:)))...
                   ./ sumMovingIntensity)) + rM.ZWorldLimits(1);
        tz = fixedZCOM - movingZCOM;
    end
case 'geo'
    if     isa(rM,'imref2d') && isa(rF,'imref2d'), dim = 2;
    elseif isa(rM,'imref3d') && isa(rF,'imref3d'), dim = 3;
    else,  error('Inconsistency in image and/or reference frame dimensionality.');
    end
    % Align centers
    fixedCenterXWorld  = mean(rF.XWorldLimits);
    fixedCenterYWorld  = mean(rF.YWorldLimits);
    movingCenterXWorld = mean(rM.XWorldLimits);
    movingCenterYWorld = mean(rM.YWorldLimits);
    tx = fixedCenterXWorld - movingCenterXWorld;
    ty = fixedCenterYWorld - movingCenterYWorld;
    if dim == 3
        fixedCenterZWorld  = mean(rF.ZWorldLimits);
        movingCenterZWorld = mean(rM.ZWorldLimits);
        tz = fixedCenterZWorld - movingCenterZWorld;
    end
end