function [uv] = rot(ang,uv)
    
    % Rotation matrix and centroid
    R     = [cosd(ang),-sind(ang);sind(ang),cosd(ang)];
    cent  = [mean(uv(:,1)),mean(uv(:,2))]             ;

    % uv
    cent2 = repmat(cent,length(uv),1);
    temp  = (uv(:,1:2)-cent2)*R+cent2;
    uv(:,1:2) = temp                 ;
    uv(:,1) = uv(:,1)-min(uv(:,1))   ;
    uv(:,2) = uv(:,2)-min(uv(:,2))   ;

end