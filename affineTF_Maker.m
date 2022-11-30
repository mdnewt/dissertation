function [TForm,TF] = affineTF_Maker(T,S,R)
% AffineTF_Maker accepts translation, scaling, and rotational components as
% vectors and combined them into an overall affine transformation matrix.
% Inputs must be given as rowwise vectors (use ',' not ';' to delimit)
%
% Assumes that initial rotations are done in this order:
% (axial,sagittal,coronal). The input angles should also be in this order,
% assuming the original stack dimensions are (coronal x sagittal x axial).
% 
% Inputs: T - Translation vector of length 3 for 3D (X,Y,Z) or 2 for 2D (X,Y)
%         S - Scaling vector of length 3 for 3D (Sx,Sy,Sz) or 2 for 2D (Sx,Sy)
%         R - For 3D, a rotation vector in degrees describing rotation 
%             around the z, x, and y axes (for 3D). For 2D, a rotation
%             scalar in degrees.
% 
% Outputs: TForm - The affine transform object which can be applied
%                  directly using imwarp
%             TF - The transformation matrix describing the combined
%                  effects of translation, rotation, and scaling

%% Check that vectors are the correct length
if length(T) == 3
    if (length(S) == 3 && length(R) == 3) ~= 1
        error(['Incorrect inputs. For a 3D transformation, T, S, and ',...
            'R are of length 3. For a 2D transformation, T and S are ',...
            'of length 2 and R is of length 1. See function notes for',...
            'a more comprehensive description.']);
    end
    d = 3;
elseif length(T) == 2
    if (length(S) == 2 && length(R) == 1) ~= 1
        error(['Incorrect inputs. For a 3D transformation, T, S, and ',...
            'R are of length 3. For a 2D transformation, T and S are ',...
            'of length 2 and R is of length 1. See function notes for',...
            'a more comprehensive description.']);
    end
    d = 2;
else
    error(['Incorrect inputs. For a 3D transformation, T, S, and ',...
        'R are of length 3. For a 2D transformation, T and S are ',...
        'of length 2 and R is of length 1. See function notes for',...
        'a more comprehensive description.']);
end

%% Build separate transformation matrices

% Creat base identity matrix
base = diag(ones(d+1,1));

% Translation and scaling
t  = base; t(end    ,1:end-1) = T      ; % Transformation
s  = base; s(1:end-1,1:end-1) = diag(S); % Scaling

% Rotation
switch d
    case 3
        rz = base; rz(1:2,1:2) = [cosd(R(1)),-sind(R(1)) ;...
                                  sind(R(1)), cosd(R(1))];
        rx = base; rx(2:3,2:3) = [cosd(R(2)),-sind(R(2)) ;...
                                  sind(R(2)), cosd(R(2))];
        ry = base; ry(1:3,1:3) = [cosd(R(3)),0,-sind(R(3)) ;...
                                  0         ,1, 0          ;...
                                  sind(R(3)),0, cosd(R(3))];
        r = rz*ry*rx;
    case 2
        r = base; r(1:2,1:2) = [cosd(R(1)),-sind(R(1)) ;...
                                sind(R(1)), cosd(R(1))];
end
        
% Combine into one transformation matrix
TF = r*s*t;

% Format as an affine transform object
switch d
    case 3, TForm = affine3d(TF);
    case 2, TForm = affine2d(TF);
end