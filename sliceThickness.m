function [sliceThickness,info] = sliceThickness()
% Gives slice thickness of DCM file

[f,p] = uigetfile('*.DCM');
info = dicominfo([p,f]);
sliceThickness = info.SliceThickness;


end

