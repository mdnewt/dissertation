function [mask,thresh] = otsu(image)

thresh = graythresh(image       );
mask   = imbinarize(image,thresh);