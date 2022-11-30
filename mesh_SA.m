function area = mesh_SA(Meshes,vox)

area = zeros(length(Meshes),1);
for ii = 1:length(Meshes)

    % Get 'verts' and 'faces'
    verts = Meshes(ii).v;
    faces = Meshes(ii).f;
   
    % Run operations on verts and faces
    a = verts(faces(:, 2), :) - verts(faces(:, 1), :);
    b = verts(faces(:, 3), :) - verts(faces(:, 1), :);
    c = cross(a, b, 2);
    area(ii) = 1/2 * sum(sqrt(sum(c.^2, 2)));
    
end

if exist('vox','var'), area = area .* vox ^2; end