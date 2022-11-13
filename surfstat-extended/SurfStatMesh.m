function surf = SurfStatMesh( mesh);

%Converts MATLAB mesh format to SurfStatMesh format
% 
% Usage: surf = SurfStatMesh( mesh);
% 
% mesh           = Matlab surface mesh format given as a a structured array:
%                  struct with fields:
%                         vertices: [#vertices×3 double]
%                            faces: [#triangles×3 double]
%
%  
% surf.coord    = 3 x v matrix of coordinates.
% surf.tri      = t x 3 matrix of triangle indices, 1-based, t=#triangles.

surf.coord = mesh.vertices';
surf.tri   = mesh.faces;
end
