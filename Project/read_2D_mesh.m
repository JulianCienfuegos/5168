%-----------------------------------------------------------
%
%  read_2D_mesh.m
%  ------------
%
%  Reads in the finite element mesh file (mesh.2d).
%
%  Variables assigned in this code:
%
%  mesh_name        = character string of mesh name
%  nelems           = number of elements in the mesh
%  nnodes           = number of nodes in the mesh
%  XNODES(nnodes,1) = node x-coordinates
%  YNODES(nnodes,1) = node y-coordinates
%  CONN(nelems,3)   = element to node connectivity
%  nbc1             = number of Dirichlet (type 1) boundary nodes
%  nbc2             = number of Neumann (type 2) boundary nodes
%  NODEBC1(nbc1,1)  = Dirichlet boundary node numbers 
%  NODEBC2(nbc2,1)  = Neumann boundary node numbers
%
%-----------------------------------------------------------

% Open mesh file

fid = fopen('mesh.2d');

% Read in mesh name

mesh_name = fgetl(fid);

% Read in number of elements and number of (end-point) nodes

nelems = fscanf(fid,'%g',1);
nnodes = fscanf(fid,'%g',1);

% Read in nodal coordinates

XNODES = zeros(nnodes,1);
YNODES = zeros(nnodes,1);
for i = 1:nnodes
    ii    = fscanf(fid,'%g',1);
    if (ii~=i)
        fclose('all')
        error('********** Node numbering in mesh file is not sequential **********')
    end
    XNODES(i) = fscanf(fid,'%g',1);
    YNODES(i) = fscanf(fid,'%g',1);
    Z(i) = fscanf(fid,'%g',1);
end

% Read in element connectivity table

CONN = zeros(nelems,3);
for j = 1:nelems
    jj = fscanf(fid,'%g',1);
    if (jj~=j)
        fclose('all')
        error('********** Element numbering in mesh file is not sequential **********')
    end
    n = fscanf(fid,'%g',1);
    for i = 1:n
        CONN(j,i) = fscanf(fid,'%g',1);
    end
end

% Read in the Dirichlet (type 1) boundary condition data

nbc1 = fscanf(fid,'%g',1);
NODEBC1 = zeros(nbc1,1);
for i = 1:nbc1
    NODEBC1(i) = fscanf(fid,'%g',1);        
end

% Read in the Neumann (type 2) boundary condition data

nbc2 = fscanf(fid,'%g',1);
NODEBC2 = zeros(nbc2,1);
for i = 1:nbc2
    NODEBC2(i) = fscanf(fid,'%g',1);        
end

% Close the mesh file

fclose(fid);

% Clear out superfluous variables

clear ans fid i ii j jj n;