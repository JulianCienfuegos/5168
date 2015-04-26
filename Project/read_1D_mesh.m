%-----------------------------------------------------------
%
%  read_mesh.m
%  ------------
%
%  Reads in the finite element mesh file (mesh.1d).
%
%  Variables assigned in this code:
%
%  mesh_name        = character string of mesh name
%  nelems           = number of elements in mesh
%  nnodes           = number of nodes in mesh
%  XNODES(nnodes,1) = x-coordinates of mesh nodes
%  CONN(nelems,2)   = element to node connectivity
%  nbc1             = number of Dirichlet (type 1) boundary nodes
%  nbc2             = number of Neumann (type 2) boundary nodes
%  NODEBC1(2,1)     = Dirichlet boundary node numbers 
%  NODEBC2(2,1)     = Neumann boundary node numbers
%
%-----------------------------------------------------------

% Open mesh file

fid = fopen('mesh.1d');

% Read in mesh name

mesh_name = fgetl(fid);

% Read in number of elements and number of nodes

nelems = fscanf(fid,'%g',1);
nnodes = fscanf(fid,'%g',1);

% Compute the degree of polynomial

p = (nnodes-1)/nelems;

% Read in element end-point coordinates

XNODES = zeros(nnodes,1);
ii = fscanf(fid,'%g',1);
XNODES(ii) = fscanf(fid,'%g',1);
for i = 1:nelems  
    jj = fscanf(fid,'%g',1);
    if ((ii + p)~=jj)
        fclose('all')
        error('********** Node numbering in mesh file is not sequential **********')
    else
        XNODES(ii+p) = fscanf(fid,'%g',1); 
        for j = 1:p-1
            XNODES(ii+j) = XNODES(ii) + j/p*(XNODES(ii+p)-XNODES(ii));
        end  
        ii = ii + p;        
    end
end

% Read in element connectivity table

CONN = zeros(nelems,2);
for j = 1:nelems
    jj = fscanf(fid,'%g',1);
    if (jj~=j)
        fclose('all')
        error('********** Element numbering in mesh file is not sequential **********')
    end
    CONN(j,1:2) = fscanf(fid,'%g %g',[1 2]);
end

% Read in the Dirichlet (type 1) boundary condition data

nbc1 = fscanf(fid,'%g',1);
if ( nbc1<0 | nbc1>2 )
    fclose('all')
    error('********** Number of boundary nodes must be 0, 1, or 2 **********')    
end
NODEBC1 = zeros(2,1);
for i = 1:nbc1 
    NODEBC1(i) = fscanf(fid,'%g',1);    
end

% Read in the Neumann (type 2) boundary condition data

nbc2 = fscanf(fid,'%g',1);
if ( (nbc1+nbc2)~=2 )
    fclose('all')
    error('********** Total number of boundary nodes must equal 2 **********')    
end
NODEBC2 = zeros(2,1);
for i = 1:nbc2
    NODEBC2(i) = fscanf(fid,'%g',1);
end

% Close the mesh file

fclose(fid);
