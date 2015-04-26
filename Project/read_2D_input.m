%---------------------------------------------------------------
%
%  read_2D_input.m
%  ------------
%
%  Reads in input file for 2D model problem (input.2d)
%
%  Variables assigned in this code:
%
%  probid          = character string identifying problem
%  KofXY(nelems,1) = material function k(x,y) evaluated at the 
%                    at the centroid of each element
%  BofXY(nelems,1) = material function b(x,y) evaluated at the 
%                    at the centroid of each element
%  FofXY(nelems,1) = forcing function f(x,y) evaluated at the 
%                    at the centroid of each element
%  VBC1(nbc1,1)    = values of Dirichlet (type 1) boundary data
%  VBC2(nbc2,1)    = values of Neumann (type 2) boundary data   
%
%---------------------------------------------------------------

% Open input file

fid = fopen('input.2d');

% Read in mesh name

probid = fgetl(fid);

% Read in material functions and forcing term values

KofXY = zeros(nelems,1);
BofXY = zeros(nelems,1);
FofXY = zeros(nelems,1);
for j = 1:nelems
    jj = fscanf(fid,'%g',1);
    if (jj~=j)
        fclose('all')
        error('********** Element numbering in input file is not sequential **********')
    end
    KofXY(j) = fscanf(fid,'%g',1);
    BofXY(j) = fscanf(fid,'%g',1);
    FofXY(j) = fscanf(fid,'%g',1);    
end

% Read in boundary condition data

VBC1 = zeros(nbc1,1);
VBC2 = zeros(nbc2,1);
for i = 1:nbc1
    nbc1_temp = fscanf(fid,'%g',1);    
    if (nbc1_temp ~= NODEBC1(i))
        fclose('all')
        error('********** Your mesh and input file boundary data are not consistent **********')        
    else
        NODEBC1(i) = nbc1_temp;
    end        
    VBC1(i) = fscanf(fid,'%g',1);
end
for i =1:nbc2
    nbc2_temp = fscanf(fid,'%g',1);    
    if (nbc2_temp ~= NODEBC2(i))
        fclose('all')
        error('********** Your mesh and input file boundary data are not consistent **********')        
    else
        NODEBC2(i) = nbc2_temp;
    end   
    VBC2(i) = fscanf(fid,'%g',1);
end

% Close the input file

fclose(fid);

% Clear out superfluous variables

clear ans fid i j jj n;