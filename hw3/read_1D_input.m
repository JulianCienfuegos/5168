%---------------------------------------------------------------
%
%  read_input.m
%  ------------
%
%  Reads in input file for 1D model problem (input.1d)
%
%  Variables assigned in this code:
%
%  probid         = character string identifying problem
%  KofX(nelems,1) = material function k(x) evaluated at the 
%                   at the midpoint of each element
%  BofX(nelems,1) = material function b(x) evaluated at the 
%                   at the midpoint of each element
%  FofX(nelems,1) = forcing function f(x) evaluated at the 
%                   at the midpoint of each element
%  VBC1(nbc1,1)   = values of Dirichlet (type 1) boundary data
%  VBC2(nbc2,1)   = values of Neumann (type 2) boundary data   
%
%---------------------------------------------------------------

% Open input file

fid = fopen('input.1d');

% Read in mesh name

probid = fgetl(fid);

% Read in material functions and forcing term values

KofX = zeros(nelems,1);
BofX = zeros(nelems,1);
FofX = zeros(nelems,1);
for j = 1:nelems
    jj = fscanf(fid,'%g',1);
    if (jj~=j)
        fclose('all')
        error('********** Element numbering in input file is not sequential **********')
    end
    KofX(j) = fscanf(fid,'%g',1);
    BofX(j) = fscanf(fid,'%g',1);
    FofX(j) = fscanf(fid,'%g',1);    
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
