function U = exact_1D(X,probid)

%-----------------------------------------------------------------
%
%  Input:
%  ------
%       
%  X      = array of points at which to evaluate the exact solution
%  probid = character string identifying problem (specified in 
%           input.1d) 
%
%  Output:
%  -------
%
%  U = exact solution evaluted at points X    
%
%-----------------------------------------------------------------

if     (probid == 'Problem 1.1a ')    
    U = log(1+X)/log(2) - X;         
elseif (probid == 'Problem 1.1b ')
    U = X - sinh(X)/cosh(1);         
elseif (probid == 'Problem 1.1c ')
    U = 5 + (1-X).*(atan(50*X-75/2) + atan(75/2));
elseif (probid == 'Problem 1.1d ')
    U = -X.^4 + 2*X + 1;
elseif (probid == 'Problem 1.2  ')    
    U = 0 + (100 - 0).*log(10./X)/log(10/1);  
else
   error('********** Invalid Problem ID!! **********')
end