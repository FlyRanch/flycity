function lsdata = legendreShapes(order,d,x) 
%
%  Returns order+1 Legendre shape functions or there derivatives at 
%  the requested points, x, which mush lie between -1 and +1. 
%  S.H. 1-06-09
%
%  order = maximum order of the shape functions
%  d     = order of derivative requested (0 = undifferentiated function)
%  h     = element length
%  x     = vector of requested points
%
%
   np     = max(size(x));
   lpoly  = zeros(order+1,np);    
   lsdata = zeros(order+1,np);    
%
%
%  --------------------------------------------------------------
%  Generate Legendre polynomial data using the recursive formula
%  --------------------------------------------------------------
   if order >= 0 
      lpoly(1,:) = 1.0;
   end 
   if order >= 1 
      lpoly(2,:) = x;
   end
   for n = 1:order-1
     ni = n + 1;
     for i = 1:np
       lpoly(ni+1,i) = ( (2*n+1)*lpoly(ni,i)*x(i) - n*lpoly(ni-1,i) ) / (n+1);
     end
   end
%
%
%  -------------------------------------------
%  Return shapes or derivatives, as requested
%  -------------------------------------------
%
   if d == 0

     lsdata = lpoly;

   elseif d == 1

    if order >= 0 
      lsdata(1,:) = 0.0;
    end
    for n = 0:order-1
      ni = n + 1;
      for i = 1:np
        lsdata(ni+1,i) =  (n+1)*lpoly(ni,i) + x(i)*lsdata(ni,i) ;
      end
    end

   else
     error('Only shapes and first derivatives available at the moment');
   end
