function compass_zeroup(u, v, varargin)
%EARTH_COMPASS   Earth-coordinate compass plot
%                 - Changes plot labels to N,S,E,W
% 
%   Syntax:
%      EARTH_COMPASS(U, V, LINESPEC)
%
%   Inputs:
%      U:   vector of u-velocity
%      V:   vector of v-velocity
%   Optional:
%      LINESPEC: line specification (ie, color)
%
%

% Cameron Sparr
% cameronsparr@gmail.com
%
% This was inspired by wind_rose.m
%
% The code to change the axis labels from a mathematical system to
% N,S,E,W were written by Jordan Rosenthal in a forum post:
%      http://www.mathworks.com/matlabcentral/newsreader/author/4601
%
%
%

compass(u, v, varargin{:});

hHiddenText = findall(gca,'type','text');
Angles = 0 : 30 : 330;
hObjToDelete = zeros( length(Angles)-4, 1 );
k = 0;
for ang = Angles
   hObj = findall(hHiddenText,'string',num2str(ang));
   switch ang
   case 0
      set(hObj,'string',' 90','HorizontalAlignment','Left');
   case 90
      set(hObj,'string','0','VerticalAlignment','Bottom');
   case 180
      set(hObj,'string','-90','HorizontalAlignment','Right');
   case 270
      set(hObj,'string','180','VerticalAlignment','Top');
   otherwise
      k = k + 1;
      hObjToDelete(k) = hObj;
   end
end
delete( hObjToDelete(hObjToDelete~=0) );


end
