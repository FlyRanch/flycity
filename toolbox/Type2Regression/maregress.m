function [b,bint,l,ang] = maregress(x,y,alpha)
%MAREGRESS Major Axis Regression (Principal Axis Regression).
%   Model II regression should be used when the two variables in the
%   regression equation are random and subject to error, i.e. not 
%   controlled by the researcher. Model I regression using ordinary least 
%   squares underestimates the slope of the linear relationship between the
%   variables when they both contain error. According to Sokal and Rohlf 
%   (1995), the subject of Model II regression is one on which research and
%   controversy are continuing and definitive recommendations are difficult
%   to make. In Sokal and Rohlf (1981, 2nd ed.), the numerical result for
%   major axis regression for the example data set is wrong. The mistake
%   has been corrected in the 1995 edition.
%
%   MAREGRESS is a Model II procedure. When both variables are in the same
%   units of measurement the slope of the major axis (principal axis) of
%   the bivariate sample has been sugested. This method requires a
%   knowledge of correlation techniques. A bivariate normal distribution
%   can be represented by means of concentric ellipses describing the
%   topography of a bell-shaped mound, and an ellipse can be described by
%   two principal axis (major and minor) at rigth angles toeach other.
%   Therefore the main task is to find the slope and equation of the major
%   axis from a sample.
%
%   [B,BINT,L,ANG] = MAREGRESS(X,Y,ALPHA) returns the vector B of 
%   regression coefficients in the linear Model II, a matrix BINT of the
%   given confidence intervals for B, the L eigenvalues of the axes, and 
%   the angle in degrees of the confidence bounds.
%
%   MAREGRESS treats NaNs in X or Y as missing values, and removes them.
%
%   Syntax: [b,bint,l,ang] = maregress(x,y,alpha)
%
%   Example. From the Box 14.12 (California fish cabezon [Scorpaenichthys 
%   marmoratus]) of Sokal and Rohlf (1995). The data are:
%
%   x=[14,17,24,25,27,33,34,37,40,41,42];
%   y=[61,37,65,69,54,93,87,89,100,90,97];
%
%   Calling on Matlab the function: 
%                [b,bint,l,ang] = maregress(x,y)
%
%   Answer is:
%
%   b =
%      6.6566    2.3017
%
%   bint =
%    -36.5408   27.8330
%      1.6043    3.7244
%
%   l =
%    494.6340   17.4933
%
%   ang =
%     66.5172
%
%   Created by A. Trujillo-Ortiz and R. Hernandez-Walls
%             Facultad de Ciencias Marinas
%             Universidad Autonoma de Baja California
%             Apdo. Postal 453
%             Ensenada, Baja California
%             Mexico.
%             atrujo@uabc.edu.mx
%
%   Copyright (C)  June 8, 2010. 
%
%   To cite this file, this would be an appropriate format:
%   Trujillo-Ortiz, A. and R. Hernandez-Walls. (2010). maregress: 
%      Major Axis Regression. A MATLAB file. [WWW document]. URL 
%      http://www.mathworks.com/matlabcentral/fileexchange/27916-maregress
%    
%   References:
%   Sokal, R. R. and Rohlf, F. J. (1995), Biometry. The principles and
%              practice of the statistics in biologicalreserach. 3rd. ed.
%              New-York:W.H.,Freeman. [Sections 14.13 and 15.7] 
%

if  nargin < 2
    error('maregress:TooFewInputs', ...
          'MAREGRESS requires at least two input arguments.');
elseif nargin == 2
    alpha = 0.05;
end

x = x(:); y = y(:);

% Check that matrix (X) and rigth hand side (Y) have compatible dimensions
[n,ncolx] = size(x);
if ~isvector(y)
    error('maregress:InvalidData', 'Y must be a vector.');
elseif numel(y) ~= n
    error('maregress:InvalidData', ...
          'The number of rows in Y must equal the number of rows in X.');
end

% Remove missing values, if any
wasnan = (isnan(y) | any(isnan(x),2));
havenans = any(wasnan);
if havenans
   y(wasnan) = [];
   x(wasnan,:) = [];
   n = length(y);
end

S = cov(x,y);
b1= ((S(2,2)-S(1,1))+sqrt(((S(2,2)-S(1,1))^2)+4*(S(1,2)^2)))/(2*S(1,2)); %slope
a = mean(y)-b1*mean(x); %intercept
b = [a,b1];
D = sqrt((S(1,1)+S(2,2))^2-4*((S(1,1)*S(2,2))-S(1,2)^2));
l1 = (S(1,1)+S(2,2)+D)/2; %eigenvalue of major (principal) axis
l2 = (S(1,1)+S(2,2)-D)/2; %eigenvalue of minor axis
H = (finv(1-alpha,1,n-2))/(((l1/l2)+(l2/l1)-2)*(n-2));
A = sqrt(H/(1-H));
Li = (b1-A)/(1+b1*A); %confidence lower limit of slope
Ls = (b1+A)/(1-b1*A); %confidence upper limit of slope
l = [l1,l2];
ai = mean(y)-Ls*mean(x); %confidence lower limit of intercept
as = mean(y)-Li*mean(x); %confidence upper limit of intercept
bint = [ai,as;Li,Ls];
ang = atand(b1); %angle in degrees

return,