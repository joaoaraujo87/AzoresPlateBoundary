function    X2v = X2v_okada(VAR,E,N,LENGTH,WIDTH,STRIKE,DIP,RAKE,Wd,data)
% chi square function for Okada dislocation
%
% VAR   vector with unknown parameters x0, y0, z0, op and phi
%
% SOURCE PARAMETERS
% OPEN        opening
% SLIP        slip
% ED          east deep motion
% ND          north deep motion
%
% CRUST PARAMETERS
% mu        shear modulus
% nu        Poisson's ratio 
%
% BENCHMARKS
% x,y       benchmark location
% z         depth within the crust (z=0 is the free surface)
%
% x0,y0     coordinates of the center of the sphere 
% z0        depth of the center of the sphere (positive downward and
%              defined as distance below the reference surface)
% nu    Poisson's ratio
% x,y   data point location
% Wd    weight matrix
% data  deformation data vector

% variables
OPEN = VAR(1); SLIP = VAR(2); DEPTH = VAR(3);

[U,V,W,~,~,~,~,~,~] = okada85_mod(E,N,DEPTH,STRIKE,DIP,LENGTH,WIDTH,RAKE,SLIP,OPEN);


model = [U V W]; model = model(:);
r = data - model;                       % residual
X2 = r'*Wd*r;                           % Chi Square - full covariance
p = length(VAR);                        % number of parameters
N = length(data);                       % number of data points
X2v = X2/(N-p);                         % Chi Square per degrees of freedom