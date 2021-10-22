function    X2v = X2v_mogi(VAR,nu,x,y,Wd,data)
% chi square function for a mogi source

% VAR   vector with unknown parameters x0, y0, z0 and dV
% x0,y0 coordinates of the center of the sphere 
% z0    depth of the center of the sphere (positive downward and
%              defined as distance below the reference surface)
% dV    volume change
% b     radius of the sphere
% nu    Poisson's ratio
% x,y   data point location
% Wd    weight matrix
% data  deformation data vector

x0 = VAR(1); y0 = VAR(2); z0 = VAR(3); dV = VAR(4);
U = mogi(x0,y0,z0,dV,nu,x,y);
u = U(1,:);
v = U(2,:);
w = U(3,:);
model = [u v w]; model = model(:);
r = data - model;                       % residual
X2 = r'*Wd*r;                           % Chi Square - full covariance
p = length(VAR);                        % number of parameters
N = length(data);                       % number of data points
X2v = X2/(N-p);                         % Chi Square per degrees of freedom