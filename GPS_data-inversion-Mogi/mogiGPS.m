function mogiGPS(~)
close all; clc;

% this line must be commented before compiling
% clear all; path(pathdef);

% add paths
addpath('C:\Program Files\MATLAB\R2017a\bin\Dmodels\B_Inverse_modeling\7_GPS_data-inversion-Mogi\functions')

% João D'Araújo Oct 2021
% This program is based on dMODELS software - Battaglia et al. (2013), J VOLCANOL GEOTH RES, 254,1-4.
% The program inverts GPS measurements to find the best fit point source (Mogi 1958).

% *** INITIALIZE PARAMETERS ===============================================
% fmisearch OPTIONS 
options = optimset('TolFun',1e-7,'TolX',1e-7,'MaxFunEvals',2000,'MaxIter',2000,'Algorithm','interior-point','Display','off');
% define figure size
scsz = get(0,'ScreenSize');
% =========================================================================


% =========================================================================
% *** READ INPUT PARAMETER FILE
fid0 = fopen('input_data_file.dat','r');
inputfile = textscan(fid0,'%s',1,'HeaderLines',2);
output_file = strcat(char(inputfile{1}),'.rsl');                               
input_file = strcat(char(inputfile{1}),'.dat');
% random walk
Mrw = textscan(fid0,'%f %f',1,'HeaderLines',5);
yrw = Mrw{1}; srw = Mrw{2};
% elastic constants  
Mec = textscan(fid0,'%f %f',1,'HeaderLines',5);
lambda = Mec{1}; mu = Mec{2}; nu=0.5*lambda/(lambda+mu);                    % lambda, shear modulus, Poisson's ratio
% number of random grid searches
Maxk = textscan(fid0,'%f',1,'HeaderLines',3); 
maxk = Maxk{1};
% sphere radius
Ma = textscan(fid0,'%f',1,'HeaderLines',3); 
a = Ma{1};
% search radius
Msr = textscan(fid0,'%f',1,'HeaderLines',3); 
SR = Msr{1};

% *** READ DATA FILE  
fprintf(1,'------ LOAD DATA -------, loading %s\n',char(input_file));
fidinp = fopen(input_file,'r');
Mco = textscan(fidinp,'%q',1); comment = char(Mco{1});
Minp = textscan(fidinp,'%s %f %f %f %f %f %f %f %f','CommentStyle','%');
site = Minp{1}; xsite = Minp{2}; ysite = Minp{3}; 
east = Minp{4}; stdes = Minp{5}; north = Minp{6}; stdno = Minp{7}; 
up   = Minp{8}; stdup = Minp{9};

% get location of vent (if available)
k = 0;
for i=1:length(site)
    if strcmp(site(i),'VENT') == 1
        Xvent = xsite(i); Yvent = ysite(i);
    else
        k=k+1;
        bm(k) = site(i);
        x(k)  = xsite(i);
        y(k)  = ysite(i);
        u(k)  = east(i);
        du(k) = stdes(i);
        v(k)  = north(i);
        dv(k) = stdno(i);
        w(k)  = up(i);
        dw(k) = stdup(i);
    end
end
% =========================================================================

% allocate variables (this variables change size during the inversion)
x0 = zeros(maxk,1); y0 = zeros(maxk,1); z0 = zeros(maxk,1); bf_chi = zeros(maxk,1);
bf_x0 = zeros(maxk,1); bf_y0 = zeros(maxk,1); bf_z0 = zeros(maxk,1);
bf_ef = zeros(maxk,1); bf_P_G = zeros(maxk,1); bf_dV = zeros(maxk,1);

% ==============================FIXED PARAMETERS===========================
% Values can be changed as needed

disp('------ INPUTS -------');

p1 = 'Type (1) for inflation or (0) for deflation: ';
i1 = input(p1);
% minimum depth
mind = 500;
% maximum depth
maxd = 5000;

% =========================== INVERSION ===================================
disp('------ DATA VECTOR -----');
data = [u v w]; data = data(:);                                             % deformation
sigma = [du dv dw]; sigma = sigma(:);                                       % uncertainty 

disp('------ COVARIACE MATRIX -----');
SigWN = diag(sigma.^2);                                                     % diagonal covariance matrix, m^2
SigRW = srw^2*yrw*eye(size(SigWN));                                         % random walk noise in m^2, (0.5 mm/sqrt(yr))
SigH  = SigWN + SigRW;                                                      % covariance matrix, m^2
WH = SigH\eye(size(SigH));                                                  % weight for least square, 1/m^2 (equivalent to inv(SigH))

% inflation case   
if i1 == 1                                                                                                                                                       
    dVmin = 0;                                                              % minumum volume m^3
    dVmax = 1E7;                                                            % maximum volume m^3
end
% deflation case                    
if i1 == 0  
    dVmin = -1E7;
    dVmax = 0;
end

disp('------ START MIN SEARCH -------');

h = waitbar(0,'Please wait ...');
k = 1;
while k <= maxk
   % initial values
   Xi = Xvent + normrnd(0,SR); Yi = Yvent + normrnd(0,SR);
   Zi = abs(normrnd(mind,maxd));
   dVi = normrnd(dVmin,dVmax);
   [BF,fval,exitflag] = fmincon(@(VAR) X2v_mogi(VAR,nu,x,y,WH,data),...
                        [Xi,Yi,Zi,dVi],[],[],[],[],...
                        [Xvent-SR,Yvent-SR,mind,dVmin],[Xvent+SR,Yvent+SR,maxd,dVmax],[],options);

    if exitflag > 0
        x0(k) = Xi; y0(k) = Yi; z0(k) = Zi; 
        bf_ef(k)  = exitflag;
        bf_chi(k) = fval;                                                   % best fit chi square
        bf_x0(k)  = BF(1);                                                  % best fit source location (East) 
        bf_y0(k)  = BF(2);                                                  % best fit source location (North)
        bf_z0(k)  = BF(3);                                                  % best fit depth, m 
        bf_dV(k) = BF(4);                                                   % best fit volume change
        k = k+1;
    end
    
    waitbar(k/maxk);
end
close(h);

% determine best fit parameters
[chi2u,Iu] = min(bf_chi);
X0 = bf_x0(Iu);                                                             % source location (east) 
Y0 = bf_y0(Iu);                                                             % source location (north)
Z0 = bf_z0(Iu);                                                             % source location (depth)
efu = bf_ef(Iu);                                                            % exit flag
dV = bf_dV(Iu);                                                             % volume change, m^3
UVW = mogi(X0,Y0,Z0,dV,nu,x,y);                                             % best fit solution
U = UVW(1,:);
V = UVW(2,:);
W = UVW(3,:);

disp('------ PLOT INVERSION STATISTICS ------');
figure('Position',[50 50 0.8*scsz(3) 0.8*scsz(4)],'Name','STATISTICS');                         
xi = linspace(min(x0),max(x0),10);                                             
yi = linspace(min(y0),max(y0),10);
zi = linspace(min(z0),max(z0),10);
subplot(2,5,1); hist(bf_ef,10); title(strcat('exitflag: ',sprintf('%2.0f',efu)));
subplot(2,5,2); hist(bf_chi,10); title(strcat('\chi^2: ',sprintf(' %4.2f',chi2u)));
subplot(2,5,3); hist(x0,xi); title('x_0');
subplot(2,5,4); hist(y0,yi); title('y_0');
subplot(2,5,5); hist(z0,zi); title('z_0');
subplot(2,5,6); hist(bf_dV,length(zi)); title(strcat('bf \DeltaV: ',sprintf(' %4.2e',dV))); 
subplot(2,5,7); hist(bf_x0,xi); title(strcat('bf x_0: ',sprintf('%8.0f',X0)));
subplot(2,5,8); hist(bf_y0,yi); title(strcat('bf y_0: ',sprintf('%8.0f',Y0)));
subplot(2,5,9); hist(bf_z0,zi); title(strcat('bf z_0: ',sprintf('%6.0f',Z0)));

set(gcf,'PaperPositionMode','auto'); print -djpeg -r300 'GPSstatistics'; 

% =========================================================================
disp('------ PRINT OUTPUT FILE ------');
fidout = fopen(output_file,'w+');                                           % output data file fid
fprintf(fidout,'%s\n',char(comment));
fprintf(fidout,'Input data file: %s\n',input_file);
fprintf(fidout,'Source parameters\n');
fprintf(fidout,'X2v: %4.2e\n',chi2u);
fprintf(fidout,'X0 Y0 Z0 radius dP(-) dV\n');
fprintf(fidout,'%8.0f %8.0f %6.0f %4.0f %4.3e\n',X0,Y0,Z0,a,dV);
bm = char(bm);
fprintf(fidout,'Site x y E mE dE N mN dN W mW dW\n');
 for i=1:length(x)
     fprintf(fidout,'%-5s %-8.0f %-8.0f %-6.4f %-6.4f %-6.4f %-6.4f %-6.4f %-6.4f %-6.4f %-6.4f %-6.4f\n',...
                       bm(i,:),x(i),y(i),u(i),U(i),du(i),v(i),V(i),dv(i),w(i),W(i),dw(i));
 end 

disp('------ PLOT BEST FIT  PROFILE------')
R = sqrt((x-X0).^2 + (y-Y0).^2);
UR = sqrt(U.^2+V.^2);                                                       % best fit radial deformation
ur = sqrt(u.^2+v.^2); dur = sqrt((u.*du./ur).^2+(v.*dv./ur).^2);            % radial deformation and error

figure('Position',[10 scsz(4)/10 scsz(3)/2 scsz(4)/2],'Name','PROFILE');
subplot(1,2,1)
errorbar(0.001*R,ur,dur,'r+'); hold on;                                     % data
plot(0.001*R,UR,'sb');                                                      % model
xlabel('R (km)'); title('Horizontal deformation');
legend('data','model'); legend('Location','best'); legend('boxoff');
subplot(1,2,2)
errorbar(0.001*R,w,dw,'r+'); hold on;                                       % data
plot(0.001*R,W,'sb');                                                       % model
xlabel('R (km)'); title('Vertical deformation');
legend('data','model'); legend('Location','best'); legend('boxoff');

set(gcf,'PaperPositionMode','auto')
print -djpeg -r300 -zbuffer GPSbestfit; 


disp('------ PLOT GPS VELOCITIES  ------')
hscale = max(sqrt(u.^2+v.^2)); vscale = max(w);                             % horizontal and vertical scale
site = char(bm);                                                            % site labels

% plot with GPS horizontal velocities
figure('Position',[10 10 0.8*scsz(4) 0.8*scsz(4)],'Name','GPS');
hold on; 
plot(X0,Y0,'or','MarkerFaceColor','y');                                     % plot source location
plot(x,y,'s','MarkerSize',4); text(x*1.0005,y,site,'FontSize',4);           % plot GPS sites
quiver(x,y,u/hscale,v/hscale,'r','LineWidth',2);                            % vector plot - data
quiver(x,y,U/hscale,V/hscale,'b');                                          % vector plot - model
title('Horizontal deformation [m] - red arrows: data');
set(gcf,'PaperPositionMode','auto')
print -djpeg -r300 -zbuffer GPShorivel; 

% plot with GPS vertical velocities
figure('Position',[10 10 0.8*scsz(4) 0.8*scsz(4)],'Name','VERT DEF');
hold on; 
plot(X0,Y0,'or','MarkerFaceColor','y');                                     % plot source location
plot(Xvent,Yvent,'or','MarkerSize',10,'MarkerFaceColor','w');               % plot crater location
plot(x,y,'s','MarkerSize',4); text(x*1.0005,y,site,'FontSize',4);           % plot GPS sites
quiver(x,y,zeros(size(W)),w/vscale,'r','LineWidth',2);                      % vector plot - data
quiver(x,y,zeros(size(W)),W/vscale,'b');                                    % vector plot - model
title('Vertical deformation [m] - red arrows: data');
set(gcf,'PaperPositionMode','auto')
print -djpeg -r300 -zbuffer GPSvertvel; 

% *************************************************************************
fclose('all'); disp('That''s all Folks!');