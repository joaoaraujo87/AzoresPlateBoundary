function buriedGPS(problem)
close all; clc;

% this line must be commented before compiling
% clear all; path(pathdef);

% Jo?o D'Ara?jo Oct 2021
% This program is based on dMODELS software - Battaglia et al. (2013), J VOLCANOL GEOTH RES, 254,1-4.
% The program inverts GPS measurements to find the best fit buried rectangular dislocation (Okada 1985)
% with a 'deep motion'. The deep motion and azimuth are the magnitude and direction of the vector
% sum of the deep slip and opening motions.
%
% Change parameters in FIXED PARAMETERS as needed
%

addpath('C:\Program Files\MATLAB\R2017a\bin\Dmodels\B_Inverse_modeling\GPS_data_inversion-Buried\functions');

if nargin == 0
    problem = 'invert';
elseif strcmpi(problem,'help')
return;
end

addpath('C:\Program Files\MATLAB\R2017a\bin\Dmodels\B_Inverse_modeling\GPS_data_inversion-Okada\functions');

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

% *** READ DATA FILE  
fprintf(1,'------ LOAD DATA -------, loading %s\n',char(input_file));
fidinp = fopen(input_file,'r');
Mco = textscan(fidinp,'%q',1); comment = char(Mco{1});
Minp = textscan(fidinp,'%s %f %f %f %f %f %f %f %f','CommentStyle','%');
site = Minp{1}; xsite = Minp{2}; ysite = Minp{3}; 
east = Minp{4}; stdes = Minp{5}; north = Minp{6}; stdno = Minp{7}; 
up   = Minp{8}; stdup = Minp{9};

k = 0;
for i=1:length(site)
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
% =========================================================================

% allocate variables (this variables change size during the inversion)

bf_ef = zeros(maxk,1);
bf_chi = zeros(maxk,1);
bf_SLIP = zeros(maxk,1);  % slip (mm)
bf_OPEN = zeros(maxk,1);  % open (mm)

% ============================= FIXED PARAMETERS ==========================
% Values can be changed as needed

disp('------ INPUTS -------');

LENGTH = 100000;                                                                        % lenght of dislocation (km)
WIDTH = 100000;                                                                         % width of dislocation (km)
E = x/1000;                                                                             % east UTM coordinate of GPS stations (km)   
N = y/1000;                                                                             % north UTM coordinate of GPS stations (km)                                                                      
STRIKE = 112.5;                                                                         % strike (degrees)
DIP = 90;                                                                               % dip (degrees)
RAKE = 180;                                                                             % rake (degrees) is 180 (right-lateral slip) or -180 (left-lateral slip)
Smin = 0;                                                               % minimum search value for displacement (m)
Smax = 4.3/1000;                                                        % maximum search value for displacement (m)
Dmin = 1;                                                               % minimum search value for buried depth of dislocation (m)
Dmax = 20;                                                              % maximum search value for buried depth of dislocation (m)
	
% =========================== INVERSION ===================================

disp('------ DATA VECTOR -----');
data = [u v w]; data = data(:);                                             % deformation
sigma = [du dv dw]; sigma = sigma(:);                                       % uncertainty 

disp('------ COVARIACE MATRIX -----');
SigWN = diag(sigma.^2);                                                     % diagonal covariance matrix, m^2
SigRW = srw^2*yrw*eye(size(SigWN));                                         % random walk noise in m^2, (0.5 mm/sqrt(yr))
SigH  = SigWN + SigRW;                                                      % covariance matrix, m^2
WH = SigH\eye(size(SigH));                                                  % weight for least square, 1/m^2 (equivalent to inv(SigH))

h = waitbar(0,'Please wait ...');    
k = 1;
while k <= maxk
    OPENi = normrnd(Smin,Smax);                                             % opening start
    SLIPi = normrnd(Smin,Smax);                                             % slip start
    DEPTH_min = WIDTH/2 + Dmin;  
    DEPTH_max = WIDTH/2 + Dmax;
    DEPTHi = normrnd(DEPTH_min,DEPTH_max);
   [BF,fval,exitflag] = fmincon(@(VAR) X2v_buried(VAR,E,N,LENGTH,WIDTH,STRIKE,DIP,RAKE,WH,data),...
                        [OPENi,SLIPi,DEPTHi],[],[],[],[],...
                        [Smin,Smin,DEPTH_min],[Smax,Smax,DEPTH_max],[],options);

    if exitflag > 0
        bf_ef(k)  = exitflag;
        bf_chi(k) = fval;                                                   
        bf_OPEN(k)  = BF(1);                                                 
        bf_SLIP(k)  = BF(2);    
        bf_DEPTH(k)  = BF(3);    
        k = k+1;
    end
    
    waitbar(k/maxk);
end
close(h);

% determine best fit parameters
[chi2u,Iu] = min(bf_chi);
efu = bf_ef(Iu);  
OPEN = bf_OPEN(Iu); % meters
SLIP = bf_SLIP(Iu); % meters        
DEPTH = bf_DEPTH(Iu); % km   
% Displacements
[U,V,W,~,~,~,~,~,~] = okada85_dm(E,N,DEPTH,STRIKE,DIP,LENGTH,WIDTH,RAKE,SLIP,OPEN); % meters
% Residuals
R_U = u - U;
R_V = v - V;


disp('------ PLOT INVERSION STATISTICS ------');
figure('Position',[50 50 0.8*scsz(3) 0.8*scsz(4)],'Name','STATISTICS');                         
subplot(1,4,1); hist(bf_ef,10); title(strcat('exitflag: ',sprintf('%2.0f',efu)));
subplot(1,4,2); hist(bf_chi,100); title(strcat('\chi^2: ',sprintf(' %4.2f',chi2u)));
subplot(1,4,3); hist(bf_OPEN,100); title(strcat('OPEN: ',sprintf(' %4.2e',OPEN)));
subplot(1,4,4); hist(bf_SLIP,100); title(strcat('SLIP: ',sprintf(' %4.2e',SLIP)));

set(gcf,'PaperPositionMode','auto'); print -djpeg -r300 'GPSstatistics'; 


% =========================================================================

disp('------ PRINT OUTPUT FILE ------');
fidout = fopen(output_file,'w+');                                           % output data file fid
fprintf(fidout,'%s\n',char(comment));
fprintf(fidout,'Input data file: %s\n',input_file);
fprintf(fidout,'Source parameters\n');
fprintf(fidout,'X2v: %4.2e\n',chi2u);
fprintf(fidout,'OPEN SLIP\n');
fprintf(fidout,'%4.2f %4.2f\n',OPEN,SLIP);
 bm = char(bm);
 fprintf(fidout,'Site x y E mE dE N mN dN W mW dW\n');
 for i=1:length(x)
     fprintf(fidout,'%-5s %-8.0f %-8.0f %-6.4f %-6.4f %-6.4f %-6.4f %-6.4f %-6.4f %-6.4f %-6.4f %-6.4f\n',...
                       bm(i,:),x(i),y(i),u(i),U(i),du(i),v(i),V(i),dv(i),w(i),W(i),dw(i));
 end

disp('------ PLOT GPS VELOCITIES  ------')
hscale = 1e-7;                                                              % horizontal scale
site = char(bm);                                                            % site labels

% plot with GPS horizontal velocities
figure('Position',[10 10 0.8*scsz(4) 0.8*scsz(4)],'Name','GPS');
hold on; 
plot(0,0,'or','MarkerSize',10,'MarkerFaceColor','w');                       % plot source location
q1 = quiver(x,y,u/hscale,v/hscale,'r','LineWidth',2);                       % vector plot - data
q1.AutoScale = 'off';
hold on;
q2 = quiver(x,y,U/hscale,V/hscale,'b');                                     % vector plot - model
q2.AutoScale = 'off';                                                       
hold on;
q3 = quiver(x,y,R_U/hscale,R_V/hscale,'g');                                 % vector plot - residuals
q3.AutoScale = 'off';                                                       
title('Horizontal deformation [m] - red arrows: data');
set(gcf,'PaperPositionMode','auto')
print -djpeg -r300 -zbuffer GPShorivel; 

% *************************************************************************
fclose('all'); disp('That''s all Folks!');