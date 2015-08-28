function [par] = GazeinSpaceCalibrateFun (eldata, vidata,opt)
% GazeinSpaceCalibrateFun estimates the geometrical parameters for the mapping between pupil
% image on the camera CCD of an Eye tracker and the relative gaze
% orientation. Additionally, the method includes a procedure to correct for possible slippages of the tracker helmet
% based on a geometrical representation  of the pupil-to-gaze mapping.
%
%
%  [par] = GazeinSpaceCalibrateFun ( eldata, vidata, tdata, opt)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INSTRUCTION FOR USERS
% Inputs should be organize as follows
%   eldata -   structure with Left and and Right eye (x,y) pupil coordinates measured with the Eye tracker.
%
%               eldata.static  - raws of signal from the Eye tracker static trial.
%                                eldata.static.data(1,:) left eye x coordinate(RAW)
%                                eldata.static.data(2,:) left eye y coordinate(RAW)
%                                eldata.static.data(3,:) right eye x coordinate(RAW)
%                                eldata.static.data(4,:) right eye y coordinate(RAW)
%                                eldata.static.time(1,:) time stamp eye data sample
%
%               eldata.cal    -  raws of signal from the Eye tracker calibration trial.
%                                eldata.cal.data(1,:) left eye x coordinate(RAW)
%                                eldata.cal.data(2,:) left eye y coordinate(RAW)
%                                eldata.cal.data(3,:) right eye x coordinate(RAW)
%                                eldata.cal.data(4,:) right eye y coordinate(RAW)
%                                eldata.cal.time(1,:) time stamp eye data sample
%
%   vidata -   structure with coordinates of the 3 Vicon markers placed on  the helmet plus the marker on the target
%
%              vidata.static   - data relative to the static trial
%                                 vidata.static.m1(1:3,:) : 3D coordinates of marker 1
%                                 vidata.static.m2(1:3,:) : 3D coordinates of marker 2
%                                 vidata.static.m3(1:3,:) : 3D coordinates of marker 3
%                                 vidata.static.g(1:3,:)  : 3D coordinates of fixation marker
%                                 vidata.static.time        time stamp eye data sample
%              vidata.cal      - data relative to the calibration trial
%                                 vidata.cal.m1(1:3,:) : 3D coordinates of marker 1
%                                 vidata.cal.m2(1:3,:) : 3D coordinates of marker 2
%                                 vidata.cal.m3(1:3,:) : 3D coordinates of marker 3
%                                 vidata.cal.g(1:3,:)  : 3D coordinates of fixation marker
%                                 vidata.cal.time        time stamp eye data sample
%
%   opt : Initialization parameters and algorthm settings. 
%         opt.par_type  -> 'standard';'drift;
%         opt.data_type -> 'CR'; 'PUPIL'
%         opt.par       -> calibration parameter
%         opt.range     -> Upper and lower boundary of variation
%         if opt.par.type    = 'standard' -> set initial parameters in opt =  initialize_options()
%         if opt.par.type    = 'drift'    -> set opt.par = initial par. Use the data in eldata.cal and vidata.cal to correct
%   
%
% Outputs :
%           par.type           : 1-> standard; 2-> drift
%           par.eye_radius     : eye radius
%           par.eye_iod        : eye interocular distance
%           par.helmet2eye_R   : helmet to eye frame rotation matrix
%           par.helmet2eye_T   : helmet to cyclopean eye translation vector
%           par.helmet2Leye_T  : helmet to L eye frame translation vector
%           par.helmet2Reye_T  : helmet to R eye frame translation vector
%           par.Leye2camera_Rf : L eye to camera rotation matrix Fick angles
%           par.Leye2camera_R  : L eye to camera rotation matrix
%           par.Leye2camera_T  : L eye to camera translation vector
%           par.helmet2Reye_T  : helmet to R eye frame translation vector
%           par.Reye2camera_Rf : R eye to camera rotation matrix Fick angles
%           par.Reye2camera_R  : R eye to camera rotation matrix
%           par.Reye2camera_T  : R eye to camera translation vector
%           par.camera(1:6)    : [k, f, Lx_off, Ly_off, Rx_off, Ry_off]
%           par.gainxy         : L and R eye aspect ration parameter
%           par.gainLR         : L and R eye scaling factor parameter 
%           par.Rf_drift       : Drift correction parameter ;
%           par.helmet2skull_T : Drift correction parameter ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors 
% B. Cesqui, R. van de Langenberg, F. Lacquaniti, and A. d'Avella.
% Modifications
% 26/03/2013 GazeinSpaceCalibrateFun is first created

if nargin<3, opt=[]; end
defopt = initialize_options;
opt=opt_setdef(opt,defopt);
switch opt.par_type
    case 'standard'
        % Initialize Calibration parameters
        par0 = Ref2par(eldata.static,vidata.static);
        inpar = opt.par;
        
        switch opt.data_type
            case 'CR'
                inpar.camera(1) = 50000;
            case 'PUPIL'
                inpar.camera(1) = 75000;
        end
        
        inpar.camera(3:6) = par0.offset;
        %Define Upper and lower boundaries:
        
        inpar.helmet2eye_Rf         = par0.helmet2eye_Rf;
        inpar.helmet2eye_R          = Fick2R(par0.helmet2eye_Rf);
        par_range.camera            = [1;1]* inpar.camera(1:2) + [-1;1]*opt.range.camera;
        par_range.gainxy            = [1;1]*inpar.gainxy + [-1;1]*opt.range.gainxy;
        par_range.gainLR            = [1;1]*inpar.gainLR + [-1;1]*opt.range.gainLR;
        par_range.eye_radius        = inpar.eye_radius + [-1 1]*opt.range.eye_radius;
        par_range.eye_iod           = inpar.eye_iod + [-1;1]*opt.range.eye_radius;
        par_range.Leye2camera_Rf    = ([1;1]*inpar.Leye2camera_Rf + [-1;1]* opt.range.Leye2camera_Rf);
        par_range.Leye2camera_T     = ([1;1]*inpar.Leye2camera_T + [-1;1]*opt.range.Leye2camera_T);
        par_range.Reye2camera_Rf    = ([1;1]*inpar.Reye2camera_Rf + [-1;1]* opt.range.Reye2camera_Rf);
        par_range.Reye2camera_T     = ([1;1]*inpar.Reye2camera_T + [-1;1]*opt.range.Reye2camera_T);
        par_range.helmet2eye_T      = ([1;1]*inpar.helmet2eye_T+ [-1;1]*opt.range.helmet2eye_T);
        

    case 'drift'
        
        inpar = opt.par;
        inpar.helmet2skull_T        = inpar.helmet2skull_T;
        par_range.Rf_drift          = ([1;1]*inpar.Rf_drift + [-1;1]* opt.range.Rf_drift);
        par_range.helmet2skull_T    = ([1;1]*inpar.helmet2skull_T + [-1;1]*opt.range.helmet2skull_T);


end

inpar.type = opt.par_type;
par =  ParExtraction(eldata.cal,vidata.cal,inpar,par_range);
par.par_range=par_range;

end

function par = ParExtraction(el,vi,par,par_range)


switch par.type
    case 'standard'
        optpar(1)       = par.eye_radius;
        optpar(2:4)     = par.helmet2eye_T;
        optpar(5:7)     = par.Leye2camera_Rf; 
        optpar(8:10)    = par.Leye2camera_T;
        optpar(11:13)   = par.Reye2camera_Rf; 
        optpar(14:16)   = par.Reye2camera_T;
        optpar(17:18)   = par.gainxy;  
        optpar(19)      = par.eye_iod;
        optpar(20:21)   = par.gainLR;

        %Upper and lower boundaries
        lb(1)           = par_range.eye_radius(1);
        lb(2:4)         = par_range.helmet2eye_T(1,:);
        lb(5:7)         = par_range.Leye2camera_Rf(1,:); 
        lb(8:10)        = par_range.Leye2camera_T(1,:);
        lb(11:13)       = par_range.Reye2camera_Rf(1,:); 
        lb(14:16)       = par_range.Reye2camera_T(1,:);
        lb(17:18)       = par_range.gainxy(1,:);
        lb(19)          = par_range.eye_iod(1);
        lb(20:21)       = par_range.gainLR(1,:);
        
        ub(1)           = par_range.eye_radius(2);
        ub(2:4)         = par_range.helmet2eye_T(2,:);
        ub(5:7)         = par_range.Leye2camera_Rf(2,:); 
        ub(8:10)        = par_range.Leye2camera_T(2,:);
        ub(11:13)       = par_range.Reye2camera_Rf(2,:); 
        ub(14:16)       = par_range.Reye2camera_T(2,:);
        ub(17:18)       = par_range.gainxy(2,:); 
        ub(19)          = par_range.eye_iod(2);
        ub(20:21)       = par_range.gainLR(2,:);
        
    case 'drift'
          
        optpar(1:3) = par.Rf_drift; 
        optpar(4:6) = par.helmet2skull_T; 
         %Upper and lower boundaries
        lb(1:3) = par_range.Rf_drift(1,:);
        lb(4:6) = par_range.helmet2skull_T(1,:);
        ub(1:3) = par_range.Rf_drift(2,:);
        ub(4:6) = par_range.helmet2skull_T (2,:);
        
        
end

options = optimset('Display','off','MaxFunEvals',1e5);
optpar = fmincon(@calibrate_fun,optpar,[],[],[],[],lb,ub,[],options,el,vi,par);

switch par.type
    case 'standard'
        
        par.eye_radius      = optpar(1);
        helmet2eye_R        = eyelink_fick2R(par.helmet2eye_Rf);
        par.helmet2eye_T    = optpar(2:4)';
        par.helmet2Leye_T   = par.helmet2eye_T +  helmet2eye_R* (optpar(19)/2*[0 1 0]');
        par.helmet2Reye_T   = par.helmet2eye_T -  helmet2eye_R* (optpar(19)/2*[0 1 0]');
        par.Leye2camera_Rf  = cast_angle(optpar(5:7)); 
        par.Leye2camera_R   = Fick2R(optpar(5:7));
        par.Leye2camera_T   = optpar(8:10)';
        par.Reye2camera_Rf  = cast_angle(optpar(11:13)); 
        par.Reye2camera_R   = Fick2R(optpar(11:13));
        par.Reye2camera_T   = optpar(14:16)';
        par.gainxy          = optpar(17:18); 
        par.eye_iod         = optpar(19);
        par.gainLR          = optpar(20:21);
      
        
    case 'drift'
        
        R_drift               = eyelink_fick2R(optpar(1:3));
        helmet2skull_T        = optpar(4:6)';
        helmet2eye_R          = R_drift*par.helmet2eye_R;
        helmet2Leye_T         = R_drift*(par.helmet2Leye_T - helmet2skull_T) + helmet2skull_T;
        helmet2Reye_T         = R_drift*(par.helmet2Reye_T - helmet2skull_T) + helmet2skull_T;
        helmet2Rcamera_T      =(par.helmet2eye_R * par.Reye2camera_T)  + par.helmet2Reye_T;
        helmet2Lcamera_T      =(par.helmet2eye_R * par.Leye2camera_T)  + par.helmet2Leye_T;
        
        
        par.Leye2camera_R     = par.Leye2camera_R * R_drift';
        par.Reye2camera_R     = par.Reye2camera_R * R_drift';
        par.Reye2camera_T     = helmet2eye_R'*(helmet2Rcamera_T - helmet2Reye_T);
        par.Leye2camera_T     = helmet2eye_R'*(helmet2Lcamera_T - helmet2Leye_T);
        par.helmet2eye_R      = helmet2eye_R ;
        par.helmet2eye_Rf     = R2fick(par.helmet2eye_R);
        par.helmet2Leye_T     = helmet2Leye_T;
        par.helmet2Reye_T     = helmet2Reye_T;
        par.Rf_drift          = R2fick(R_drift);
        par.helmet2skull_T    = helmet2skull_T;
        
end

par.rms            = calibrate_fun(optpar,el,vi,par);
[rms]              = ReconstructionErrorPlot(par,el,vi,1);

end

function rms = calibrate_fun(optpar,el,vi,par)
%calibrate2_fun comupute error in pupil projected position (Eye tracker RAW data) 

switch par.type
    
    case 'standard'

        par.eye_radius       = optpar(1);
        par.helmet2Leye_T    = optpar(2:4)'+  par.helmet2eye_R* (optpar(19)/2*[0 1 0]');
        par.Leye2camera_Rf   = optpar(5:7);
        par.Leye2camera_R    = Fick2R(optpar(5:7));
        par.Leye2camera_T    = optpar(8:10)';
        par.helmet2Reye_T    = optpar(2:4)'- par.helmet2eye_R*(optpar(19)/2*[0 1 0]');
        par.Reye2camera_Rf   = optpar(11:13);
        par.Reye2camera_R    = Fick2R(optpar(11:13));
        par.Reye2camera_T    = optpar(14:16)';
        par.gainxy           = optpar(17:18);
        par.gainLR           = optpar(20:21);
        
    case 'drift'
              
        R_drift               = eyelink_fick2R(optpar(1:3));
        helmet2skull_T        = optpar(4:6)';
        helmet2eye_R          = R_drift*par.helmet2eye_R;
        helmet2Leye_T         = R_drift*(par.helmet2Leye_T - helmet2skull_T) + helmet2skull_T;
        helmet2Reye_T         = R_drift*(par.helmet2Reye_T - helmet2skull_T) + helmet2skull_T;
        helmet2Rcamera_T      =(par.helmet2eye_R * par.Reye2camera_T)  + par.helmet2Reye_T;
        helmet2Lcamera_T      =(par.helmet2eye_R * par.Leye2camera_T)  + par.helmet2Leye_T;
        
        
        par.Leye2camera_R     = par.Leye2camera_R * R_drift';
        par.Reye2camera_R     = par.Reye2camera_R * R_drift';
        par.Reye2camera_T     = helmet2eye_R'*(helmet2Rcamera_T - helmet2Reye_T);
        par.Leye2camera_T     = helmet2eye_R'*( helmet2Lcamera_T - helmet2Leye_T);
        par.helmet2eye_R      = helmet2eye_R ;
        par.helmet2eye_Rf     = R2fick(par.helmet2eye_R);
        par.helmet2Leye_T     = helmet2Leye_T;
        par.helmet2Reye_T     = helmet2Reye_T;
        par.Rf_drift          = R2fick(R_drift);
        par.helmet2skull_T    = helmet2skull_T;
               
end

[rms] = ReconstructionErrorPlot(par,el,vi,0);

end

function [rms,resid,datahat] = ReconstructionErrorPlot(par,el,vi,figureplot)

if nargin<4 figureplot=1;end
% Assume that the same number of eldata and gazedata samples are passed.
eldata = el.data;
ntime = size(eldata,2);
gazedata =  Vi2gh(vi);

resid = zeros(4,ntime);
CM = [0 1 0; 0 0 1]; % camera projection matrix
datahat = nan(4,length(eldata(1,:))) ;

k              = par.camera(1);
f              = par.camera(2);
Lx_off         = par.camera(3);
Ly_off         = par.camera(4);
Rx_off         = par.camera(5);
Ry_off         = par.camera(6);

for i=1:ntime
    gh = gazedata(:,i);
    % compute gaze target in eye frame
    ge = par.helmet2eye_R'*gh;  % rotation only
    gLe = ge - par.helmet2eye_R'*par.helmet2Leye_T;
    gRe = ge - par.helmet2eye_R'*par.helmet2Reye_T;
    
    % compute pupil center in eye frame
    pLe = gLe/norm(gLe)*par.eye_radius;
    pRe = gRe/norm(gRe)*par.eye_radius;
    
    % compute pupil center in camera frame
    pLc = par.Leye2camera_R'*(pLe-par.Leye2camera_T);
    pRc = par.Reye2camera_R'*(pRe-par.Reye2camera_T);
    
    % Use the central projection
    kL = k / (f - pLc(1));
    kR = k / (f - pRc(1));
    pLp = kL*CM*pLc;
    pRp = kR*CM*pRc;
    
    
    eldatahat(1:2,1) = par.gainLR(1)*(pLp.*[1;par.gainxy(1)])+ [Lx_off;Ly_off]; % Left  eye
    eldatahat(3:4,1) = par.gainLR(2)*(pRp.*[1;par.gainxy(2)])+ [Rx_off;Ry_off]; % Right eye
    datahat(1:2,i)   = par.gainLR(1)*(pLp.*[1;par.gainxy(1)])  + [Lx_off;Ly_off];
    datahat(3:4,i)   = par.gainLR(2)*(pRp.*[1;par.gainxy(2)])  + [Rx_off;Ry_off];
    resid(:,i)       = eldata(:,i) - eldatahat;
    
end
rms = sqrt(sum(resid(:).^2)/ntime/4);

if figureplot
    figure
    subplot(2,2,1)
    plot(el.time, eldata(1,:),'k');
    hold on
    plot(el.time,datahat(1,:),'color',[0.5 0.5 0.5]);
    legend('registred', 'esteemed');
    title('x left');
    
    
    subplot(2,2,2)
    plot(el.time, eldata(2,:));
    hold on
    plot(el.time,datahat(2,:),'color',[0.5 0.5 0.5]);
    legend('registred', 'esteemed');
    title('y left');
    
    subplot(2,2,3)
    plot(el.time, eldata(3,:));
    hold on
    plot(el.time,datahat(3,:),'color',[0.5 0.5 0.5]);
    legend('registred', 'esteemed');
    title('x right');
    
    subplot(2,2,4)
    plot(el.time, eldata(4,:));
    hold on
    plot(el.time,datahat(4,:),'color',[0.5 0.5 0.5]);
    legend('registred', 'esteemed');
    title('y right')
end

end

function R = Fick2R(fickangles) 
%EYELINK_FICK2R compute rotation matrix R from Fick angles (XYZ)
% 
%   R = eyelink_fick2R(fickangles)
%
%    fickangles: theta (Z), phi (Y) ,psi (X)
%
theta = fickangles(1);
phi = fickangles(2);
psi = fickangles(3);

% Z
c = cos(theta);
s = sin(theta);
R1 = [c -s 0; s c 0; 0 0 1];

% Y
c = cos(phi);
s = sin(phi);
R2 = [c 0 s; 0 1 0; -s 0 c];

% X
c = cos(psi);
s = sin(psi);
R3 = [1 0 0; 0 c -s; 0 s c];

R = R1 * R2 * R3;  % Fick convention
end

function fickangles = R2fick(R) 
%EYELINK_R2FICK compute Fick angles (XYZ) from rotation matrix R


if abs(R(3,1))==1 
  
  theta = 0;
  delta = atan2(R(1,2),R(1,3));
  if R(3,1)==-1
    phi = pi/2;
    psi = theta + delta;
  else
    phi = -pi/2;
    psi = -theta+delta;
  end
  
else
  
  phi = -asin(R(3,1));
  %  phi2 = pi - phi;
  psi = atan2(R(3,2)/cos(phi),R(3,3)/cos(phi));
  % psi2 = atan2(R(3,2)/cos(phi2),R(3,3)/cos(phi2));
  theta = atan2(R(2,1)/cos(phi),R(1,1)/cos(phi));
  % theta2 = atan2(R(2,1)/cos(phi2),R(1,1)/cos(phi2));
  
end

fickangles = [theta, phi, psi];
end

function gh = Vi2gh(vi) 
%Vi2gh compute gaze target (vi.g) in helmet coordinate frame (using helmet markers vi.m1, vi.m2, vi.m3)
% compute helmet frame


[T,X,Y,Z] = Vi2helmetframe(vi);

ntime = length(vi.time);
gh = zeros(3,ntime);

for i=1:ntime
  
  % compute gaze target in helmet coordinate frame
  Ri = [X(:,i),Y(:,i),Z(:,i)];
  gh(:,i) = Ri'*(vi.g(:,i)-T(:,i));
  
end
end

function par = Ref2par(el,vi,imrk) 
% Ref2par compute RAW offsets and helmet2eye_R from static posture 
 if nargin<3 || isempty(imrk), imrk = [2 3]; end % Assume m2 and m3 the two markers placed at the extremities of the helmetband.

% Assuming that v1.m1 .m2 .m3 exist
ntime = size(vi.m1,2);
% Computes the average value over the time interval
viav.time = 0;
for j=1:3
  M = zeros(3,ntime);
  for i=1:ntime
    M(:,i) = vi.(sprintf('m%i',j))(:,i);
  end
  T(:,j) = mean(M,2);
  viav.(sprintf('m%i',j)) = T(:,j);
end


E = (T(:,imrk(1))+T(:,imrk(2)))/2;  % Compute mid point between the two markers on the helmetband of the helmet    
G = mean(vi.g,2);
viav.g = G;
elav.data = mean(el.data,2);
% Compute eye frame assuming Z axis is vertical and marker 2

axisx = G-E; %x-axis is the line from E to gazed target, placed at eye height.
axisx(3) = 0; %force the x-axis to be horizontal
axisx = axisx/norm(axisx);
axisz = [0 0 1]';
axisy = cross(axisz,axisx);
% rotation matrix consisting of eye frame in vicon coordinates
Rve = [axisx axisy axisz];

% Get helmet frame in ref posture and gaze direction
[T,X,Y,Z] = Vi2helmetframe(viav); 

Rvh = [X,Y,Z];
Rhv = Rvh';
par.helmet2eye_R = Rhv*Rve;  % helmet to eye rotation
par.helmet2eye_Rf = eyelink_R2fick(par.helmet2eye_R);
par.vicon2helmet_R0 = Rvh;
% offset is mean raw pupil data at ref, also use default k and f
par.offset = [elav.data'];  % [ Lx offset, L y offset, R x offset, R y offset] 
end

function [T,X,Y,Z] = Vi2helmetframe(vi)
%Vi2helmetframe compute helmet reference frame based on three Vicon markers


% Assuming that v1.m1 .m2 .m3 exist
ntime = size(vi.m1,2);

% helmet reference origin (T) in Vicon coordinates -> marker 1
T = vi.m1;  % [3,ntime]

% helmet reference frame (X,Y,Z) in Vicon coordinates
% x axis -> marker 2 - marker 1
X = vi.m2 - vi.m1;  % [3,ntime]
Xnorm = sqrt(sum(X.^2));
X = X ./ (ones(3,1)*Xnorm);

% z axis -> cross product of x axis and (marker 3 - marker 1)
Z = cross(X,vi.m3-vi.m1,1);
Znorm = sqrt(sum(Z.^2));
Z = Z ./ (ones(3,1)*Znorm);

% y axis -> cross product of z axis and x axis
Y = cross(Z,X);

for i=1:ntime
    Ri = [X(:,i),Y(:,i),Z(:,i)];
    Ri = Ri;
    X(:,i) = Ri(:,1);
    Y(:,i) = Ri(:,2);
    Z(:,i) = Ri(:,3);
end
end

function x = cast_angle(x,th0)
% force angle x into th0 + [0 2*pi]
if nargin<2, th0 = -pi; end
x = mod(x - th0,2*pi) + th0;
end

function opt = opt_setdef(opt,defopt)
%OPT_SETDEF set opt field corresponding to defopt fields. 
%   If present opt fields are used.
%
%   opt = opt_setdef(opt,defopt)
%
%   opt_setdef.m is currently locked by $Locker$

if nargin<2, return, end

optfield = fieldnames(defopt);
for i=1:length(optfield)
    if isfield(opt,optfield{i})
        optfield_i = getfield(opt,optfield{i});
        if isstruct(optfield_i) % recursion in subfields
            defoptfield_i = getfield(defopt,optfield{i});
            optfield_i = opt_setdef(optfield_i,defoptfield_i);
        end
        defopt = setfield(defopt,optfield{i},optfield_i);
    end
end
opt = defopt;
end

function opt = initialize_options()

% Set parameters initial values and contraints. Insert measured and
% estimated values. Angles expressed in radians, distances expressed in
% meters.

opt.par_type ='standard'; % other options -> 'drift';
opt.data_type = 'CR'; % other options -> 'PUPIL';
opt.par.camera(1:2)              = [50000 .0075]; %[c.u, m]
opt.par.gainxy                  = [.85 .85];
opt.par.gainLR                  = [1 1];
opt.par.proj_type = 1;
opt.par.Kp =  650000;
opt.par.eye_radius              = 0.012; %[m]
opt.par.eye_iod                 = 0.065; %[m]
opt.par.helmet2eye_R            = [];
opt.par.helmet2eye_Rf           = [];%[rad]
opt.par.Leye2camera_Rf          = [0 30 180]*pi/180;%[rad]
opt.par.Reye2camera_Rf          = [0 30 180]*pi/180;%[rad]
opt.par.Leye2camera_R           = Fick2R(opt.par.Leye2camera_Rf);
opt.par.Reye2camera_R           = Fick2R(opt.par.Reye2camera_Rf);
opt.par.Leye2camera_T           = [.06 0 -.035]; %[m]
opt.par.Reye2camera_T           = [.06 0 -.035]; %[m]
opt.par.helmet2eye_T            = [.08 .11 -.04];%[m]
opt.par.helmet2Leye_T           = [];
opt.par.helmet2Reye_T           = [];
opt.range.camera                = [5000 .001 ]; % [c.u, m];
opt.range.gainxy                = [.5 .5];  % u;
opt.range.gainLR                = [.5 .5];  % u;
opt.range.eye_radius            = 0.003; % [m];
opt.range.eye_iod               = 0.005; % [m];
opt.range.Leye2camera_Rf        = [20 20 20]*pi/180;  % rad
opt.range.Leye2camera_T         = [.02 .02 .02];  % [m]
opt.range.Reye2camera_Rf        = [20 20 20]*pi/180;  % rad
opt.range.Reye2camera_T         = [.02 .02 .02];  % [m]
opt.range.helmet2eye_T          = [.02 .02 .02]; % [m]
opt.range.helmet2eye_Rf         = [20 20 20]*pi/180; % rad
opt.par.Rf_drift                = [0 0 0];
opt.par.helmet2skull_T          = [0 .12 -.05];  % [m] rotation center approx in the center of skull par range
opt.range.Rf_drift              = [10 10 10]*pi/180; % [m]
opt.range.helmet2skull_T        = [.005 .005 .005]; %[m]
end