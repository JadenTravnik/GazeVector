function [vL,vR] = Rawdata2W(el,vi,par) 
%EYELINK_RAWDATA2V map pupil positions into viewing direction in Vicon
%frame coordinates (centered in the eye)
% 
%   [vL,vR] = eyelink_rawdata2v(el,vi,par,type) 
%
%   type : 1-> gaze direction in Vicon coordinates; 2-> gaze direction in reference
%          eye frame coordinates
% 
%   eyelink_rawdata2v.m is currently locked by $Locker$ 
  
 %   by $Author$ Benedetta Cesqui
%   $Date$ 17-Apr-2012 14:59:21 

% gaze direction in head frame
[vhL,vhR] = Rawdata2Head(el.data,par);

% vicon frame
[T,X,Y,Z] = Vi2helmetframe(vi);

% gaze direction in vicon frame
ntime = length(el.time);
for i=1:ntime
    R = [X(:,i),Y(:,i),Z(:,i)];
    
    vL(:,i) = R*vhL(:,i);
    vR(:,i) = R*vhR(:,i);
    
end
end

function [vL,vR] = Rawdata2Head(eldata,par)
%EYELINK_RAWDATA2VH map pupil positions into viewing direction in head frame coordinates
%
%   [vL,vR] = eyelink_rawdata2vh(eldata,par)
%
%   eyelink_rawdata2vh.m is currently locked by $Locker$


%   by $Author$ Benedetta Cesqui
%   $Date$ 17-Apr-2012 14:59:21 

ntime = size(eldata,2);

k = par.camera(1);
f = par.camera(2);
Lx_off         = par.camera(3);
Ly_off         = par.camera(4);
Rx_off         = par.camera(5);
Ry_off         = par.camera(6);

for i=1:ntime
    
    % inverse projection of pupil center on camera plane
   
    % central projection
            
            % left eye
            xy = (eldata(1:2,i) - par.camera(3:4)')./([1;par.gainxy(1)].*par.gainLR(1));
            xy2 = sum(xy.^2);
            R = par.Leye2camera_R;
            T = -R'*par.Leye2camera_T;
            T1 = T(1);
            xyT = T(2:3)'*xy;
            T2r2 = sum(T.^2) - par.eye_radius^2;
            C(1) = 1 + xy2/k^2;
            C(2) = -2*T1 -2*f/k^2*xy2 + 2/k*xyT;
            C(3) = T2r2 + f^2/k^2*xy2 - 2*f/k*xyT;
            if C(2)^2 - 4*C(1)*C(3)>=0
                p1 = roots(C);
                p_1 = [p1(1);(f - p1(1))*xy/k];
                p_2 = [p1(2);(f - p1(2))*xy/k];
                n1 = sum(p_1.^2);
                n2 = sum(p_2.^2);
                if n1<n2
                    ve = (R*p_1+par.Leye2camera_T)/par.eye_radius;
                else
                    ve = (R*p_2+par.Leye2camera_T)/par.eye_radius;
                end
                vL(:,i) = par.head2eye_R*ve;
            else
                vL(:,i) = NaN*[1 1 1]';
            end
            
            % right eye
            
            xy = (eldata(3:4,i) - par.camera(5:6)')./([1;par.gainxy(2)].*par.gainLR(2));
            
            xy2 = sum(xy.^2);
            R = par.Reye2camera_R;
            T = -R'*par.Reye2camera_T;
            T1 = T(1);
            xyT = T(2:3)'*xy;
            T2r2 = sum(T.^2) - par.eye_radius^2;
            C(1) = 1 + xy2/k^2;
            C(2) = -2*T1 -2*f/k^2*xy2 + 2/k*xyT;
            C(3) = T2r2 + f^2/k^2*xy2 - 2*f/k*xyT;
            if C(2)^2 - 4*C(1)*C(3)>=0
                p1 = roots(C);
                p_1 = [p1(1);(f - p1(1))*xy/k];
                p_2 = [p1(2);(f - p1(2))*xy/k];
                n1 = sum(p_1.^2);
                n2 = sum(p_2.^2);
                if n1<n2
                    ve = (R*p_1+par.Reye2camera_T)/par.eye_radius;
                else
                    ve = (R*p_2+par.Reye2camera_T)/par.eye_radius;
                end
                vR(:,i) = par.head2eye_R*ve;
            else
                vR(:,i) = NaN*[1 1 1]';
            end
                
            
            
    end
end


function [T,X,Y,Z] = Vi2helmetframe(vi)
%Vi2helmetframe compute head reference frame based on three Vicon
%markers
%
%    [T,X,Y,Z] = Vi2helmetframe(vi,type,imrk,R0,T0)
%
%   T            : origin of helmet frame (marker 1)
%   X            : x axis of helmet frame (marker 2 - marker 1)
%   Y            : y axis of helmet frame (cross product of z axis and x axis)
%   Z            : z axis of helmet frame (cross product of x axis and (marker 3 - marker 1))))
%
%   vi           : structure with coordinates of 3 Vicon markers on the helmet
%   vi.m1(1:3,:) : 3D coordinates of marker 1
%   vi.m2(1:3,:) : 3D coordinates of marker 2
%   vi.m3(1:3,:) : 3D coordinates of marker 3
%   vi.time(1,:) : time stamp (Vicon clock) of marker data samples
%

% Assuming that v1.m1 .m2 .m3 exist
ntime = size(vi.m1,2);

% head reference origin (T) in Vicon coordinates -> marker 1
T = vi.m1;  % [3,ntime]

% head reference frame (X,Y,Z) in Vicon coordinates
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