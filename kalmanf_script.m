% Kalman filter using simple example

%Build the input structure
clear s;
[zs, T] = size(np);
nstates = 4;
nob = 2;
s.A = eye(nstates);
s.A = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1]; 
s.z = np(:,1);
s.x = nan*zeros(nstates,1); 
s.H = eye(zs);
s.R = eye(zs); %measurement error covariance
%s.P = cov(np');
s.P = eye(zs);
s.u = 0;
s.B = eye(nstates);
s.Q = eye(nstates);

for t=1:T
    if t<3 
        s(end).u = zeros(zs,1);
    else 
        s(end).u = np(:,t-1) - np(:,t-2);
    end
    s(end).z = np(:,t);
    s(end+1)=kalmanf(s(end)); % perform a Kalman filter iteration
end

% figure
% hold on
% grid on
% % plot measurement data:
% hz=plot([s(1:end-1).z],'r.');
% % plot a-posteriori state estimates:
% hk=plot([s(2:end).x],'b-');
% ht=plot(tru,'g-');
% legend([hz hk ht],'observations','Kalman output','true voltage',0)
% title('Automobile Voltimeter Example')
% hold off

figure
hold on
grid on
% plot measurement data:
hz=plot(np(1,:), np(2,:),'rx');
% plot a-posteriori state estimates:
x = [s(2:end).x];
hk=plot(x(1,:), x(2,:),'b-o');

legend([hz hk],'observations','Kalman output',0)
title('Kalman filtering of positions')
hold off

