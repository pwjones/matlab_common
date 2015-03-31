function [bootMat, muStar]= stationaryBootstrap(N, B, data)
% Fucntion to sample for a stationary bootstrap
% B = 10; % Average block size
% N = 10; % Loop over N bootstraps
T = length(data(:));
yRepl = [data(:);data(:)];
u = zeros(T,1);
muStar = zeros(N,1);
bootMat = NaN*zeros(T,N);
for n=1:N
	u(1) = ceil(T*rand);
    for t=2:T
        if rand<1/B
        	u(t) = ceil(T*rand);
		else
			u(t) = u(t-1)+1;
		end 
	end
    % y-star sample simulation
    yStar = yRepl(u);
    bootMat(:,n) = yStar;
    % Mean of y-star
    muStar(n) = mean(yStar); 
end