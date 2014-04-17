function [ x,y,s ] = plotPolicy( policy, traj )
%PLOTPOLICY Summary of this function goes here
%   Detailed explanation goes here
mqlength = 5;
nu = 3;
nx = (mqlength+1)^3;

x = 1;
y = policy(traj(1));
s = traj(1);
for i=2:length(traj)
    if traj(i)~=traj(i-1)
        y = [y;policy(traj(i-1));policy(traj(i))];
        x = [x;i;i];
        s = [s;traj(i-1);traj(i)];
    else
        y = [y;policy(traj(i))];
        x = [x;i];
        s = [s;traj(i)];
    end    
end

figure
subplot(2,1,1)
plot(x(1:1000),y(1:1000))
subplot(2,1,2)
plot(x(1:1000),s(1:1000))
end

