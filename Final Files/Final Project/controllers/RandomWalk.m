%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D discrete random walk by Aniket N Prabhu %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear variables;
close all;

WalkStride = 100;  % Number of time steps
WalksNum = 100;      % Number of robots
x = zeros(WalkStride,1);
y = zeros(WalkStride,1);

for m=1:WalksNum
    for n=1:WalkStride
        STEP = sign(randn);  % generates +/-1 depending on the
                             % sign of the pseudo-random value drawn
                             % from a standard normal distribution
        if rand > 0.5  % pseudorandom value drawn from standard
                       % uniform distribution on the open interval
                       % (0,1). Reflects probability of
                       % walk along x (>0.5) or y (<0.5)
            x(n+1) = x(n) + STEP;
            y(n+1) = y(n);
        else
            x(n+1) = x(n);
            y(n+1) = y(n) + STEP;
        end
    end
    plot(x,y,'LineWidth',2)
    hold on
    grid on
end
axis square  % Keeps axis square even upon maximaizing