close all;
clear; clc;

n = 201;
nstep = 1000;
length = 1.0;
h = length/(n-1);
dt = 0.001;
v = 0.01;
f = zeros(n, 1);
y = zeros(n, 1);

for i=1:n, f(i) = sin(2*pi*h*(i-1)) + 1.0; end % initial conditions

for m=1:nstep
    figure(1);
    plot(f, 'LineWidth', 2); ylim([-0.5 2.5]); xlim([0 201]);
    Title = sprintf('Time t = %5.3f s (Nonconservative)',dt*m); title(Title);
    y = f; % Store the solution
    for i=2:n-1
        f(i) = y(i) + dt*(v*(y(i+1)-2*y(i)+y(i-1))/h^2 - y(i)*(y(i+1)-y(i-1))/(2*h));
    end
    f(n) = y(n) + dt*(v*(y(2)-2*y(n)+y(n-1))/h^2 - y(n)*(y(2)-y(n-1))/(2*h));
    f(1) = f(n);
    pause(0.01);
end
