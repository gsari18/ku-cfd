clear; clc; close all;

nx=128;ny=32;
dt=0.001;nstep=4000;
mu=0.1;
maxit=5000;
beta=1.8;
h=1/ny;
u=zeros(nx+1,ny+2);v=zeros(nx+2,ny+1);p=zeros(nx+2,ny+2);
ut=zeros(nx+1,ny+2);vt=zeros(nx+2,ny+1);c=zeros(nx+2,ny+2)+0.25;
uu=zeros(nx+1,ny+1);vv=zeros(nx+1,ny+1);w=zeros(nx+1,ny+1);
force = zeros(1,nstep);
T_flux = zeros(1,nstep);
nu = zeros(1,nstep);

% Initial boundary conditions
u(1, :) = 1; % inlet velocity
T = zeros(nx+2,ny+2);
T(34:65, 1:17) = 1; % Temperature of the body
c(2,3:ny)=1/3;
c(nx+1,3:ny)=1/4;
c(3:nx,2)=1/3;
c(3:nx,ny+1)=1/3;
c(2,2)=1/2;c(2,ny+1)=1/2;
c(nx+1,2)=1/3;c(nx+1,ny+1)=1/3;

% c values for the body
c(33,2:17) = 1/3;
c(66,1:17) = 1/3;
c(34:65,18) = 1/3;
c(33,2) = 1/2;
c(66,2) = 1/2;
c(34:65,2:17) = 0;
time = 0;
for is = 1:nstep
    % Symmetry line (Full-slip boundary condition)
    u(1:nx+1,1) = u(1:nx+1,2);
    v(:,1) = 0;
    T(:,1) = T(:,2);

    % Top wall
    u(1:nx+1,ny+2) = -u(1:nx+1,ny+1);
    v(:,end) = 0;
    T(:,end) = T(:,end-1); % well insulated wall

    % Outlet
    u(end,:)=u(end-1,:);
    v(end,:)=v(end-1,:);
    T(end,:) = T(end-1,:);

    % Body (Top)
    u(34:64,17) = -u(34:64,18);
    v(34:65,17) = 0;

    % Body (Left)
    u(33,1:17) = 0;
    v(34,1:16) = -v(33,1:16);

    % Body (Right)
    u(65,1:17) = 0;
    v(64,1:16) = -v(65,1:16);

    ut = u;
    vt = v;

    for i=2:nx
        for j=2:ny+1 % temporary u-velocity
            ut(i,j) = u(i,j) + dt*(-(0.25/h)*((u(i+1,j) + u(i,j))^2 - (u(i,j)+...
                u(i-1,j))^2 + (u(i,j+1) + u(i,j))*(v(i+1,j) +...
                v(i,j)) - (u(i,j) + u(i,j-1))*(v(i+1,j-1) + v(i,j-1))) +...
                (mu/h^2)*(u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - 4*u(i,j)));
        end
    end

    for i=2:nx+1
        for j=2:ny % temporary v-velocity
            vt(i,j)=v(i,j) + dt*(-(0.25/h)*((u(i,j+1) + u(i,j))*(v(i+1,j) +...
                v(i,j)) - (u(i-1,j+1) + u(i-1,j))*(v(i,j) + v(i-1,j)) +...
                (v(i,j+1) + v(i,j))^2-(v(i,j) + v(i,j-1))^2) +...
                (mu/h^2)*(v(i+1,j) + v(i-1,j) + v(i,j+1) + v(i,j-1)-4*v(i,j)));
        end
    end

    for it=1:maxit % solve for pressure
        p_old = p;
        for i=2:nx+1
            for j=2:ny+1
                p(i,j)=beta*c(i,j)*(p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1)-...
                    (h/dt)*(ut(i,j)-ut(i-1,j)+vt(i,j)-vt(i,j-1)))+(1-beta)*p(i,j);
                
                p(34:65,1:17) = 0;
                p(end,:) = -p(end-1,:);
            end
        end
        
        if norm(p-p_old) < 0.0001
            break
        end
    end

    % correct the velocity
    u(2:nx,2:ny+1)=...
        ut(2:nx,2:ny+1)-(dt/h)*(p(3:nx+1,2:ny+1)-p(2:nx,2:ny+1));

    v(2:nx+1,2:ny)=...
        vt(2:nx+1,2:ny)-(dt/h)*(p(2:nx+1,3:ny+1)-p(2:nx+1,2:ny));

    time=time+dt;

    % Body (Top)
    u(34:64,17) = -u(34:64,18);
    v(34:65,17) = 0;

    % Body (Left)
    u(33,1:17) = 0;
    v(34,1:16) = -v(33,1:16);

    % Body (Right)
    u(65,1:17) = 0;
    v(64,1:16) = -v(65,1:16);

    % Body (in)
    u(34:64,1:17) = 0;
    v(35:64,1:16) = 0; 

    % Temperature
    T2 = T;
    for i=2:nx+1
        for j=2:ny+1
            T2(i,j) = T(i,j) + dt*( -(1/h)*(-(T(i-1,j)+T(i,j))*u(i-1,j)/2 + ...
                (T(i+1,j)+T(i,j))*u(i,j)/2 + (T(i,j+1)+T(i,j))*v(i,j)/2 + ...
                -(T(i,j-1)+T(i,j))*v(i,j-1)/2) + ...
                (0.1/h^2)*(T(i+1,j)+T(i-1,j)+T(i,j-1)+T(i,j+1)-4*T(i,j)));
        end
    end

    T = T2;
    T(34:65,1:17) = 1;

    % plot results
    uu(1:nx+1,1:ny+1)=0.5*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
    vv(1:nx+1,1:ny+1)=0.5*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
    w(1:nx+1,1:ny+1)=(u(1:nx+1,2:ny+2)-u(1:nx+1,1:ny+1)-...
        v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1))/(2*h);

    hold off;
    contourf(T', 'LineColor','none');
    colormap jet
    hold on;
    quiver(uu',vv')
    %contour(w',20);
    axis equal;
    title(num2str(is*dt))
    pause(0.001)

    % Force
    f_drag = sum(mu*(uu(33:65,18)-u(33:65,17)));
    f_pres = sum(p(33,2:17)*h - p(66,2:17)*h);
    force(is) = f_drag + f_pres;

    % Temperature flux
    T_flux(is) = (0.1/h^2)*(sum((T(66,1:17)-T(65,1:17))) + ...
        sum((T(33,1:17)-T(34,1:17))) +...
        sum((T(34:65,18)-T(34:65,17))));

    % Nusselt number
    nu(is) = 2*(sum(T(34,1:17)-T(33,1:17))/0.5+...
        sum(T(34:65,17)-T(34:65,18))/1+...
        sum(T(65,1:17)-T(66,1:17))/0.5);

end

figure

for i = 1:2:33
    strm = streamline(uu',vv',1,i);
    set(strm, 'Color', [14 173 190]/255)
    set(strm, 'LineWidth', 4)
    hold on
end

plot([1,nx+1],[ny+1,ny+1],'k-', 'LineWidth',2)
plot([1,nx+1],[1,1],'k--')
axis equal