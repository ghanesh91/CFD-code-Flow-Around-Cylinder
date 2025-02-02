clear
%Grid definitions
Nx=129;Ny=65;
Lx=10;Ly=5;
kfunc=1e12;
rad=.5;
xo=2;yo=2.5;
%parpool('local',2)
%%%%%%%%%%%nodal points%%%%%%%%%%%%%%%%%

x=linspace(0,Lx,Nx);
y=linspace(0,Ly,Ny);
dx=x(2)-x(1);dy=y(2)-y(1);
  for j=1:Ny
    for i=1:Nx
    X(j,i)=x(i);
    Y(j,i)=y(Ny+1-j);
    end
  end
  
dt=1e-2;
t=0:dt:.1;
ntime=1;
Re=150;
if ntime==1
%define velocities at nodal points
u(1:Ny,1:Nx,ntime)=0;
v(1:Ny,1:Nx,ntime)=0;
p(1:Ny,1:Nx,ntime)=0;

omegaz(1:Ny,1:Nx,ntime)=0;
%u(1,:,ntime)=1;
%define cylinder based on nodal points
mn(Ny,Nx)=0;
for i=1:Nx
  for j=1:Ny
    radius=sqrt((x(i)-xo)^2+(y(j)-yo)^2);
    if radius < rad
        mn(j,i)=kfunc;
        u(j,i)=0;
    end 
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define u staggered grid
xu=linspace(0,Lx,Nx);
yu=linspace(0-dy/2,Ly+dy/2,Ny+1);
us_star(1:Ny+1,1:Nx)=0;
us(1:Ny+1,1:Nx,ntime)=1;%initial condition

for j=1:Ny+1
    for i=1:Nx
    XU(j,i)=xu(i);
    YU(j,i)=yu(Ny+2-j);
    end
end
%define cylinder based on u staggered grid
mu(Ny+1,Nx)=0;
for i=1:Nx
    for j=1:Ny+1
    radius=sqrt((xu(i)-xo)^2+(yu(j)-yo)^2);
    if radius < rad
        mu(j,i)=kfunc;
        us(j,i)=0;
    end 
    end
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define v staggered grid
xv=linspace(0-dx/2,Lx+dx/2,Nx+1);
yv=linspace(0,Ly,Ny);
%[XV,YV]=meshgrid(xv,yv);
for j=1:Ny
    for i=1:Nx+1
    XV(j,i)=xv(i);
    YV(j,i)=yv(Ny+1-j);
    end
end

vs_star(1:Ny,1:Nx+1)=0;
vs(1:Ny,1:Nx+1,ntime)=0;
%define cylinder based on v staggered grid
mv(Ny,Nx+1)=0;
for i=1:Nx+1
    for j=1:Ny
    radius=sqrt((xv(i)-xo)^2+(yv(j)-yo)^2);
    if radius < rad
        mv(j,i)=kfunc;
        vs(j,i)=0;
    end 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define pressure and velocities at cell center of the grid 
xc=linspace(0-dx/2,Lx+dx/2,Nx+1);
yc=linspace(0-dy/2,Ly+dy/2,Ny+1);
%[XC,YC]=meshgrid(xc,yc);
for j=1:Ny+1
    for i=1:Nx+1
    XC(j,i)=xc(i);
    YC(j,i)=yc(Ny+2-j);
    end
end

pc(1:Ny+1,1:Nx+1,ntime)=0;
divu(1:Ny+1,1:Nx+1,ntime)=0;
uc(1:Ny+1,1:Nx+1,ntime)=0;
vc(1:Ny+1,1:Nx+1,ntime)=0;
%define cylinder based on cell center
mc(Ny+1,Nx+1)=0;
for i=1:Nx+1
    for j=1:Ny+1
    radius=sqrt((xc(i)-xo)^2+(yc(j)-yo)^2);
    if radius < rad
        mc(j,i)=kfunc;  
    end 
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time loop
for n=1:numel(t)
    tic
%**************% u momentum equation on u staggered grid %******************%
%advection term du^2/dx and duv/dy
[uc(:,:,n),vc(:,:,n),u(:,:,n),v(:,:,n)]=velocities(Nx,Ny,us,vs,n);

Hux(1:Ny+1,1:Nx)=0;
Huy(1:Ny+1,1:Nx)=0;
for i=2:Nx-1
    for j=2:Ny
    Hux(j,i)=-(1./dx).*(uc(j,i+1,n)^2-uc(j,i,n)^2);   
    Huy(j,i)=-(1./dy).*(u(j-1,i,n)*v(j-1,i,n)-u(j,i,n)*v(j,i,n));     
    end
end

%diffusion term in u momentum equation
alphau = 1+mu+(2*dt/Re)*(dx^(-2)+dy^(-2));
beta  = dt/Re/dx^2;
gamma = dt/Re/dy^2;

%RHS term in u momentum equation
RHSu(1:Ny+1,1:Nx)=0;
for i=2:Nx-1
    for j=2:Ny
        RHSu(j,i)=us(j,i,n)+dt*Hux(j,i)+dt*Huy(j,i);
    end
end
errora=1e-8;
us_star=ugauss_seidel(RHSu,Nx,Ny,alphau,beta,gamma,errora);
us_star=mass_correction(us_star,dy,1,Ly,Nx,Ny);

%**************% v momentum equation on v staggered grid %******************%
alphav = 1+mv+(2*dt/Re)*(dx^(-2)+dy^(-2));

%advection term dv^2/dy and duv/dx
Hvx(1:Ny,1:Nx+1)=0;
Hvy(1:Ny,1:Nx+1)=0;
for i=2:Nx
    for j=2:Ny-1
    Hvy(j,i)=-(1./dy).*(vc(j,i,n)^2-vc(j+1,i,n)^2);   
    Hvx(j,i)=-(1./dx).*(u(j,i,n)*v(j,i,n)-u(j,i-1,n)*v(j,i-1,n));     
    end
end

%RHS term in v momentum equation
RHSv(1:Ny,1:Nx+1)=0;
for i=2:Nx
    for j=2:Ny-1
        RHSv(j,i)=vs(j,i,n)+dt*Hvx(j,i)+dt*Hvy(j,i);
    end
end
errora=1e-8;
vs_star=vgauss_seidel(RHSv,Nx,Ny,alphav,beta,gamma,errora);

%**************% Pressure Poisson equation%******************%
alphap = 2*((dt/dx^2)+(dt/dy^2));
betap  = dt/dx^2;
gammap = dt/dy^2;

%divergence of ustar%
divustar(1:Ny+1,1:Nx+1)=0;
for i=2:Nx
    for j=2:Ny
        divustar(j,i)=((us_star(j,i)-us_star(j,i-1))/dx)+((vs_star(j-1,i)-vs_star(j,i))/dy);
    end
end

errora=1e-5;omg=1.98;
%[pc(:,:,n+1),residue]=pgauss_seidel(divustar,Nx,Ny,alphap,betap,gammap,errora,omg);
[pc(:,:,n+1),residue]=Copy_of_pgauss_seidel(divustar,Nx,Ny,alphap,betap,gammap,errora,omg);

for i=1:Nx
    for j=1:Ny
        p(j,i,n+1)=.25*(pc(j,i,n+1)+pc(j,i+1,n+1)+pc(j+1,i,n+1)+pc(j+1,i+1,n+1));
    end
end

%**************%velocity correction%******************%
for i=1:Nx+1
    for j=1:Ny
        vs(j,i,n+1)=vs_star(j,i)-(dt/dy)*(pc(j,i,n+1)-pc(j+1,i,n+1));
    end
end

for i=1:Nx
    for j=1:Ny+1
        us(j,i,n+1)=us_star(j,i)-(dt/dx)*(pc(j,i+1,n+1)-pc(j,i,n+1));
    end
end
[uc(:,:,n+1),vc(:,:,n+1),u(:,:,n+1),v(:,:,n+1)]=velocities(Nx,Ny,us,vs,n+1); 

%divergence of u
%divu(1:Ny,1:Nx,n+1)=0;
for i=2:Nx
    for j=2:Ny
        divu(j,i,n+1)=((us(j,i,n+1)-us(j,i-1,n+1))/dx)+((vs(j-1,i,n+1)-vs(j,i,n+1))/dy);
        %divu(j,i,n+1)=((u(j,i+1,n+1)-u(j,i-1,n+1))/2/dx)*1+1*((v(j-1,i,n+1)-v(j+1,i,n+1))/2/dy);
    end
end
%z vorticity
for i=1:Nx
    for j=1:Ny
        omegaz(j,i,n+1)=((vs(j,i+1,n+1)-vc(j,i,n+1))/dx)-((us(j+1,i,n+1)-us(j,i,n+1))/dy);
    end
end

%Display output on screen
maxu=max(max(u(:,:,n)));
maxv=max(max(v(:,:,n)));
cfl=max(maxu,maxv)*dt/dx;
Xd = ['iteration ',num2str(n) ,' time: ',num2str(t(n)),' CFL: ',num2str(cfl),' maxvel: ',num2str(max(maxu,maxv)),' divu: ',num2str(norm(divu(:,:,n)))];
disp(Xd)
MyVariable(1,1:numel(evalc('disp(Xd)')),n)=evalc('disp(Xd)');
toc
if mod(n,10000)==0
save('129x65_Re150_t0_1.mat','-v7.3')
end

end


%Uncomment the following lines for post processing.
% %%%%%%%%%%%%%%%% POST PROCESSING %%%%%%%%%%%%%%%%%%%%%%%
% 
% % residue
% for i=2:2202
% Ru(i-1)=rms(rms(u(:,:,i)-u(:,:,i-1)));
% end
% for i=2:2202
% Rv(i-1)=rms(rms(v(:,:,i)-v(:,:,i-1)));
% end
% 
% 
% 
% load('reference.mat')
% loc=41;file=2201;
% %u at x=0.5
% plot(Y(1:5:end,loc),u(1:5:end,loc,file),'r');hold on
% plot(ref(:,1),ref(:,2),'ko')
% 
% %percentage difference in u
% for i=2:16
%     pdiffu(i-1)=100*abs(ref(i,2)-u_vel(i))/abs(ref(i,3));
% end
% 
% %v at y=0.5
% plot(X(loc,1:5:end),v(loc,1:5:end,file),'r');hold on
% plot(ref(:,1),ref(:,3),'ko')
% 
% for i=2:16
%     pdiffv(i-1)=100*abs(ref(i,3)-v_vel(i))/abs(ref(i,3));
% end
% 
% %velocity contour
% %contourf(X,Y,p(:,:,2000),50,'linestyle','none');
% 
% [sx,sy] = meshgrid(-.2:.095:1,-.5:.095:1);
% streamline(stream2(X,Y,u(:,:,2201),v(:,:,2201),sx,sy));
% box on
% set(gca,'DataAspectratio',[1 1 1])
% xlabel('$x$','interpreter','latex','fontsize',18)
% ylabel('$y$','interpreter','latex','fontsize',18,'rot',0)
% h=title('Streamlines at $t=11$');
% set(h,'interpreter','latex','fontsize',14);
% set(gcf,'Color','w')
% %colorbar
% %colormap(bwr)
% set(gca,'XLim',[0 1])
% set(gca,'YLim',[0 1])
% set(gca,'fontsize',12,'fontname','times');
% 
% New_folder='/home/ghanesh/Desktop/Courses/export_fig-master';
% current_folder=pwd;
% filename='residue.eps';
% cd(New_folder)
% print2eps(filename,gcf)
% movefile(filename,current_folder)
% cd(current_folder)
% 
