function [uc,vc,u,v]=velocities(Nx,Ny,us,vs,n)
uc(Ny+1,Nx+1)=0;
vc(Ny+1,Nx+1)=0;
u(Ny,Nx)=0;
v(Ny,Nx)=0;


for j=2:Ny
  for i=2:Nx
     uc(j,i)=(us(j,i-1,n)+us(j,i,n))/2;% centerline u velocity interpolated from adjacent staggered u vlocities
     vc(j,i)=(vs(j-1,i,n)+vs(j,i,n))/2;% centerline v velocity interpolated from adjacent staggered v vlocities  
  end
end

uc(2:Ny,1)=2-uc(2:Ny,2);
uc(2:Ny,Nx+1)=uc(2:Ny,Nx);
uc(1,2:Nx)=uc(2,2:Nx);
uc(Ny+1,2:Nx)=uc(Ny,2:Nx);

vc(2:Ny,1)=-vc(2:Ny,2);
vc(2:Ny,Nx+1)=vc(2:Ny,Nx);
vc(1,2:Nx)=-vc(2,2:Nx);
vc(Ny+1,2:Nx)=-vc(Ny,2:Nx);

% uc(1,1)=1.5;%uc(2,2);
% vc(1,1)=0;%vc(2,2);
% uc(Ny+1,1)=1.5;%uc(Ny,2);
% vc(Ny+1,1)=0;%vc(Ny,2);
% uc(1,Nx+1)=0;%uc(2,Nx);
% vc(1,Nx+1)=0;%vc(2,Nx);
% uc(Ny+1,Nx+1)=0;%uc(Ny,Nx);
% vc(Ny+1,Nx+1)=0;%vc(Ny,Nx);


for i=1:Nx
    for j=1:Ny
     v(j,i)=(vs(j,i+1,n)+vs(j,i,n))/2;
     u(j,i)=(us(j+1,i,n)+us(j,i,n))/2;
    end
end

end