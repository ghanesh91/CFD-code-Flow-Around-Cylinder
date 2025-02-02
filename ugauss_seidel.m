function ustemp=ugauss_seidel(RHS,Nx,Ny,alpha,beta,gamma,errora)

u(1:Ny+1,1:Nx,1)=0;
ustemp(1:Ny+1,1:Nx)=0;
u(:,:,1)=.1;
error=errora;
E=1;counter1=0;k=1;
while E > error
counter1=counter1+1;
for i=2:Nx-1
  for j=2:Ny
      if mod(i+j,2)==0
      u(j,i,k+1)=(1/alpha(j,i))*(RHS(j,i)+beta*(u(j,i+1,k)+u(j,i-1,k))+gamma*(u(j+1,i,k)+u(j-1,i,k)));        
      end
      if mod(i+j,2)>0
      u(j,i,k+1)=(1/alpha(j,i))*(RHS(j,i)+beta*(u(j,i+1,k+1)+u(j,i-1,k+1))+gamma*(u(j+1,i,k)+u(j-1,i,k+1)));        
      end      
  end
end

%u velocity boundary conditions
u(2:Ny,1,k+1)=1;
u(2:Ny,Nx,k+1)=u(2:Ny,Nx-1,k+1);%(4*u(2:Ny,Nx-1,k+1)-u(2:Ny,Nx-2,k+1))/3;
u(Ny+1,1:Nx,k+1)=u(Ny,1:Nx,k+1);
u(1,1:Nx,k+1)=u(2,1:Nx,k+1);
E=norm(u(:,:,k+1)-u(:,:,k));

X = ['u counter:',num2str(counter1), ' u error:',num2str(E)];
%disp(X)
%pause
%close all
if E<error
ustemp=u(:,:,k+1);
end
u(:,:,k)=u(:,:,k+1);
end