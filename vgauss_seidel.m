function vstemp=vgauss_seidel(RHS,Nx,Ny,alpha,beta,gamma,errora)

v(1:Ny,1:Nx+1,1)=0;
vstemp(1:Ny,1:Nx+1)=0;
v(:,:,1)=.1;
error=errora;
E=1;k=1;counter2=0;
while E > error
counter2=counter2+1;
for i=2:Nx
  for j=2:Ny-1
      if mod(i+j,2)==0
      v(j,i,k+1)=(1/alpha(j,i))*(RHS(j,i)+beta*(v(j,i-1,k)+v(j,i+1,k))+gamma*(v(j+1,i,k)+v(j-1,i,k)));        
      end
      if mod(i+j,2)>0
      v(j,i,k+1)=(1/alpha(j,i))*(RHS(j,i)+beta*(v(j,i-1,k+1)+v(j,i+1,k+1))+gamma*(v(j+1,i,k+1)+v(j-1,i,k+1)));        
      end
  end
end
%Boundary conditions
v(1:Ny,1,k+1)=-v(1:Ny,2,k+1);
v(1:Ny,Nx+1,k+1)=v(1:Ny,Nx,k+1);
v(1,2:Nx,k+1)=0;
v(Ny,2:Nx,k+1)= 0;

X = ['v counter:',num2str(counter2), ' v error:',num2str(E)];
%disp(X)

E=norm(v(:,:,k+1)-v(:,:,k));


%pause
%close all
if E<error
vstemp=v(:,:,k+1);
end
v(:,:,k)=v(:,:,k+1);
end