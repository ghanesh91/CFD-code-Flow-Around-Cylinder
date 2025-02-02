function [pc,residue]=Copy_of_pgauss_seidel(divustar,Nx,Ny,alphap,betap,gammap,errora,omg)

diff(1:Ny+1,1:Nx+1)=0;
p(1:Ny+1,1:Nx+1,1)=0;
pc(1:Ny+1,1:Nx+1)=0;
p(:,:,1)=1;

E=1;k=1;counter3=0;
while E > errora
counter3=counter3+1;
for i=2:Nx
  for j=2:Ny      
       if mod(i+j,2)==0
      p(j,i,k+1)=omg*(1/alphap)*(-divustar(j,i)+betap*(p(j,i-1,k)+p(j,i+1,k))+gammap*(p(j+1,i,k)+p(j-1,i,k)))...
          +(1-omg)*p(j,i,k);        
       end
  end
end
 for i=2:Nx
   for j=2:Ny       
       if mod(i+j,2)>0
       p(j,i,k+1)=omg*(1/alphap)*(-divustar(j,i)+betap*(p(j,i-1,k+1)+p(j,i+1,k+1))+gammap*(p(j+1,i,k+1)+p(j-1,i,k+1)))...
           +(1-omg)*p(j,i,k);        
       end
   end
 end
%Pressure boundary conditions
p(2:Ny,1,k+1)=p(2:Ny,2,k+1);
p(2:Ny,Nx+1,k+1)=p(2:Ny,Nx,k+1);
p(Ny+1,2:Nx,k+1)  = p(Ny,2:Nx,k+1);
p(1,2:Nx,k+1)=p(2,2:Nx,k+1);
% p(1,1,k+1)=p(2,2,k+1);
% p(1,Nx+1,k+1)=p(2,Nx,k+1);
% p(Ny+1,1,k+1)=p(Ny,2,k+1);
% p(Ny+1,Nx+1,k+1)=p(Ny,Nx,k+1);
for i=2:Nx
     for j=2:Ny
     %diff(j,i)=abs(p(j,i,k+1)*alphap-(-divustar(j,i)+betap*(p(j,i-1,k+1)+p(j,i+1,k+1))+gammap*(p(j+1,i,k+1)+p(j-1,i,k+1))));
     diff(j,i)=abs(p(j,i,k+1)*alphap-(-divustar(j,i)+betap*(p(j,i-1,k+1)+p(j,i+1,k+1))+gammap*(p(j+1,i,k+1)+p(j-1,i,k+1))));
     end
end

E=norm(diff);
X = ['p counter:',num2str(counter3),' p error:',num2str(E)];
%disp(X)
%contourf(p(:,:,2),22)
residue(counter3)=E;
if E< errora
    pc=p(:,:,k+1);
end
p(:,:,k)=p(:,:,k+1);
end

