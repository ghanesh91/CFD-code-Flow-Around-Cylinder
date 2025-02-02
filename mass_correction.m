function [us]=mass_correction(us,dy,n,Ly,Nx,Ny)

massin=0;massout=0;
for j=2:Ny    
massin=massin+abs(us(j,1,n))*dy;
end

for j=2:Ny
massout=massout+abs(us(j,Nx,n))*dy;
end
deltau=(massout-massin)/(Ly);
deltam=massout-massin;
ratio=massin/massout;

X1 = ['before mass correction dm: ',num2str(deltam)];
disp(X1)
%  for j=2:Ny
%  us(j,Nx,n)=us(j,Nx,n)*ratio;
%  end

for j=2:Ny
if massin>massout    
us(j,Nx,n)=us(j,Nx,n)+(massin-massout)/Ly;
else
us(j,Nx,n)=us(j,Nx,n)-(massout-massin)/Ly;
end
end

% for j=2:Ny
% deltau(j)=us(j,Nx,n)-us(j,1,n);
% if deltau(j) > 0
% us(j,Nx,n)=us(j,Nx,n)-abs(deltau(j));
% else
% us(j,Nx,n)=us(j,Nx,n)+abs(deltau(j));    
% end
% end



us(:,Nx-1,n)=us(:,Nx,n);
massin=0;massout=0;
for j=2:Ny
massin=massin+abs(us(j,1,n))*dy;
end

for j=2:Ny
massout=massout+abs(us(j,Nx,n))*dy;
end
deltau=(massout-massin)/(Ly);
deltam=massout-massin;

X2 = ['after mass correction dm: ',num2str(deltam), ' ratio: ', num2str(ratio)];
disp(X2)

end