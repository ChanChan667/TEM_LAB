function [monoSEobjectfunction]=MonoSEobjectfunction(Z,theta_0,params,Nx,Ny,Lx,Ly,r0)
% generate SE object function from Zhu's Ultramicroscope paper
% z is atomic number
E=6.75*Z;%energy loss in ev,E=6.75*z z is atomic number
lambda=HighEnergyWavLen_X(params.KeV);
k=2*pi/lambda;
thetaE=(params.KeV+511)/(params.KeV+1022)*E/params.KeV/1000;
dx=Lx/Nx;
dy=Ly/Ny;
x=-Lx/2:dx:Lx/2-dx;
y=-Ly/2:dy:Ly/2-dy;
[X,Y]=meshgrid(x,y);
fx=InitFreqAxis(Lx, Nx);
fy=InitFreqAxis(Ly, Ny);
[Fx,Fy]=meshgrid(fx,fy);

% %using numerical integral to calculate A
% r=(X-r0(1)).^2+(Y-r0(2)).^2;
% fun=@(theta,c) besselj(0,k*c*sin(theta)).*sin(2*theta).*...
%     sqrt(thetaE^2./(thetaE^2+theta.^2)./(theta.^2+theta_0^2));
% A=zeros(Nx,Ny);
% for i=1:Nx
%     for j=1:Ny
%         A(i,j)=integral(@(theta) fun(theta,r(i,j)),0,pi);
%     end
% end
% monoSEobjectfunction=(k^2/(4*pi)*A).^2;

% using fft to calculate A
q=sqrt(Fx.^2+Fy.^2);
theta=asin(q/k);
A=sqrt(thetaE^2./(thetaE^2+theta.^2)./(theta.^2+theta_0^2));%.*...
    %CircApert_X(Lx, Ly, Nx, Ny, lambda,2*pi*1e3);
A=ifft2(fftshift(A.*exp(1i*(Fx*r0(1)+Fy*r0(2)))));
monoSEobjectfunction=abs(ifftshift(A));




end