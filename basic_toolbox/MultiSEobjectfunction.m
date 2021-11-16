function [multiSEobjectfunction]=MultiSEobjectfunction(theta_0,params,Nx,Ny,Lx,Ly,fracTypeCoord, expanNum, lattConst)
dx=Lx/Nx;
dy=Ly/Ny;
x=-Lx/2:dx:Lx/2-dx;
y=-Ly/2:dy:Ly/2-dy;
[X,Y]=meshgrid(x,y);
fx=InitFreqAxis(Lx, Nx);
fy=InitFreqAxis(Ly, Ny);
[Fx,Fy]=meshgrid(fx,fy);
Z=fracTypeCoord(1,:);
atomnumber=length(Z);
multiSEobjectfunction=zeros(Nx,Ny);
a=lattConst(1);
b=lattConst(2);
for n=1:atomnumber
    for i= 1:expanNum(1)+1
        for j=1:expanNum(2)+1 %+1是为了防止周围空白一行
            rx=a*fracTypeCoord(3,n)+(i-1)*a;
            ry=b*fracTypeCoord(4,n)+(j-1)*b;
            r0=[rx,ry];
            A=MonoSEobjectfunction(Z(n),theta_0,params,Nx,Ny,Lx,Ly,r0);
            multiSEobjectfunction=multiSEobjectfunction+A;
        end
    end
end
figure
imagesc(multiSEobjectfunction)

end