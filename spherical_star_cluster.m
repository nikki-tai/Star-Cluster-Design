clear all
close all
G=1
R=1
a=1.5*R
n=1000
nred=50
mstar=1/n
m=mstar*ones(n,1);
x=zeros(n,1);
y=zeros(n,1);
z=zeros(n,1);
r=zeros(n,1);
s=zeros(n,1);
[A,map] = rgb2ind(frame2im(getframe),256);
imwrite(A,map,'starcluster.gif','LoopCount',65535,'DelayTime',0.01);
for i=1:n
    r(i)=2*R;
    while r(i) > R
        x(i) = R*2*(rand-0.5);
        y(i) = R*2*(rand-0.5);
        z(i) = R*2*(rand-0.5);
        r(i) = sqrt((x(i))^2 + (y(i))^2 + (z(i))^2);
    end
end
xcm=sum(x.*m)/sum(m);
ycm=sum(y.*m)/sum(m);
zcm=sum(z.*m)/sum(m);
x=x-xcm;
y=y-ycm;
z=z-zcm;

r=sqrt(x.^2 + y.^2 + z.^2);
mass = n * mstar * (r/R).^3;
s = sqrt (G * mass ./ r);

u = randn(n,1);
v = randn(n,1);
w = randn(n,1);

uu = y.*w - z.*v;
vv = z.*u - x.*w;
ww = x.*v - y.*u;

ss = sqrt(uu.^2 + vv.^2 + ww.^2);

u = s.*uu./ss;
v = s.*vv./ss;
w = s.*ww./ss;

ucm=sum(u.*m)/sum(m);
vcm=sum(v.*m)/sum(m);
wcm=sum(w.*m)/sum(m);

u = u - ucm;
v = v - vcm;
w = w - wcm;

hxyz=plot3(x,y,z,'b.')
hold on
hxyzred=plot3(x(1:nred),y(1:nred),z(1:nred),'ro')
hxyztrj= ... 
    plot3([x(1:nred),x(1:nred)]' ...
         ,[y(1:nred),y(1:nred)]' ...
         ,[z(1:nred),z(1:nred)]')
hold off
axis([-a,a,-a,a,-a,a])
axis equal
axis manual
tmax=20
clockmax = 1000
dt = tmax/clockmax
dtGm = dt*G*m;
xsave=zeros(nred,clockmax);
ysave=zeros(nred,clockmax);
zsave=zeros(nred,clockmax);
for clock=1:clockmax
    for j=1:n
        dx = x(j) - x;
        dy = y(j) - y;
        dz = z(j) - z;
        rr=sqrt(dx.^2 + dy.^2 + dz.^2);
        rr(j)=1;
        dxorr3 = dx ./ rr.^3;
        dyorr3 = dy ./ rr.^3;
        dzorr3 = dz ./ rr.^3;
        u = u + dtGm(j)*dxorr3;
        v = v + dtGm(j)*dyorr3;
        w = w + dtGm(j)*dzorr3;
    end
    x = x + dt * u;
    y = y + dt * v;
    z = z + dt * w;
    xsave(:,clock)=x(1:nred);
    ysave(:,clock)=y(1:nred);
    zsave(:,clock)=z(1:nred);    
    hxyz.XData = x; hxyzred.XData = x(1:nred);
    hxyz.YData = y; hxyzred.YData = y(1:nred);
    hxyz.ZData = z; hxyzred.ZData = z(1:nred);
    for i=1:nred
         hxyztrj(i).XData = xsave(i,1:clock);
         hxyztrj(i).YData = ysave(i,1:clock);
         hxyztrj(i).ZData = zsave(i,1:clock);
    end
    drawnow
    if(mod(j,20)==0)
        [A,map] = rgb2ind(frame2im(getframe),256);
        imwrite(A,map,'2.gif','WriteMode','append','DelayTime',0.01);
    end
end



