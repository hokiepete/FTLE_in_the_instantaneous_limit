
clear all
close all
clc
filename='3D_OECS.avi'
v = VideoWriter(filename,'Uncompressed AVI');
v.FrameRate=10;
opengl('software')
open(v)
load OECS_DATA
xp=x;
yp=y;
zp=z;
d3ballsize=10
radius = 0.25
xorg= 1.2*pi
zorg= 1.75*pi
yorg= pi
[x,y,z]=sphere(d3ballsize);
x=radius*x;
y=radius*y;
z=radius*z;
x = reshape(x,[],1) + xorg;
y = reshape(y,[],1) + yorg;
z = reshape(z,[],1) + zorg;

%{
%{
dirdiv2 = dirdiv;
concav2 = concav;
s12 = s1;
for i=1:3
    dirdiv2 = smooth3(dirdiv2);
    concav2 = smooth3(concav2);
    s12 = smooth3(s12);
end
%}

figure
model.cdata = s1
model.alpha = [];
model.xdata = [0,2*pi];
model.ydata = [0,2*pi];
model.zdata = [0,2*pi];
model.parent = [];
model.handles = [];
model.texture = '3D';
colormap(parula);
vol3d_v2(model)
colorbar
title('s1')
view(3)
%

%

dirdiv(s1>0)=NaN;
dirdiv(concav<=0)=NaN;
dirdiv(abs(dirdiv)>0.1)=NaN;
%dirdiv = abs(dirdiv)<0.2;
figure
hold on
model.cdata = 1./(abs(dirdiv)+1)
model.alpha = [];%1./abs(dirdiv);
model.xdata = [0,2*pi];
model.ydata = [0,2*pi];
model.zdata = [0,2*pi];
model.parent = [];
model.handles = [];
model.texture = '3D';
colormap(parula);
vol3d_v2(model)
colorbar
title('dirdiv')
xlabel('x')
ylabel('y')
zlabel('z')
scatter3(x,y,z)
view(3)
%}


n=25;
tend=2*pi
twant = linspace(0,tend,n);
for i = 1:length(x)
    y0=[x(i),y(i),z(i)];
    [t,yout] = ode45(@abc_int,[0,tend],y0);
    xx(:,i) = interp1(t,yout(:,1),twant);
    yy(:,i) = interp1(t,yout(:,2),twant);
    zz(:,i) = interp1(t,yout(:,3),twant);
end

figure
hold on
for i =1:n
    scatter3(xx(i,:),yy(i,:),zz(i,:),'.');
    %axis([0,2*pi,0,2*pi,0,2*pi])
    view(3)
    %pause(0.5)
end
axis equal


%
for i=1:3
    dirdiv = smooth3(dirdiv);
    concav = smooth3(concav);
    s1 = smooth3(s1);
end
%}
dirdiv(s1>0)=NaN;
dirdiv(concav<=0)=NaN;
dirdiv(xp<2.5)=NaN;
dirdiv(xp>5)=NaN;
dirdiv(yp<2)=NaN;
dirdiv(yp>4.5)=NaN;
dirdiv(zp<4)=NaN;
dirdiv(abs(dirdiv)>0.2)=NaN;
FV=isosurface(xp,yp,zp,dirdiv,0);
%FV = smoothpatch(FV)


%{
figure
hold on
p=patch(FV,'facecolor','blue','edgecolor','none')
scatter3(x,y,z,'r')
axis([0,2*pi,0,2*pi,0,2*pi])
xlabel('x')
ylabel('y')
zlabel('z')
view(3)
camlight 
lighting gouraud
%}
%
figure
hold on
splitpatch = splitFV(FV);
index = 1;
len = length(splitpatch(1).faces)
for i=2:length(splitpatch)
    if length(splitpatch(i).faces)>len
        index = 1
        len = length(splitpatch(i).faces)
    end
end
FV=splitpatch(index)
patch(FV,'facecolor','blue','edgecolor','none');
scatter3(x,y,z,'r')
axis([0,2*pi,0,2*pi,0,2*pi])
xlabel('x')
ylabel('y')
zlabel('z')
camlight
lighting gouraud

points = FV.vertices;
n=25;
tend=2*pi
twant = linspace(0,tend,n);
for i = 1:length(points)
    y0=[points(i,1),points(i,2),points(i,3)];
    [t,yout] = ode45(@abc_int,[0,tend],y0);
    xvert(:,i) = interp1(t,yout(:,1),twant);
    yvert(:,i) = interp1(t,yout(:,2),twant);
    zvert(:,i) = interp1(t,yout(:,3),twant);
end
verts = cat(3,xvert,yvert,zvert);
fig=figure
for i =1:25
    clf
    hold on
    FV.vertices=squeeze(verts(i,:,:));
    patch(FV,'facecolor','blue','edgecolor','none');
    scatter3(xx(i,:),yy(i,:),zz(i,:),'r.');
    xlabel('x')
    ylabel('y')
    zlabel('z')
    camlight
    lighting gouraud
    axis([-15,5,-10,15,-5,15])
    %view(3)
    view(-124,32)
    pause(0.25)
     g=getframe(fig);
    writeVideo(v,g)
end
close(v)


%}
%{
dirdiv2(s12>0)=NaN;
dirdiv2(concav2<=0)=NaN;
dirdiv2(abs(dirdiv2)>0.2)=NaN;
figure
isosurface(x,y,z,dirdiv2,0)
view(3)
%}