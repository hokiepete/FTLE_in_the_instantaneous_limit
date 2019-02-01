
clear all
close all
clc
filename='attracting_ilcs.avi'
v = VideoWriter(filename,'Uncompressed AVI');
v.FrameRate=10;
opengl('software')
open(v)
load OECS_DATA
xp=x;
yp=y;
zp=z;
d3ballsize=151
radius = 0.1
xorg= 1
zorg= 0.5
yorg= 0.8
zorg2= zorg%-2*radius%3.3
[x,y,z]=sphere(d3ballsize);
x=radius*x;
y=radius*y;
z=radius*z;
x1 = reshape(x,[],1) + xorg;%+radius;
y1 = reshape(y,[],1) + yorg;
z1 = reshape(z,[],1) + zorg;

%x1 = vertcat(x1,reshape(x,[],1) + xorg-radius);
%y1 = vertcat(y1,reshape(y,[],1) + yorg);
%z1 = vertcat(z1,reshape(z,[],1) + zorg);


lims = [0,0,0;
        0,1,0;
        0,1,1;
        2,0,0;
        2,0,1;
        2,1,0
];
%0,0,1;
%2,1,1
n=141;
tend=8
twant = linspace(0,tend,n);
for i = 1:length(x1)
    y01=[x1(i),y1(i),z1(i)];
    [t1,yout1] = ode45(@dg_int,[0,tend],y01);
    xx1(:,i) = interp1(t1,yout1(:,1),twant);
    yy1(:,i) = interp1(t1,yout1(:,2),twant);
    zz1(:,i) = interp1(t1,yout1(:,3),twant);
end

%
for i=1:3
    dirdiv1 = smooth3(dirdiv1);
    concav1 = smooth3(concav1);
    s1 = smooth3(s1);
end
%}
dirdiv1(s1>0)=NaN;
dirdiv1(concav1<=0)=NaN;
%dirdivn(abs(dirdivn)>0.2)=NaN;
%

%{
dirdivn(yp<2.5)=NaN;
%dirdivn(xp>3.5)=NaN;
dirdivn(xp<3)=NaN;
dirdivn(zp>5)=NaN;
dirdivn(zp<3)=NaN;
FV=isosurface(xp,yp,zp,dirdivn,0);
%FV = smoothpatch(FV)
%}
%{
splitpatch = splitFV(FV);
index = 1;
len = length(splitpatch(1).faces)
for i=2:length(splitpatch)
    if length(splitpatch(i).faces)>len
        index = i
        len = length(splitpatch(i).faces)
    end
end
FV=splitpatch(index)
%}
FV=isosurface(xp,yp,zp,dirdiv1,0);
points = FV.vertices;
%n=25;
twant = linspace(0,tend,n);
for i = 1:length(points)
    y0=[points(i,1),points(i,2),points(i,3)];
    [t,yout] = ode45(@dg_int,[0,tend],y0);
    xvert(:,i) = interp1(t,yout(:,1),twant);
    yvert(:,i) = interp1(t,yout(:,2),twant);
    zvert(:,i) = interp1(t,yout(:,3),twant);
end
verts = cat(3,xvert,yvert,zvert);
%
fig=figure
for i =1:n
    clf
    hold on
    FV.vertices=squeeze(verts(i,:,:));
    patch(FV,'facecolor','blue','edgecolor','none','FaceAlpha',0.4);
    scatter3(lims(:,1),lims(:,2),lims(:,3),1,'w.');
    %scatter3(xx1(i,:),yy1(i,:),zz1(i,:),'g.');
    dt1 = delaunayTriangulation(xx1(i,:)',yy1(i,:)',zz1(i,:)');
    triboundary1 = convexHull(dt1);
    trisurf(triboundary1,xx1(i,:)',yy1(i,:)',zz1(i,:)','FaceColor', 'green','edgecolor','none','FaceAlpha',1.0);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title(sprintf('time = %d',twant(i)))
    camlight
    lighting gouraud
    axis equal tight
    view(3)
    g=getframe(fig);
    writeVideo(v,g)
end

close(v)
%{
az2 = 120%28%47
el2 = 15%29%35
az = -60%-190%-25
el =15% 21
alpT = 1.0
alpS = 0.6
font = 'cmr'
width = 5+3/8
fig=figure('units','inch','position',[0,0,width,width],'DefaultTextFontName', font, 'DefaultAxesFontName', font);
i = 1
subplot(2,2,1)
hold on
FV.vertices=squeeze(verts(i,:,:));
patch(FV,'facecolor','blue','edgecolor','none','FaceAlpha',alpS);
scatter3(lims(:,1),lims(:,2),lims(:,3),1,'w.');
dt1 = delaunayTriangulation(xx1(i,:)',yy1(i,:)',zz1(i,:)');
triboundary1 = convexHull(dt1);
trisurf(triboundary1,xx1(i,:)',yy1(i,:)',zz1(i,:)','FaceColor', 'green','edgecolor','none','FaceAlpha',1.0)
xlabel('x')
ylabel('y')
zlabel('z')
%set(gca, 'XTick',[2,3])
%set(gca, 'YTick',[3,4])
camlight
lighting gouraud
%axis([0,2,0,1,0,1])
axis equal tight;
view(az,el)

subplot(2,2,2)
hold on
FV.vertices=squeeze(verts(i,:,:));
patch(FV,'facecolor','blue','edgecolor','none','FaceAlpha',alpS);
scatter3(lims(:,1),lims(:,2),lims(:,3),1,'w.');
dt1 = delaunayTriangulation(xx1(i,:)',yy1(i,:)',zz1(i,:)');
triboundary1 = convexHull(dt1);
trisurf(triboundary1,xx1(i,:)',yy1(i,:)',zz1(i,:)','FaceColor', 'green','edgecolor','none','FaceAlpha',1.0)
xlabel('x')
ylabel('y')
zlabel('z')
%set(gca, 'XTick',[2,3])
%set(gca, 'YTick',[3,4])
camlight
lighting gouraud
axis equal tight;
view(az2,el2)

i = 4
subplot(2,2,3)
hold on
FV.vertices=squeeze(verts(i,:,:));
patch(FV,'facecolor','blue','edgecolor','none','FaceAlpha',alpS);
scatter3(lims(:,1),lims(:,2),lims(:,3),1,'w.');
dt1 = delaunayTriangulation(xx1(i,:)',yy1(i,:)',zz1(i,:)');
triboundary1 = convexHull(dt1);
trisurf(triboundary1,xx1(i,:)',yy1(i,:)',zz1(i,:)','FaceColor', 'green','edgecolor','none','FaceAlpha',1.0)
xlabel('x')
ylabel('y')
zlabel('z')
camlight
lighting gouraud
%axis([0,2,0,1,0,1])
axis equal tight;
view(az,el)

subplot(2,2,4)
hold on
FV.vertices=squeeze(verts(i,:,:));
patch(FV,'facecolor','blue','edgecolor','none','FaceAlpha',alpS);
scatter3(lims(:,1),lims(:,2),lims(:,3),1,'w.');
dt1 = delaunayTriangulation(xx1(i,:)',yy1(i,:)',zz1(i,:)');
triboundary1 = convexHull(dt1);
trisurf(triboundary1,xx1(i,:)',yy1(i,:)',zz1(i,:)','FaceColor', 'green','edgecolor','none','FaceAlpha',1.0)
xlabel('x')
ylabel('y')
zlabel('z')
camlight
lighting gouraud
%axis([0,2,0,1,0,1])
axis equal tight;
view(az2,el2)

saveas(fig,'attracting_ilcs_v4.eps','epsc')
%}
%{
i = n
subplot(3,2,5)
hold on
FV.vertices=squeeze(verts(i,:,:));
patch(FV,'facecolor','blue','edgecolor','none');
dt = delaunayTriangulation(xx(i,:)',yy(i,:)',zz(i,:)');
triboundary = convexHull(dt);
trisurf(triboundary,xx(i,:)',yy(i,:)',zz(i,:)','FaceColor', 'green','edgecolor','none','FaceAlpha',1.0)
xlabel('x')
ylabel('y')
zlabel('z')
title(sprintf('time = %f',twant(i)))
camlight
lighting gouraud
%axis([-15,5,-10,15,-5,15])
axis equal tight
%view(3)
view(-124,32)
%view(20,20)

subplot(3,2,6)
dt = delaunayTriangulation(xx(i,:)',yy(i,:)',zz(i,:)');
triboundary = convexHull(dt);
trisurf(triboundary,xx(i,:)',yy(i,:)',zz(i,:)','FaceColor', 'green','edgecolor','none','FaceAlpha',1.0)
xlabel('x')
ylabel('y')
zlabel('z')
title(sprintf('time = %f',twant(i)))
camlight
lighting gouraud
%axis([-15,5,-10,15,-5,15])
axis equal tight
%view(3)
view(-124,32)
%view(20,20)


%}


