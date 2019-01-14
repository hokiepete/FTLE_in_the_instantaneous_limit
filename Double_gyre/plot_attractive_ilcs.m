
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

%
for i=1:3
    dirdiv = smooth3(dirdiv);
    concav = smooth3(concav);
    s1 = smooth3(s1);
end
%}
dirdiv(s1>0)=NaN;
dirdiv(concav<=0)=NaN;
%dirdiv(xp<2.5)=NaN;
%dirdiv(xp>5)=NaN;
%dirdiv(yp<2)=NaN;
%dirdiv(yp>4.5)=NaN;
%dirdiv(zp<4)=NaN;
dirdiv(abs(dirdiv)>0.2)=NaN;
FV=isosurface(xp,yp,zp,dirdiv,0);
az = 47
el = 35
fig=figure('units','inch','position',[0,0,6,6]);
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
axis([0,2*pi,0,2*pi,0,2*pi])
xlabel('x')
ylabel('y')
zlabel('z')
view(az,el)
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
for i =1:n
    clf
    hold on
    FV.vertices=squeeze(verts(i,:,:));
    patch(FV,'facecolor','blue','edgecolor','none');
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
    pause(0.25)
    g=getframe(fig);
    writeVideo(v,g)
end
close(v)


%}


