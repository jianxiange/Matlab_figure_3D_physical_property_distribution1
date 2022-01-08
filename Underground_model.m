% Draw the 3D underground model.
% Developed by Xiange Jian. China University og Geosciences (Wuhan).
% Email: jian_xiange@cug.edu.cn
% Last revision in 25/12/2021.
%% Input Data:
% C: A three-dimensional matrix corresponding the physical parameter to each underground grid.
% xmin and xmax: Underground grid x coordinate range.
% ymin and ymax: Underground grid y coordinate range.
% zmin and zmax: Underground grid z coordinate range.
% The size of each underground grid is equal.
%% Input Data:
clear;close all;clc

DATA1=load('model3.dat');
DATA2=load('deltaT.dat');

deltaT = reshape(DATA2(:,4),21,21);  % the size of ¦¤T
xmin = 0; xmax = 1000;               % underground x limitation
ymin = 0; ymax = 1000;               % underground y limitation
zmin = 0; zmax =  500;               % underground z limitation

xn = DATA1(:,1);                     % underground x coordinate 
yn = DATA1(:,2);                     % underground y coordinate
zn = -DATA1(:,3);                    % z coordinate (upward positive)
C = reshape(DATA1(:,4),20,20,10);    % the size of underground grid
CC = DATA1(:,5);                     % the invesion data
%%
[nx,ny,nz] = size(C);
a = (xmax-xmin)/(nx);
b = (ymax-ymin)/(ny);
c = (zmax-zmin)/(nz);
C = C(:);
set (gcf,'Position',[00,00,1200,900], 'color','w')
Mmin = min(C);
Mmax = max(C);
limitC = 0.2*Mmax;                   % data under LimitC wolld not be shown
Azimuth=45;
Elevation=15;
%% Forward model and corresponding anomaly
subplot(2,2,1);
hold on;
for i = 1:length(xn)
    if CC(i) < limitC
        continue
    end
    
    [X,Y] = meshgrid(xn(i)-a/2:a:xn(i)+a/2,yn(i)-b/2:b:yn(i)+b/2);
    Z = ones(size(X))*(zn(i)-c/2);
    patch([X(1) X(2) X(4) X(3)],[Y(1) Y(2) Y(4) Y(3)],Z(:)  ,CC(i))
    patch([X(1) X(2) X(4) X(3)],[Y(1) Y(2) Y(4) Y(3)],Z(:)+c,CC(i))
    
    [X,Z] = meshgrid(xn(i)-a/2:a:xn(i)+a/2,zn(i)-c/2:c:zn(i)+c/2);
    Y = ones(size(X))*(yn(i)-b/2);
    patch([X(1) X(2) X(4) X(3)],Y(:)  ,[Z(1) Z(2) Z(4) Z(3)],CC(i))
    patch([X(1) X(2) X(4) X(3)],Y(:)+b,[Z(1) Z(2) Z(4) Z(3)],CC(i))
    
    [Y,Z] = meshgrid(yn(i)-b/2:b:yn(i)+b/2,zn(i)-c/2:c:zn(i)+c/2);
    X = ones(size(X))*(xn(i)-a/2);
    patch(X(:)  ,[Y(1) Y(2) Y(4) Y(3)],[Z(1) Z(2) Z(4) Z(3)],CC(i))
    patch(X(:)+a,[Y(1) Y(2) Y(4) Y(3)],[Z(1) Z(2) Z(4) Z(3)],CC(i))
    
end
cl1 = colorbar;
cl1.Label.String = 'Magnetization (A/m)';
set(cl1,'fontsize',10,'fontname','Times');
ax1 = gca;
axpos = ax1.Position;
cl1.Position(1) = 1.04*cl1.Position(1);% x offset of colorbar
cl1.Position(2) = 1.20*cl1.Position(2);% y offset of colorbar
cl1.Position(3) = 0.50*cl1.Position(3);% width of colorbar
cl1.Position(4) = 0.35*cl1.Position(4);% length of colorbar
ax1.Position = axpos;

colormap(ax1,jet);
caxis([Mmin,Mmax])

deltaT1 = deltaT*(Mmax-Mmin)/(max(max(deltaT))-min(min(deltaT)));
deltaT1 = deltaT1+Mmax-max(max(deltaT1));
[xx,yy]=meshgrid(xmin:a:xmax,ymin:b:ymax);
total=contourf(xx,yy,fliplr(deltaT1),10);
alpha(0.7)

grid on
axis equal
axis([xmin xmax ymin ymax zmin zmax])
xlabel('Easting(m)','rotation',0);
ylabel('Northing(m)','rotation',0);
zlabel('Depth(m)','rotation',90);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1,'LineWidth',0.5);
set(gca,'fontsize',10,'fontname','Times');
set(gca,'ZDir','reverse')
view([Azimuth,Elevation]+10)

% set the label to be parallel to the coordinate axis
h_a = gca;            
[thx,thy,thz] = axislabel_rotation_angle(h_a);
set(get(h_a,'xlabel'),'rotation',thx(1)); 
% if the label is upside down,plus 180 degree :thx(1)+180
set(get(h_a,'ylabel'),'rotation',thy(1));
set(get(h_a,'zlabel'),'rotation',thz(1));
%% Slice
f2 = subplot(2,2,2);

xlice = 10;   % x slice position
ylice = 10;   % x slice position
zlice = 5;    % x slice position

Xlicen=find(xn==xlice*a-a/2);
Ylicen=find(yn==ylice*b-b/2);
Zlicen=find(zn==zlice*c-c/2);

for n = 1:length(Xlicen)
    i=Xlicen(n);
    [Y,Z] = meshgrid(yn(i)-b/2:b:yn(i)+b/2,zn(i)-c/2:c:zn(i)+c/2);
    X = ones(size(Y))*xn(i);
    patch(X(:),[Y(1) Y(2) Y(4) Y(3)],[Z(1) Z(2) Z(4) Z(3)],C(i),'EdgeColor','k','LineWidth',1)
end

for n = 1:length(Ylicen)
    i=Ylicen(n);
    [X,Z] = meshgrid(xn(i)-a/2:a:xn(i)+a/2,zn(i)-c/2:c:zn(i)+c/2);
    Y = ones(size(X))*yn(i);
    patch([X(1) X(2) X(4) X(3)],Y(:),[Z(1) Z(2) Z(4) Z(3)],C(i),'EdgeColor','k','LineWidth',1)
end

for n = 1:length(Zlicen)
    i=Zlicen(n);
    [X,Y] = meshgrid(xn(i)-a/2:a:xn(i)+a/2,yn(i)-b/2:b:yn(i)+b/2);
    Z = ones(size(X))*zn(i);
    patch([X(1) X(2) X(4) X(3)],[Y(1) Y(2) Y(4) Y(3)],Z(:),C(i),'EdgeColor','k','LineWidth',1)
end

alpha(f2,1)

cl1 = colorbar;
cl1.Label.String = 'Magnetization (A/m)';
set(cl1,'fontsize',10,'fontname','Times');
ax = gca;
axpos = ax.Position;
cl1.Position(1) = 1.04*cl1.Position(1);% x offset
cl1.Position(2) = 1.20*cl1.Position(2);% y offset
cl1.Position(3) = 0.50*cl1.Position(3);% width
cl1.Position(4) = 0.35*cl1.Position(4);% length
ax.Position = axpos;

colormap jet
caxis([Mmin,Mmax])
grid on
axis equal
axis([xmin xmax ymin ymax zmin zmax])
xlabel('Easting(m)','rotation',0);
ylabel('Northing(m)','rotation',0);
zlabel('Depth(m)','rotation',90);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1,'LineWidth',0.5);
set(gca,'fontsize',10,'fontname','Times');
set(gca,'ZDir','reverse')
view([Azimuth,Elevation]+10)

h_a = gca;
[thx,thy,thz] = axislabel_rotation_angle(h_a);
set(get(h_a,'xlabel'),'rotation',thx(1));
set(get(h_a,'ylabel'),'rotation',thy(1));
set(get(h_a,'zlabel'),'rotation',thz(1));
%% inversed model
subplot(2,2,3);
for i = 1:length(xn)
    if C(i) < limitC
        continue
    end
    
    [X,Y] = meshgrid(xn(i)-a/2:a:xn(i)+a/2,yn(i)-b/2:b:yn(i)+b/2);
    Z = ones(size(X))*(zn(i)-c/2);
    patch([X(1) X(2) X(4) X(3)],[Y(1) Y(2) Y(4) Y(3)],Z(:)  ,C(i))
    patch([X(1) X(2) X(4) X(3)],[Y(1) Y(2) Y(4) Y(3)],Z(:)+c,C(i))
    
    [X,Z] = meshgrid(xn(i)-a/2:a:xn(i)+a/2,zn(i)-c/2:c:zn(i)+c/2);
    Y = ones(size(X))*(yn(i)-b/2);
    patch([X(1) X(2) X(4) X(3)],Y(:)  ,[Z(1) Z(2) Z(4) Z(3)],C(i))
    patch([X(1) X(2) X(4) X(3)],Y(:)+b,[Z(1) Z(2) Z(4) Z(3)],C(i))
    
    [Y,Z] = meshgrid(yn(i)-b/2:b:yn(i)+b/2,zn(i)-c/2:c:zn(i)+c/2);
    X = ones(size(X))*(xn(i)-a/2);
    patch(X(:)  ,[Y(1) Y(2) Y(4) Y(3)],[Z(1) Z(2) Z(4) Z(3)],C(i))
    patch(X(:)+a,[Y(1) Y(2) Y(4) Y(3)],[Z(1) Z(2) Z(4) Z(3)],C(i))
    
end
alpha(0.8)
cl2=colorbar;
cl2.Label.String = 'Magnetization (A/m)';
set(cl2,'fontsize',10,'fontname','Times');
ax1 = gca;
axpos = ax1.Position;
cl2.Position(1) = 1.04*cl2.Position(1);% x offset
cl2.Position(2) = 2.0*cl2.Position(2);% y offset
cl2.Position(3) = 0.5*cl2.Position(3);% width
cl2.Position(4) = 0.35*cl2.Position(4);% length
ax1.Position = axpos;

colormap jet
caxis([Mmin,Mmax])
grid on
axis equal
axis([xmin xmax ymin ymax zmin zmax])
xlabel('Easting(m)','rotation',0);
ylabel('Northing(m)','rotation',0);
zlabel('Depth(m)','rotation',90);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1,'LineWidth',0.5);
set(gca,'fontsize',10,'fontname','Times');
set(gca,'ZDir','reverse')
view([Azimuth,Elevation+10])
cmap = colormap;

h_a = gca;
[thx,thy,thz] = axislabel_rotation_angle(h_a);
set(get(h_a,'xlabel'),'rotation',thx(1));
set(get(h_a,'ylabel'),'rotation',thy(1));
set(get(h_a,'zlabel'),'rotation',thz(1));
%% inversied vector
subplot(2,2,4);
a = a/20;
b = b/20;
c = c/20;
for i = 1:length(xn)
    if C(i) < limitC
        continue
    end
    
    [X,Y] = meshgrid(xn(i)-a/2:a:xn(i)+a/2,yn(i)-b/2:b:yn(i)+b/2);
    Z = ones(size(X))*(zn(i)-c/2);
    patch([X(1) X(2) X(4) X(3)],[Y(1) Y(2) Y(4) Y(3)],Z(:)  ,C(i))
    patch([X(1) X(2) X(4) X(3)],[Y(1) Y(2) Y(4) Y(3)],Z(:)+c,C(i))
    
    [X,Z] = meshgrid(xn(i)-a/2:a:xn(i)+a/2,zn(i)-c/2:c:zn(i)+c/2);
    Y = ones(size(X))*(yn(i)-b/2);
    patch([X(1) X(2) X(4) X(3)],Y(:)  ,[Z(1) Z(2) Z(4) Z(3)],C(i))
    patch([X(1) X(2) X(4) X(3)],Y(:)+b,[Z(1) Z(2) Z(4) Z(3)],C(i))
    
    [Y,Z] = meshgrid(yn(i)-b/2:b:yn(i)+b/2,zn(i)-c/2:c:zn(i)+c/2);
    X = ones(size(X))*(xn(i)-a/2);
    patch(X(:)  ,[Y(1) Y(2) Y(4) Y(3)],[Z(1) Z(2) Z(4) Z(3)],C(i))
    patch(X(:)+a,[Y(1) Y(2) Y(4) Y(3)],[Z(1) Z(2) Z(4) Z(3)],C(i))
    
end
alpha(1)
ax4 = gca;
panduan=find(DATA1(:,4) < limitC);

X = DATA1(:,1);X(panduan) = [];
Y = DATA1(:,2);Y(panduan) = [];
Z = -DATA1(:,3);Z(panduan) = [];
M = DATA1(:,4);M(panduan) = [];
U = DATA1(:,6);U(panduan) = [];
V = DATA1(:,7);V(panduan) = [];
W = DATA1(:,8);W(panduan) = [];

Color =[ones(length(X),1) zeros(length(X),1) zeros(length(X),1)];

hold on
for i = 1:length(X)
    cm = floor((M(i)-Mmin)/((Mmax-Mmin)/256));
    quiver3(X(i),Y(i),Z(i),U(i),V(i),W(i),100,'color',cmap(cm,:),'MaxHeadSize',5);
end

cl2=colorbar;
cl2.Label.String = 'Magnetization (A/m)';
set(cl2,'fontsize',10,'fontname','Times');
ax1 = gca;
axpos = ax1.Position;
cl2.Position(1) = 1.03*cl2.Position(1);% x offset
cl2.Position(2) = 2.0*cl2.Position(2);% y offset
cl2.Position(3) = 0.5*cl2.Position(3);% width
cl2.Position(4) = 0.35*cl2.Position(4);% length
ax1.Position = axpos;

colormap(ax1,jet)
caxis([Mmin,Mmax])
grid on
axis equal
axis([xmin xmax ymin ymax zmin zmax])
xlabel('Easting(m)','rotation',0);
ylabel('Northing(m)','rotation',0);
zlabel('Depth(m)','rotation',90);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1,'LineWidth',0.5);
set(gca,'fontsize',10,'fontname','Times');
set(gca,'ZDir','reverse')
view([Azimuth,Elevation+10])

h_a = gca;
[thx,thy,thz] = axislabel_rotation_angle(h_a);
set(get(h_a,'xlabel'),'rotation',thx(1));
set(get(h_a,'ylabel'),'rotation',thy(1));
set(get(h_a,'zlabel'),'rotation',thz(1));