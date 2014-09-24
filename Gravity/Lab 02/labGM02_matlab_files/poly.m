function poly(inpf);
%
% function poly(inpf)
%
% function to calculate the gravity field due to one
% or more polygonal 2D causative bodies.
%
% Require Input file "poly.inp" in the following format:
%-------------------------------------------------------
% title            <---- Title/Name of the model
% plot.ps          <---- Output file name for the plot
% xmin xmax delx   <---- end and interval of profile
% npoly            <---- # of polygons
%
% nc1 rho1         <---- # of vertices and density
% xc(1) zc(1)      <---- vertices ordered clockwise (IMPORTANT)
% xc(2) zc(2)
% ...
% xc(nc) zc(nc)
%
% nc2 rho2
% xc(1) zc(1)
% xc(2) zc(2)
% ...
% xc(nc2) zc(nc2)
%
% ....
%--------------------------------------------------------
%
% Author: Yaoguo Li
% Date:   27/09/99
% Written for the GPGN303 lab
%

% open input file and read in the profile position
if nargin <1,
   inpf=uigetfile('*.inp','Input parameter file');
end;
fi=fopen(inpf,'rt');

model=fgetl(fi);
outf=fgetl(fi);
xmin=fscanf(fi,'%f',1);
xmax=fscanf(fi,'%f',1);
delx=fscanf(fi,'%f',1);
depthMax=0;

% set up the x-location, set z-location to zero
xo=xmin:delx:xmax;
zo=0.0*xo;

% set up the data array
gdat=zo.*0.0;
gtmp=gdat;

close;
figure(1);
subplot(3,1,3),

% begin computing the response from each polygon and
% sum up to produce the total gravity anomaly
npoly=fscanf(fi,'%d',1);
for ip=1:1:npoly,
    nc=fscanf(fi,'%d',1);
    rho=fscanf(fi,'%f',1);
    for ic=1:1:nc,
        xc(ic)=fscanf(fi,'%f',1);
        zc(ic)=fscanf(fi,'%f',1);
    end;
    
    % call the function gpoly which calculates the
    % gravity anomaly along a profile due to a single
    % polygon.
    gtmp=gpoly(xo, zo, xc, zc, nc, rho);
    gdat=gdat+gtmp;
    
    % plot the polygon
    xp=xc;
    zp=zc;
    xp(nc+1)=xc(1);
    zp(nc+1)=zc(1);
    plot(xp,zp);
    hold on;
    if max(zc)>depthMax
        depthMax=max(zc);
    end;

end;
axis([xmin xmax 0 depthMax*1.5]);
set(gca, 'Ydir', 'reverse');
xlabel(['x (m)']);
ylabel(['z (m)']);

%
st=fclose(fi);
%
% plot the gravity profile
%
subplot(2,1,1), plot(xo,gdat);
title(model);
xlabel(['x (m)']);
ylabel(['g (mGal)']);
hold off;

%
% end of the program
%

function g=gpoly(xo,zo,xc,zc,nc,rho);
%
% function g=gpoly(xo,zo,xc,zc,nc,rho)
%
% Calculate the gravity field of a 2D prism
% with polygonal cross-section.
%
% input:
%
%    xo, zo: coordinates of the observation point (arrays)
%    xc, zc: coordinates of the vertices of the polygon (arrays)
%            ordered clockwise
%    nc:     number of vertices
%    rho:    density
%
% output:
%
%    g:        calculated gravity data
%
gcons=6.672E-03;

sum=xo.*0.0;

for ic=1:1:nc,

   if ic==nc;
      ic2=1;
   else
      ic2=ic+1;
   end;

   xcb=xc(ic);
   xce=xc(ic2);

   dx=abs(xc(ic2)-xc(ic));
   dz=abs(zc(ic2)-zc(ic));
   if dz > 1.0E-6, 
      zcb=zc(ic);
      zce=zc(ic2);
   else
      zcb=zc(ic)+0.00001*dx;
      zce=zc(ic2);
   end;

   x1=xcb-xo;
   x2=xce-xo;

   z1=zcb-zo;
   z2=zce-zo;

   rt1=x1.*x1 + z1.*z1;
   rt2=x2.*x2 + z2.*z2;

   bot=z2-z1;

   alpha=(x2-x1)./bot;
   beta=(x1.*z2 - x2.*z1)./bot;

   factor=beta./(1.0+alpha.*alpha);

   term1=log(rt2./rt1)./2;
   term2=atan2(z2,x2)-atan2(z1, x1);

   sum=sum + factor.*(term1-alpha.*term2);

end;
     
g=2.0*rho*gcons*sum;

%
% end of function
%