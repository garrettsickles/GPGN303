function polymod(action);
%
% function polymod(inpf)
%
% function to model a profile of gravity data by a single
% polygonal shaped body.
%
% The program allows user to pick the coordinates of the
% polygonal vertices by clicking the left mouse button
% in the model plot section of the plot window.
%
% A separate control button also allows the user to select
% pre-set density contrast values.
%
% An option is also available for inputing a input file
% that contains a user defined model.
%
% Upon starting the program, the user must supply a data
% file that contains the gravity profile to be modelled.
%
% The following are the format for each of the two input files:
%
% Data file: contains the gravity data.
%-------------------------------------------------------
% ndat                       <--- number of data points
% x(1), g(1)                 <--- x-location and value of the gravity
% .
% .
% x(ndat), g(ndat)           <--- x-location and valu eof the gravity
%--------------------------------------------------------
%
% Model file: This file
%--------------------------------------------------------
% nc, rho                    <--- # of vertices, density
% xc(1), zc(1)               <--- coord of vertices, must be order
% .                               clockwise!
% .
% xc(nc), zc(nc)
%--------------------------------------------------------
%
% Author: Yaoguo Li
% Date:   18/10/99
% Written for the GPGN303
%
% Expanded: 08/09/2011 by Richard Krahenbuhl
%
%
% open input file and read in the data
global nc xc  zc  rho xo  zo  gdat  gpre stn ndat ...
   xmin xmax;
%
%
%
if nargin <1,
   action='init';
end;

if strcmp(action,'init'),
   inpf=input(' data file: ');
   
   %=============================================
   oldFigNumber=watchon;
   %scrsz=get(0,'ScreenSize');
   %'Position',[scrsz(3)/8 scrsz(4)/4 2.5*scrsz(3)/4 scrsz(4)/2],...
      fig = figure( ...
   	'Visible','off', ...
      'NumberTitle','off', ...
      'Name','Gravity Modelling by a Polygon');
   colordef(fig,'white');
   axes( ...
      'Units','normalized', ...
      'Position',[0.04 0.17 0.70 0.70]);
 
	%====================================   
	% The CONSOLE frame
	uicontrol( ...
		'Style','frame', ...
   	'Units','normalized', ...
   	'Position',[0.85 0.17 0.14 0.70], ...
      'BackgroundColor',[0.45 0.45 0.45]);
   %=============================================   
   % The Modify button.
	uicontrol( ...
		'Style','pushbutton', ...
   	'Units','normalized', ...
   	'Position',[0.86 0.69 0.12 0.06], ...
   	'String','New Polygon', ...
      'Callback','polymod modify');
   %=============================================
   % The DENSITY command popup button
	% Generic label information
	uicontrol( ...
	 	'Style','text', ...
    	'Units','normalized', ...
 	 	'Position',[0.86 0.57 0.12 0.06], ...
 	 	'BackgroundColor',[0.8 0.8 0.8], ...
	 	'HorizontalAlignment','center', ...
	 	'String','rho (g/cm^3)');
	% Generic popup button information
	dns=uicontrol( ...
    	'Style','popup', ...
    	'Units','normalized', ...
    	'Position',[0.86 0.53 0.12 0.06], ...
    	'String','-3.0|-2.8|-2.6|-2.4|-2.2|-2.0|-1.8|-1.6|-1.4|-1.2|-1.0|-0.8|-0.6|-0.4|-0.2|0.0|0.2|0.4|0.6|0.8|1.0|1.2|1.4|1.6|1.8|2.0|2.2|2.4|2.6|2.8|3.0', ...
       'Callback','polymod density');
       %'String','0.2|0.4|0.6|0.8|1.0|1.2|1.4|1.6|1.8|2.0', ...

   %=============================================
   % The model command popup button
   % Generic label information
   uicontrol( ...
		'Style','pushbutton', ...
   	'Units','normalized', ...
   	'Position',[0.86 0.41 0.12 0.06], ...
   	'String','User Model', ...
      'Callback','polymod model');

	%=============================================
   %The save button.
	uicontrol( ...
		'Style','pushbutton', ...
   	'Units','normalized', ...
   	'Position',[0.86 0.29 0.12 0.06], ...
   	'String','Save', ...
      'Callback','polymod save');
	%=============================================   
	% The CLOSE button.
	uicontrol( ...
		'Style','pushbutton', ...
   	'Units','normalized', ...
   	'Position',[0.86 0.19 0.12 0.06], ...
   	'String','Close', ...
      'Callback','close(gcf)');
	%=============================================
	% Uncover the figure
    hndlList=[dns];
	watchoff(oldFigNumber);
	set(fig,'Visible','on','UserData',hndlList);

   fi=fopen(inpf,'rt');

   ndat=fscanf(fi,'%d',1);
   for ii=1:ndat,
      xo(ii)=fscanf(fi,'%f',1);
      gdat(ii)=fscanf(fi,'%f',1);
      stn(ii)=fscanf(fi,'%f',1);
   end;
   %gdat = gdat/100;
   %stn = stn/100;
   fclose(fi);

   xmin=min(min(xo));
   xmax=max(max(xo));

   % set up the x-location, set z-location to zero
   zo=0.0*xo;

   % set up the data array
   gpre=gdat*0.0;
   rho=0.2;
   rho = -3.0;

   % an 'initial mode'
   nc=4;
   xc(1)=-100.0;
   zc(1)=50.0;
   xc(2)=100.0;
   zc(2)=50.0;
   xc(3)=100.0;
   zc(3)=500.0;
   xc(4)=-100.0;
   zc(4)=500.0;
   polymod('calc');

elseif strcmp(action,'calc'),
   % call the function gpoly which calculates the
   % gravity anomaly along a profile due to a single
   % polygon.
   gpre=gpoly(xo, zo, xc, zc, nc, rho);
   %
   % plot the polygon
   modplot(xc,zc,nc,xo,zo,gdat,gpre,stn,ndat);
      
elseif strcmp(action,'modify'),
   nc=0;
   idone=0;
   while idone==0,
      nc=nc+1;
      [xc(nc),zc(nc),button]=ginput(1);
      if button==3; idone=1; end;
   end;
   polymod('calc');
elseif strcmp(action,'density'),
   hndlList=get(gcf,'UserData');
   tmp = (get(hndlList(1),'Value'));  % added for expansion to negatives
   range = -3.0:0.2:3.0;
   rho=range(tmp);
   %rho=0.2*(get(hndlList(1),'Value'));
   rho;
   polymod('calc');
elseif strcmp(action,'model'),
   usermod=input(' Input user model file: ');
   fin=fopen(usermod,'rt');
   nc=fscanf(fin,'%d',1);
   rho=fscanf(fin,'%f',1);
   for i=1:nc,
      xc(i)=fscanf(fin,'%f',1);
      zc(i)=fscanf(fin,'%f',1);
   end;
   fclose(fin);
   polymod('calc');
elseif strcmp(action,'save'),
   fo=fopen('polymod.mod','wt');
   fprintf(fo,'%6d',nc);
   fprintf(fo,'%9.2f\n',rho);
   for i=1:nc,
      fprintf(fo,'%8.2f',xc(i));
      fprintf(fo,'%8.2f\n',zc(i));
   end;
   fclose(fo);
   figure(2);
   
   modplot(xc,zc,nc,xo,zo,gdat,gpre,stn,ndat);
   %print('-dps', 'polymod.ps');

end;
%
% end of the function polymod
function modplot(xc,zc,nc,xo,zo,gdat,gpre,stn,ndat)
   xp=xc;
   zp=zc;
   xp(nc+1)=xc(1);
   zp(nc+1)=zc(1);
   
   xmin=min(min(xo));
   xmax=max(max(xo));

   %subplot(3,1,3),
   subplot('position',[0.1 0.15 0.72 0.2])
   plot(xp(1:nc+1),zp(1:nc+1));
   axis([xmin xmax 0 2*max(max(zp))]);
   %axis([xmin xmax 0 1000]);

   set(gca, 'Ydir', 'reverse');
   xlabel(['x (m)']);
   ylabel(['z (m)']);

   %
   % plot the gravity profile
   %
   %subplot(2,1,1)
   subplot('position',[0.1 0.5 0.72 0.4])
   errorbar(xo,gdat,stn); hold on
   plot(xo,gpre,'g--');
   gmax=1.1*max(max(gpre),max(gdat));
   gmin=1.1*min(min(gpre),min(gdat));
   axis([xmin xmax gmin gmax]);
   xlabel(['x (m)']);
   ylabel(['g (mGal)']);
   hold off
   %end of the function modplot
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
