function [hfig,ax] = plotInflatedMap_SL_new(data,minv,maxv,cmap,options)
% IN
% - data = vector [2447,1]
% - minv = minimum in colorbar
% - maxv = maximum in colorbar
% - cmap = diff (*Spectral), diff-inv (Spectral), parula, parula-inv, stat (*RdBu), stat-inv (RdBu), pval (*OrRd), grade (*YlOrRd),
% blues (*Blues)
% - options:
% -- nc = number of grade in colorbar (useful to define for latency maps // default = 64
% -- subplot_dist (shrinked or normal) = distance between subplots, shrinked = all SM close in together // default = normal
% -- text (with or without) = LL, LM, RL, RM text in subplots // default = without
% -- colorbar (with or without) 
% -- fig_bg = color of background ('k' is useful for cool videos) // default = none 
% OUT
% - hfig and ax to be able to further manipulate the figure and axes properties 


%% Load dependencies ----------------
Stuff4SM_dir = '\\filebck\dmi\CIRS_BACKUP\Functions\SourceModel\Stuff4SM\'; % ---- adapt to your PC
currdir = pwd;
cd(Stuff4SM_dir)
load([Stuff4SM_dir 'InflatedMapStuff.mat'])

%% Set --------------

% SM greys and empty vertices setup -------------
darkgrey = .75; % .3
lightgrey = .99; % .6
VertexColor(ccvverts,:)=repmat([darkgrey darkgrey darkgrey],length(ccvverts),1);
VertexColor(cvxverts,:)=repmat([lightgrey lightgrey lightgrey],length(cvxverts),1);

% options ------------
% number of color grade - useful to define for latency maps
if ~isfield(options,'nc')
    nc=64;
else
    nc = options.nc;
end
% subplots distance 
if ~isfield(options,'subplot_dist')
    options.subplot_dist='shrinked';
end
% text on subplots
if ~isfield(options,'text')
    options.text='without';
end
% colorbar
if ~isfield(options,'colorbar')
    options.colorbar='without';
end

% colormap -----------
switch cmap
    case 'diff'
        cm = brewermap(nc, '*Spectral');
    case 'diff-inv'
        cm = brewermap(nc, 'Spectral');
        %cm = colormap('parula');
    case 'parula'
        cm = colormap('parula');
    case 'parula-inv'
        cm = flipud(colormap('parula'));
    case 'stat'
        cm = brewermap(nc, '*RdBu');
    case 'stat-inv'
        cm = brewermap(nc, 'RdBu');
    case 'pval'
        cm = brewermap(nc, '*OrRd');
    case 'grade'
        cm = brewermap(nc, '*YlOrRd');
        %         cm = autumn;
    case 'blues'
        cm = brewermap(nc, '*Blues');
end
% cm(end-2:end,:) = repmat([1,0,0],3,1);
% nc=size(colormap,1); ORIGINAL CODE - MODIFIED BY AMS ON 05 OCTOBER 2021
% cm=flipud(colormap(cmap));
%cm=colormap(cmap);%cmap
%cm = spectral;
cgdark=rgb2hsl(cm);
cgdark(:,3)=cgdark(:,3)-.2;%-.2
cgdark(:,2)=cgdark(:,2)-.2;%-.2
cgdark=hsl2rgb(cgdark);
cg=rgb2hsl(cm);
cg(:,2)=cg(:,2)-.4;%-.4
cg=hsl2rgb(cg);
% close gcf;


%% Create data vertices --------------
figcol=transpose((minv:((maxv-minv)/2446):maxv));
[n bin]=hist(figcol,nc-2);
[n2 databin]=histc(data,[-Inf bin Inf]);
% [n2 databin]=histc(data,bin);

dataverts=VertexColor;
for dip=1:length(data)
    verti=find(cDtV==dip);
    if ~isempty(verti) && ~isnan(data(dip))
        ccv=find(dataverts(verti,1)==darkgrey);
        cvx=find(dataverts(verti,1)==lightgrey);
        if ~isempty(ccv);dataverts(verti(ccv),:)=repmat(cgdark(databin(dip),:),length(ccv),1);end
        if ~isempty(cvx);dataverts(verti(cvx),:)=repmat(cg(databin(dip),:),length(cvx),1);end
    end
end


%% Figure --------------
hfig = figure(); 
set(gcf,'position',[100 100 700 550])
if isfield(options,'fig_bg')
    set(gcf,'color',options.fig_bg)
end

% LM -------------------
ax(3)=subplot(2,2,3); 
h=patch('Faces',TriangleVertex(goodtris(1:length(goodtris)/2),:),'Vertices',VertexCoordinate([1:size(VertexCoordinate,1)/2],:),...
    'FaceVertexcData',dataverts([1:size(VertexCoordinate,1)/2],:),'FaceColor','interp','EdgeColor','none');% figure;
axis off;
az=0;el=90;
view(az,el);
set(gca,'ydir','reverse');
set(gca,'xdir','reverse');
if strcmp(options.text, 'with')
text(.1,.15,'LM','units','normalized','FontWeight','bold','FontSize',13,'FontName','Arial')
end
% mask
pos = [.24 .187 .128 .122];
ax(5)=axes('pos',pos);
[a,~,c]=imread('mask_LM_trans_100pt.png');
image(a,'AlphaData',c); set(gca,'Color','none');axis off;
clear a b c

% LL -------------------
ax(1)=subplot(2,2,1); 
h=patch('Faces',TriangleVertex(goodtris(1:length(goodtris)/2),:),'Vertices',VertexCoordinate([1:size(VertexCoordinate,1)/2],:),...
    'FaceVertexcData',dataverts([1:size(VertexCoordinate,1)/2],:),'FaceColor','interp','EdgeColor','none');% figure;
axis off;
el=270;
view(az,el);
set(gca,'ydir','normal');
set(gca,'xdir','normal');
if strcmp(options.text, 'with')
text(.18,.15,'LL','units','normalized','FontWeight','bold','FontSize',13,'FontName','Arial')
end

% RL -------------------
ax(2)=subplot(2,2,2); 
h=patch('Faces',TriangleVertex(goodtris([(length(goodtris)/2+1):length(goodtris)]),:)-(size(VertexCoordinate,1)/2),'Vertices',...
    VertexCoordinate([(size(VertexCoordinate,1)/2+1):size(VertexCoordinate,1)],:),'FaceVertexcData',dataverts([(size(VertexCoordinate,1)/2+1):size(VertexCoordinate,1)],:),'FaceColor','interp','EdgeColor','none');
axis off;
az=0;el=90;
view(az,el);
set(gca,'ydir','reverse');
set(gca,'xdir','reverse');
if strcmp(options.text, 'with')
text(.08,.15,'RL','units','normalized','FontWeight','bold','FontSize',13,'FontName','Arial')
end

% RM -------------------
ax(4)=subplot(2,2,4); 
h=patch('Faces',TriangleVertex(goodtris([(length(goodtris)/2+1):length(goodtris)]),:)-(size(VertexCoordinate,1)/2),'Vertices',...
    VertexCoordinate([(size(VertexCoordinate,1)/2+1):size(VertexCoordinate,1)],:),'FaceVertexcData',dataverts([(size(VertexCoordinate,1)/2+1):size(VertexCoordinate,1)],:),'FaceColor','interp','EdgeColor','none');
axis off;
el=270;
view(az,el);
set(gca,'ydir','normal');
set(gca,'xdir','normal');
% mask
pos = [.664 .183 .129 .124];
ax(6)=axes('pos',pos);
[a,~,c]=imread('mask_RM_trans_100pt.png');
image(a,'AlphaData',c); set(gca,'Color','none');axis off;
clear a b c
if strcmp(options.text, 'with')
text(-.32,-.15,'RM','units','normalized','FontWeight','bold','FontSize',13,'FontName','Arial')
end

%% Set colorbar ------------------

if strcmp(options.colorbar,'with')
hbar=colorbar('Location','SouthOutside');
xlim1 = hbar.Limits(1);
xlim2 = hbar.Limits(2);
set(hbar,'XTick', [xlim1, xlim2], 'xticklabel', {round(minv,4), round(maxv,4)}, 'Position', [0.35 0.5 0.3 0.06],'FontWeight','bold', 'tickdirection', 'out');
end
colormap(cm); %cmap

%% Set subplots distance 
if strcmp(options.subplot_dist,'shrinked')
    % subplot 1 - top left - LL
    pos=get(ax(1),'Position');
    pos(1)=pos(1)+.0675;%0.1975
    pos(2)=pos(2)-0.1938;%0.39;
    set(ax(1),'Position',pos);
    
    % subplot 2 - top right - RL
    pos=get(ax(2),'Position');
    pos(1)=pos(1)-0.0503;%0.52;
    pos(2)=pos(2)-0.1888;%0.395;
    set(ax(2),'Position',pos);
    
    % subplot 3 - bottom left - LM
    pos=get(ax(3),'Position');
    pos(1)=pos(1)+0.09;%0.22;
    set(ax(3),'Position',pos);
    % mask
    pos=get(ax(5),'Position');
    pos(1)=pos(1)+0.09;%0.22;
    set(ax(5),'Position',pos);
    
    % subplot 4 - bottom right - RM
    pos=get(ax(4),'Position');
    pos(1)= pos(1)-.0728;%0.4975;
    set(ax(4),'Position',pos);
    % mask
    pos=get(ax(6),'Position');
    pos(1)= pos(1)-.0728;%0.4975;
    set(ax(6),'Position',pos);
    
end

disp('figure done')
cd(currdir)



%% Functions -----------------

function hsl=rgb2hsl(rgb_in)
%Converts Red-Green-Blue Color value to Hue-Saturation-Luminance Color value
%
%Usage
%       HSL = rgb2hsl(RGB)
%
%   converts RGB, a M [x N] x 3 color matrix with values between 0 and 1
%   into HSL, a M [x N] X 3 color matrix with values between 0 and 1
%
%See also hsl2rgb, rgb2hsv, hsv2rgb

% (C) Vladimir Bychkovsky, June 2008
% written using: 
% - an implementation by Suresh E Joel, April 26,2003
% - Wikipedia: http://en.wikipedia.org/wiki/HSL_and_HSV

rgb=reshape(rgb_in, [], 3);

mx=max(rgb,[],2);%max of the 3 colors
mn=min(rgb,[],2);%min of the 3 colors

L=(mx+mn)/2;%luminance is half of max value + min value
S=zeros(size(L));

% this set of matrix operations can probably be done as an addition...
zeroidx= (mx==mn);
S(zeroidx)=0;

lowlidx=L <= 0.5;
calc=(mx-mn)./(mx+mn);
idx=lowlidx & (~ zeroidx);
S(idx)=calc(idx);

hilidx=L > 0.5;
calc=(mx-mn)./(2-(mx+mn));
idx=hilidx & (~ zeroidx);
S(idx)=calc(idx);

hsv=rgb2hsv(rgb);
H=hsv(:,1);

hsl=[H, S, L];

hsl=round(hsl.*100000)./100000; 
hsl=reshape(hsl, size(rgb_in));

function rgb=hsl2rgb(hsl_in)
%Converts Hue-Saturation-Luminance Color value to Red-Green-Blue Color value
%
%Usage
%       RGB = hsl2rgb(HSL)
%
%   converts HSL, a M [x N] x 3 color matrix with values between 0 and 1
%   into RGB, a M [x N] X 3 color matrix with values between 0 and 1
%
%See also rgb2hsl, rgb2hsv, hsv2rgb

% (C) Vladimir Bychkovsky, June 2008
% written using: 
% - an implementation by Suresh E Joel, April 26,2003
% - Wikipedia: http://en.wikipedia.org/wiki/HSL_and_HSV

hsl=reshape(hsl_in, [], 3);

H=hsl(:,1);
S=hsl(:,2);
L=hsl(:,3);

lowLidx=L < (1/2);
q=(L .* (1+S) ).*lowLidx + (L+S-(L.*S)).*(~lowLidx);
p=2*L - q;
hk=H; % this is already divided by 360

t=zeros([length(H), 3]); % 1=R, 2=B, 3=G
t(:,1)=hk+1/3;
t(:,2)=hk;
t(:,3)=hk-1/3;

underidx=t < 0;
overidx=t > 1;
t=t+underidx - overidx;
    
range1=t < (1/6);
range2=(t >= (1/6) & t < (1/2));
range3=(t >= (1/2) & t < (2/3));
range4= t >= (2/3);

% replicate matricies (one per color) to make the final expression simpler
P=repmat(p, [1,3]);
Q=repmat(q, [1,3]);
rgb_c= (P + ((Q-P).*6.*t)).*range1 + ...
        Q.*range2 + ...
        (P + ((Q-P).*6.*(2/3 - t))).*range3 + ...
        P.*range4;
       
rgb_c=round(rgb_c.*10000)./10000; 
rgb=reshape(rgb_c, size(hsl_in));
