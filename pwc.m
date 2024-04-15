% This is a Matlab code (pwc.m) for calculating partial wavelet coherency. 
 
function varargout=pwc(X,varargin)
%  Partial Wavelet coherency
% Creates a figure of partial wavelet coherency
% USAGE: [Rsq,period,scale,coi,sig95]=pwc(X,[,settings])
%
% Input: X: a matrix of multiple variables equally distributed in space 
%             or time. The first column corresponds to the dependent variable,
%          the second column corresponds to the independent variable, 
%         and the third and subsequent columns are excluding variables.
% 
% Settings: Pad: pad the time series with zeros? 
% .         Dj: Octaves per scale (default: '1/12')
% .         S0: Minimum scale
% .         J1: Total number of scales
% .         Mother: Mother wavelet (default 'morlet')
% .         MaxScale: An easier way of specifying J1
% .         MakeFigure: Make a figure or simply return the output.
% .         BlackandWhite: Create black and white figures
% .         AR1: the ar1 coefficients of the series 
% .              (default='auto' using a naive ar1 estimator. See ar1nv.m)
% .         MonteCarloCount: Number of surrogate data sets in the 
%             significance calculation. (default=300)
 
% Settings can also be specified using abbreviations. e.g. ms=MaxScale.
% For detailed help on some parameters type help wavelet.
% Example:
%    t=[1:200]';
%    pwc([sin(t),sin(t.*cos(t*.01)),cos(t.*sin(t*.01))]) % for one
%    excluding variable
%    pwc([sin(t),sin(t.*cos(t*.01)),cos(t),cos(t.*sin(t*.01))])  % for two
%    excluding variables
% Please acknowledge the use of this software package in any publications,
% by including text such as:
 
%   "The software for the partial wavelet coherency was provided by Wei Hu
%    and is available from https://doi.org/10.6084/m9.figshare.13031123."
%   and cite the paper:
% "Hu, W., and Si,B (2021), Technical Note: Improved partial wavelet coherency 
% for understanding scale-specific and localized bivariate relationships in geosciences,
% Hydrol. Earth Syst. Sci., 25, 321-331.
% 
%  
%  (C) W. Hu 2020
%
% -----------------------------------------------------------
%
%   Copyright (C) 2020, W. Hu 2020
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.
%
%    Wavelet software was provided by C. Torrence and G. Compo,
%       and is available at URL: http://paos.colorado.edu/research/wavelets/.
%
%    Crosswavelet and wavelet coherence software were provided by
%      Aslak Grinsted and is available at URL:
%    http://www.glaciology.net/wavelet-coherence
%     
%        Multiple wavelet coherence software are provided by
%           Hu and Si (2016) and is available at URL:
%      (https://www.hydrol-earth-syst-sci.net/20/3183/2016/hess-20-3183-2016-
%     supplement.pdf)
% The updated version for faster calculation is vailable at URL:
%         (https://doi.org/10.6084/m9.figshare.13031123)
 
%        Partial wavelet coherency software are provided by
%           Hu and Si (2020) and is available at URL:
%         (https://doi.org/10.6084/m9.figshare.13031123)
%   We acknowledge Aslak Grinsted for his wavelet coherency code (wtc.m) on
% which this code builds. 
% 
%---------parse function arguments------------------------------------------
 
[row,col]=size(X);
[y1,dt1]=formatts(X(:,1));
mm=y1(1,1);
nn=y1(end,1);
 
[y2,dt2]=formatts(X(:,2));
mm2=y2(1,1);
nn2=y2(end,1);
 
for i=3:col;
[x,dtx]=formatts(X(:,i));
 
if (dt1~=dtx |dt2~=dtx)
    error('timestep must be equal between time series');
end
 
mm1=x(1,1);
nn1=x(end,1);
 
 
mm=max([mm,mm1,mm2]);
nn=min([nn,nn1,nn2]);
 
x1(:,(i-2))=x(:,1);
x2(:,(i-2))=x(:,2);
 
end
 
t=(mm:dt1:nn)';
 
 
%common time period
if length(t)<4
    error('The three time series must overlap.');
end 
  
n=length(t);
 
%----------default arguments for the wavelet transform-----------
Args=struct('Pad',1,...      % pad the time series with zeroes (recommended)
            'Dj',1/12, ...    % this will do 12 sub-octaves per octave
            'S0',2*dt1,...    % this says start at a scale of 2 years
            'J1',[],...
            'Mother','Morlet', ...
            'MaxScale',[],...   %a more simple way to specify J1
            'MakeFigure',(nargout==0),...
            'MonteCarloCount',300,...
            'BlackandWhite',0,...
            'AR1','auto',...
            'ArrowDensity',[30 30],...
            'ArrowSize',1,...
            'ArrowHeadSize',1);
 
Args=parseArgs(varargin,Args,{'BlackandWhite'});
 
if isempty(Args.J1)
    if isempty(Args.MaxScale)
        Args.MaxScale=(n*.17)*2*dt1; %auto maxscale;
    end
    Args.J1=round(log2(Args.MaxScale/Args.S0)/Args.Dj);
end
 
ad=mean(Args.ArrowDensity);
Args.ArrowSize=Args.ArrowSize*30*.03/ad;
%Args.ArrowHeadSize=Args.ArrowHeadSize*Args.ArrowSize*220;
Args.ArrowHeadSize=Args.ArrowHeadSize*120/ad;
 
if ~strcmpi(Args.Mother,'morlet')
    warning('MWC:InappropriateSmoothingOperator','Smoothing operator is designed for morlet wavelet.');
end
 
if strcmpi(Args.AR1,'auto')
      for i=1:col
arc(i)= ar1nv(X(:,i));
       end
Args.AR1=arc
    if any(isnan(Args.AR1))
        error('Automatic AR1 estimation failed. Specify it manually (use arcov or arburg).');
    end
end
 
%-----------:::::::::::::--------- ANALYZE ----------::::::::::::------------
 
%Calculate and smooth wavelet spectrum Y1, Y2, and X
  
[Y1,period,scale,coiy1] = wavelet(y1(:,2),dt1,Args.Pad,Args.Dj,Args.S0,Args.J1,Args.Mother);
sinv=1./(scale');
smY1=smoothwavelet(sinv(:,ones(1,n)).*(abs(Y1).^2),dt1,period,Args.Dj,scale);
dte=dt1*.01; 
idx=find((y1(:,1)>=(t(1)-dte))&(y1(:,1)<=(t(end)+dte)));
Y1=Y1(:,idx);
smY1=smY1(:,idx);
coiy1=coiy1(idx);
coi1=coiy1;
 
 
[Y2,period,scale,coiy2] = wavelet(y2(:,2),dt1,Args.Pad,Args.Dj,Args.S0,Args.J1,Args.Mother);
sinv=1./(scale');
smY2=smoothwavelet(sinv(:,ones(1,n)).*(abs(Y2).^2),dt1,period,Args.Dj,scale);
dte=dt2*.01; 
idx=find((y2(:,1)>=(t(1)-dte))&(y2(:,1)<=(t(end)+dte)));
Y2=Y2(:,idx);
smY2=smY2(:,idx);
coiy2=coiy2(idx);
coi2=coiy2;
 
 
 
for  i=3:col
 [XS,period,scale,coix] = wavelet(x2(:,(i-2)),dt1,Args.Pad,Args.Dj,Args.S0,Args.J1,Args.Mother);
 
idx=find((x1(:,(i-2))>=(t(1))-dte)&(x1(:,(i-2))<=(t(end)+dte)));
XS=XS(:,idx);
coix=coix(idx);
 
XS1(:,:,(i-2))=XS;
coi0=min(coi1,coi2);
coi=min(coi0,coix);
end
 
% -------- Calculate Cross Wavelet Spectra----------------------------
 
% -- between dependent variable (or independent variables) and excluding
% factors------
 
for i=1:(col-2)
Wy1x=Y1.*conj(XS1(:,:,i));
sWy1x=smoothwavelet(sinv(:,ones(1,n)).*Wy1x,dt1,period,Args.Dj,scale);
sWy1x1(:,:,i)=sWy1x;
 
Wy2x=Y2.*conj(XS1(:,:,i));
sWy2x=smoothwavelet(sinv(:,ones(1,n)).*Wy2x,dt1,period,Args.Dj,scale);
sWy2x1(:,:,i)=sWy2x;
end
 
% ---- between dependent variable and independent variables------
Wy1y2=Y1.*conj(Y2);
sWy1y2=smoothwavelet(sinv(:,ones(1,n)).*Wy1y2,dt1,period,Args.Dj,scale);
 
% ----between excluding variables and excluding variables------
for i=1:(col-2);
for j=1:(col-2);
Wxx=XS1(:,:,i).*conj(XS1(:,:,j));
sWxx=smoothwavelet(sinv(:,ones(1,n)).*Wxx,dt1,period,Args.Dj,scale);
sWxx1(:,:,i,j)=sWxx;
end
end
 
% --------------- Partial wavelet coherence ---------------------------------
% calculate the partial wavelet coherence 
m=length(scale);
Rsq=zeros(m,n);
Wy1y2x=zeros(m,n);
nofe=col-2; %number of excluding factors
 
if nofe==1;    
% ------Partial wavelet coherence for one excluding variable ------------------
%%%%%%%%% |1-R^2 y,x·Z (s,t)|^2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:m;
    for jj=1:n;
        sWxx1_i(ii,jj)=inv(sWxx1(ii,jj)); % inversed sWxx1
    end
end
 
sWy1x1_sWxx1=bsxfun(@times,sWy1x1,sWxx1_i); % sWy1x1*sWxx1_i
sWy2x1_c=conj(sWy2x1); %conjugate of smoothed cross-wavelet power spectra
% between independent variable and excluding variables
sWy1x1_sWxx1_sWy2x1=bsxfun(@times,sWy1x1_sWxx1,sWy2x1_c); 
% sWy1x1*sWxx1_i*sWy2x1_c
R2y1y2x=sWy1x1_sWxx1_sWy2x1./sWy1y2; %(Ry,x.Z)^2
abs_one_minus_R2y1y2x=(abs(1-R2y1y2x)).^2; %squared absolute of 1-(Ry,x.Z)^2
%%%%%%%%%%%%%%%% R^2 y,x (s,t)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R2y1y2=(sWy1y2.*conj(sWy1y2))./(smY1.*smY2); %(Ry,x)^2
%%%%%%%%%%%%%%%%%%% 1-R^2 y,Z(s,t)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sWy1x1_c=conj(sWy1x1);
sWy1x1_sWxx1_sWy1x1=bsxfun(@times,sWy1x1_sWxx1,sWy1x1_c);
R2y1x=sWy1x1_sWxx1_sWy1x1./smY1; %(Ry,Z)^2
one_minus_R2y1x=(1-R2y1x);
 
%%%%%%%%%%%%%%%%%%%%%  1-R^2 x,Z(s,t) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sWy2x1_sWxx1=bsxfun(@times,sWy2x1,sWxx1_i);
sWy2x1_sWxx1_sWy2x1=bsxfun(@times,sWy2x1_sWxx1,sWy2x1_c);
R2y2x=sWy2x1_sWxx1_sWy2x1./smY2; %(Rx,Z)^2
one_minus_R2y2x=(1-R2y2x);
%%%%%%%%%%%%%%%%%%  R^2 y,x·Z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rsq_num=roundn(bsxfun(@times,abs_one_minus_R2y1y2x,R2y1y2),-16); 
%numerator part of the equation for squared pwc, if the value is less than 
%10^(-16), 0 is assigned
Rsq_den=real(bsxfun(@times,one_minus_R2y1x,one_minus_R2y2x));  
%denominator part of the equation for squared pwc
Rsq= bsxfun(@rdivide,Rsq_num,Rsq_den); %squared partial wavelet coherence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
Wy1y2x=1-R2y1y2x;
 
else
 
% -----Partial wavelet coherence for two or more excluding variables ----------
 
%%%%%%%% |1-R^2 y,x·Z (s,t)|^2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sWy1x1_p=permute(sWy1x1,[3,1,2]); %re-structure sWy1x1
sWy1x1_rp=reshape(sWy1x1_p,[1,nofe,m,n]);%change 3D to 4D matrix
sWxx1_p=permute(sWxx1,[4,3,1,2]); % re-structure sWxx1 (smoothed auto-or cross-
% wavelet power spectra for excluding factors)
for jj=1:n;
    for ii=1:m;
        sWxx1_ip(:,:,ii,jj)=inv(sWxx1_p(:,:,ii,jj)); % inversed sWxx1
    end
end
sWy1x1_sWxx1=bsxfun(@times,sWy1x1_rp,sWxx1_ip); % sWy1x1_rp*sWxx1_ip
sWy1x1_sWxx1_s=sum(sWy1x1_sWxx1,2);
sWy1x1_sWxx1_ps=permute(sWy1x1_sWxx1_s,[2,1,3,4]);  % re-structure 
sWy1x1_sWxx1_sps=squeeze(sWy1x1_sWxx1_ps); % 4D to 3D
 
sWy2x1_c=conj(sWy2x1); %conjugate of smoothed cross-wavelet power spectra
% between independent variable and excluding variables
sWy2x1_pc=permute(sWy2x1_c,[3,1,2]); % re-structure
 
sWy1x1_sWxx1_sWy2x1=bsxfun(@times,sWy1x1_sWxx1_sps,sWy2x1_pc); 
% sWy1x1_rp*sWxx1_ip*sWy2x1_pc
sWy1x1_sWxx1_sWy2x1_s=sum(sWy1x1_sWxx1_sWy2x1,1);
sWy1x1_sWxx1_sWy2x1_ss=squeeze(sWy1x1_sWxx1_sWy2x1_s);
 
R2y1y2x=sWy1x1_sWxx1_sWy2x1_ss./sWy1y2; %(Ry,x.Z)^2
abs_one_minus_R2y1y2x=(abs(1-R2y1y2x)).^2; %squared absolute of 1-(Ry,x.Z)^2
%%%%%%%%%%%%%%%% R^2 y,x (s,t)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R2y1y2=(sWy1y2.*conj(sWy1y2))./(smY1.*smY2); %(Ry,x)^2
%%%%%%%%%%%%%%%%%%% 1-R^2 y,Z(s,t)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sWy1x1_c=conj(sWy1x1);
sWy1x1_pc=permute(sWy1x1_c,[3,1,2]);
sWy1x1_sWxx1_sWy1x1=bsxfun(@times,sWy1x1_sWxx1_sps,sWy1x1_pc);
sWy1x1_sWxx1_sWy1x1_s=sum(sWy1x1_sWxx1_sWy1x1,1);
sWy1x1_sWxx1_sWy1x1_ss=squeeze(sWy1x1_sWxx1_sWy1x1_s);
R2y1x=sWy1x1_sWxx1_sWy1x1_ss./smY1; %(Ry,Z)^2
one_minus_R2y1x=(1-R2y1x);
 
%%%%%%%%%%%%%%%%%%%%%  1-R^2 x,Z(s,t) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sWy2x1_p=permute(sWy2x1,[3,1,2]);
sWy2x1_rp=reshape(sWy2x1_p,[1,nofe,m,n]);
sWy2x1_sWxx1=bsxfun(@times,sWy2x1_rp,sWxx1_ip);
sWy2x1_sWxx1_s=sum(sWy2x1_sWxx1,2);
sWy2x1_sWxx1_ps=permute(sWy2x1_sWxx1_s,[2,1,3,4]);
sWy2x1_sWxx1_sps=squeeze(sWy2x1_sWxx1_ps);
sWy2x1_sWxx1_sWy2x1=bsxfun(@times,sWy2x1_sWxx1_sps,sWy2x1_pc);
sWy2x1_sWxx1_sWy2x1_s=sum(sWy2x1_sWxx1_sWy2x1,1);
sWy2x1_sWxx1_sWy2x1_ss=squeeze(sWy2x1_sWxx1_sWy2x1_s);
R2y2x=sWy2x1_sWxx1_sWy2x1_ss./smY2; %(Rx,Z)^2
one_minus_R2y2x=(1-R2y2x);
%%%%%%%%%%%%%%%%%%  R^2 y,x·Z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rsq_num=roundn(bsxfun(@times,abs_one_minus_R2y1y2x,R2y1y2),-16); 
%numerator part of the equation for squared pwc, if the value is less than 
% 10^(-16), 0 is assigned
Rsq_den=real(bsxfun(@times,one_minus_R2y1x,one_minus_R2y2x));  
%denominator part of the equation for squared pwc
Rsq= bsxfun(@rdivide,Rsq_num,Rsq_den); %squared partial wavelet coherence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
Wy1y2x=1-R2y1y2x;
end
 
 
% --------------- make figure--------------------------------------------
if (nargout>0)||(Args.MakeFigure)
    pwcsig=pwcsignif(Args.MonteCarloCount,Args.AR1,dt1,length(t)*2,Args.Pad,Args.Dj,Args.S0,Args.J1,Args.Mother,.6);
    pwcsig=(pwcsig(:,2))*(ones(1,n));
    pwcsig=Rsq./pwcsig;
end
 
if Args.MakeFigure
   
    Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
    
    if Args.BlackandWhite
        levels = [0 0.5 0.7 0.8 0.9 1];
        [cout,H]=safecontourf(t,log2(period),Rsq,levels);
 
        colorbarf(cout,H)
        cmap=[0 1;.5 .9;.8 .8;.9 .6;1 .5];
        cmap=interp1(cmap(:,1),cmap(:,2),(0:.1:1)');
        cmap=cmap(:,[1 1 1]);
        colormap(cmap)
        set(gca,'YLim',log2([min(period),max(period)]), ...
            'YDir','reverse', 'layer','top', ...
            'YTick',log2(Yticks(:)), ...
            'YTickLabel',num2str(Yticks'), ...
            'layer','top');
        ylabel('Period');
        hold on
 
        %phase plot
       aWxy=angle(Wy1y2)+angle(Wy1y2x);
               aaa=aWxy;
        aaa(Rsq<.5)=NaN; %remove phase indication where Rsq is low
           aaa(isnan(Rsq))=NaN; %remove phase indication where Rsq is NaN
           %[xx,yy]=meshgrid(t(1:5:end),log2(period));
 
        phs_dt=round(length(t)/Args.ArrowDensity(1)); tidx=max(floor(phs_dt/2),1):phs_dt:length(t);
        phs_dp=round(length(period)/Args.ArrowDensity(2)); pidx=max(floor(phs_dp/2),1):phs_dp:length(period);
        phaseplot(t(tidx),log2(period(pidx)),aaa(pidx,tidx),Args.ArrowSize,Args.ArrowHeadSize);
 
        if ~all(isnan(pwcsig))
            [c,h] = contour(t,log2(period),pwcsig,[1 1],'k');%#ok
            set(h,'linewidth',2);
        end
        %suptitle([sTitle ' coherence']);
        %plot(t,log2(coi),'k','linewidth',2)
                tt=[t([1 1])-dt1*.5;t;t([end end])+dt1*.5];
        %hcoi=fill(tt,log2([period([end 1]) coi period([1 end])]));
        %hatching- modified by Ng and Kwok
        hcoi=fill(tt,log2([period([end 1]) coi period([1 end])]),'w');
  
        hatch(hcoi,45,[0 0 0]);
        hatch(hcoi,135,[0 0 0]);
        set(hcoi,'alphadatamapping','direct','facealpha',.5);
        plot(t,log2(coi),'color','black','linewidth',1.5);
        hold off
    else
        H=imagesc(t,log2(period),Rsq);%#ok
        %[c,H]=safecontourf(t,log2(period),Rsq,[0:.05:1]);
        %set(H,'linestyle','none')
        
        set(gca,'clim',[0 1]);
        
        HCB=safecolorbar;%#ok
        
        set(gca,'YLim',log2([min(period),max(period)]), ...
            'YDir','reverse', 'layer','top', ...
            'YTick',log2(Yticks(:)), ...
            'YTickLabel',num2str(Yticks'), ...
            'layer','top');
        ylabel('Period');
        hold on
        
  %phase plot
       aWxy=angle(Wy1y2)+angle(Wy1y2x);
               aaa=aWxy;
        aaa(Rsq<.5)=NaN; %remove phase indication where Rsq is low
        aaa(isnan(Rsq))=NaN; %remove phase indication where Rsq is NaN
        %[xx,yy]=meshgrid(t(1:5:end),log2(period));
 
        phs_dt=round(length(t)/Args.ArrowDensity(1)); tidx=max(floor(phs_dt/2),1):phs_dt:length(t);
        phs_dp=round(length(period)/Args.ArrowDensity(2)); pidx=max(floor(phs_dp/2),1):phs_dp:length(period);
        phaseplot(t(tidx),log2(period(pidx)),aaa(pidx,tidx),Args.ArrowSize,Args.ArrowHeadSize);
        
        if ~all(isnan(pwcsig))
            [c,h] = contour(t,log2(period),pwcsig,[1 1],'k');%#ok
            set(h,'linewidth',2);
        end
        %suptitle([sTitle ' coherence']);
        tt=[t([1 1])-dt1*.5;t;t([end end])+dt1*.5];
        hcoi=fill(tt,log2([period([end 1]) coi period([1 end])]),'w');
        set(hcoi,'alphadatamapping','direct','facealpha',.5);
        hold off
    end
end
%---------------------------------------------------------------%
 
varargout={Rsq,period,scale,coi,pwcsig};
varargout=varargout(1:nargout);
 
function [cout,H]=safecontourf(varargin)
vv=sscanf(version,'%i.');
if (version('-release')<14)|(vv(1)<7)
    [cout,H]=contourf(varargin{:});
else
    %[cout,H]=contourf('v6',varargin{:});
    [cout,H]=contourf(varargin{:});
end
 
function hcb=safecolorbar(varargin)
vv=sscanf(version,'%i.');
 
if (version('-release')<14)|(vv(1)<7)
    hcb=colorbar(varargin{:});
else
    %hcb=colorbar('v6',varargin{:});
    hcb=colorbar(varargin{:});
end 
