% This is a Matlab file (pwcsignif.m) for calculating significance tests on partial wavelet coherence. 
 

function pwcsig=pwcsignif(mccount,ar1,dt,n,pad,dj,s0,j1,mother,cutoff)
% Partial Wavelet Coherence Significance Calculation (Monte Carlo)
%
% pwcsig=pwcsignif(mccount,ar1,dt,n,pad,dj,s0,j1,mother,cutoff)
%
% mccount: number of time series generations in the monte carlo run 
%(the greater the better)
% ar1: a vector of AR1 coefficients. 
% dt,pad,dj,s0,j1,mother: see wavelet help... 
% n: length of each generated timeseries. (obsolete) 
%
% cutoff: (obsolete)
%
% RETURNED
% pwcsig: the 95% significance level as a function of scale... (scale,sig95level)
% -----------------------------------------------------------
% Please acknowledge the use of this software package in any publications,
% by including text such as:
%
%   "The software for the partial wavelet coherency was provided by Wei Hu
%    and is available from https://doi.org/10.6084/m9.figshare.13031123."
%   and cite the paper:
% "Hu, W., and Si,B (2021), Technical Note: Improved partial wavelet coherency 
% for understanding scale-specific and localized bivariate relationships in geosciences,
% Hydrol. Earth Syst. Sci., 25, 321-331.
%  
%  (C) W. Hu  2020
%
% -----------------------------------------------------------
%  W. Hu  2020
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.
%
%    Wavelet software was provided by C. Torrence and G. Compo,
%       and is available at URL: http://paos.colorado.edu/research/wavelets/.
%
%    Crosswavelet and wavelet coherence software were provided by
%      A. Grinsted and is available at URL:
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

% We acknowledge Aslak Grinsted for his code (wtcsignif.m) on
% which this code builds. 
% 
%---------------------------------------------------------------------
cachedir=fileparts(mfilename('fullpath'));
cachedir=fullfile(cachedir,'.cache');

%we don't need to do the monte carlo if we have a cached
%siglevel for ar1s that are almost the same. (see fig4 in Grinsted et al., 2004)
aa=round(atanh(ar1(:)')*4); %this function increases the sensitivity near 1 & -1
aa=abs(aa)+.5*(aa<0); %only positive numbers are allowed in the checkvalues (because of log)

% do a check that it is not the same as last time... (for optimization purposes)
checkvalues=single([aa dj s0/dt j1 double(mother)]); %n & pad are not important.
%also the resolution is not important.

checkhash=['' mod(sum(log(checkvalues+1)*127),25)+'a' mod(sum(log(checkvalues+1)*54321),25)+'a'];

cachefilename=fullfile(cachedir,['pwcsignif-cached-' checkhash '.mat']);

%the hash is used to distinguish cache files.
try
    last=load(cachefilename);
    if (last.mccount>=mccount) && (isequal(checkvalues,last.checkvalues))
        pwcsig=last.pwcsig;
        return
    end
catch
end

%choose a n so that largest scale have atleast some part outside the coi
ms=s0*(2^(j1*dj))/dt; %maxscale in units of samples
n=ceil(ms*6);

warned=0;
%precalculate stuff that's constant outside the loop
%d1=ar1noise(n,1,ar1(1),1);
d1=rednoise(n,ar1(1),1);
[W1,period,scale,coi] = wavelet(d1,dt,pad,dj,s0,j1,mother);
outsidecoi=zeros(size(W1));
for s=1:length(scale)
    outsidecoi(s,:)=(period(s)<=coi);
end
sinv=1./(scale');
sinv=sinv(:,ones(1,size(W1,2)));

if mccount<1
    mwcsig=scale';
    mwcsig(:,2)=.71; %pretty good 
    return
end

sig95=zeros(size(scale));

maxscale=1;
for s=1:length(scale)
    if any(outsidecoi(s,:)>0)
        maxscale=s;
    else
        sig95(s)=NaN;
        if ~warned
warning('Long wavelengths completely influenced by COI. (suggestion: set n higher, or j1 lower)'); 
         warned=1;
        end
    end
end

%PAR1=1./ar1spectrum(ar1(1),period');
%PAR1=PAR1(:,ones(1,size(W1,2)));
%PAR2=1./ar1spectrum(ar1(2),period');
%PAR2=PAR2(:,ones(1,size(W1,2)));

nbins=1000;
wlc=zeros(length(scale),nbins);

wbh = waitbar(0,['Running Monte Carlo (significance)... (H=' checkhash ')'],'Name','Monte Carlo (PWC)');

for ii=1:mccount
    waitbar(ii/mccount,wbh);

dy1=rednoise(n,ar1(1),1);
[Wdy1,period,scale,coiy] = wavelet(dy1,dt,pad,dj,s0,j1,mother);
sinv=1./(scale');
smdY1=smoothwavelet(sinv(:,ones(1,n)).*(abs(Wdy1).^2),dt,period,dj,scale);


dy2=rednoise(n,ar1(2),1);
[Wdy2,period,scale,coiy] = wavelet(dy2,dt,pad,dj,s0,j1,mother);
sinv=1./(scale');
smdY2=smoothwavelet(sinv(:,ones(1,n)).*(abs(Wdy2).^2),dt,period,dj,scale);

col=size(ar1,2);

for  i=3:col
dx=rednoise(n,ar1(i),1);
 [Wdx,period,scale,coix] = wavelet(dx,dt,pad,dj,s0,j1,mother);
Wdx1(:,:,(i-2))=Wdx;
end

% -------- Calculate Cross Wavelet Spectra----------------------------

% 
% ---- between dependent variable (or independent variables) and excluding factors------

for i=1:(col-2)
Wdy1x=Wdy1.*conj(Wdx1(:,:,i));
sWdy1x=smoothwavelet(sinv(:,ones(1,n)).*Wdy1x,dt,period, dj,scale);
sWdy1x1(:,:,i)=sWdy1x;

Wdy2x=Wdy2.*conj(Wdx1(:,:,i));
sWdy2x=smoothwavelet(sinv(:,ones(1,n)).*Wdy2x,dt,period, dj,scale);
sWdy2x1(:,:,i)=sWdy2x;
end

% ---- between dependent variable and independent variables------
Wdy1y2=Wdy1.*conj(Wdy2);
sWdy1y2=smoothwavelet(sinv(:,ones(1,n)).*Wdy1y2,dt,period,dj,scale);

% ----between excluding variables and excluding variables-
for i=1:(col-2);
for j=1:(col-2);
Wdxx=Wdx1(:,:,i).*conj(Wdx1(:,:,j));
sWdxx=smoothwavelet(sinv(:,ones(1,n)).*Wdxx,dt,period,dj,scale);
sWdxx1(:,:,i,j)=sWdxx;
end
end

% --------------- Partial wavelet coherence  ---------------------------------
% calculate the partial wavelet coherence 
m=length(scale);
Rsq=zeros(m,n);
nofe=col-2; %number of excluding factors

if nofe==1;    
% --------------- Partial wavelet coherence for one excluding variable ---------------------------------
%%%%%%%%%%%%%%%%%%%%% |1-R^2 y,x?Z (s,t)|^2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:m;
    for jj=1:n;
        sWdxx1_i(ii,jj)=inv(sWdxx1(ii,jj)); % inversed sWxx1
    end
end

sWdy1x1_sWdxx1=bsxfun(@times,sWdy1x1,sWdxx1_i); % sWy1x1*sWxx1_i
sWdy2x1_c=conj(sWdy2x1); %conjugate of smoothed cross-wavelet power spectra between independent variable and excluding variables
sWdy1x1_sWdxx1_sWdy2x1=bsxfun(@times,sWdy1x1_sWdxx1,sWdy2x1_c); % sWy1x1*sWxx1_i*sWy2x1_c
R2y1y2x=sWdy1x1_sWdxx1_sWdy2x1./sWdy1y2; %(Ry,x.Z)^2
abs_one_minus_R2y1y2x=(abs(1-R2y1y2x)).^2; %squared absolute of 1-(Ry,x.Z)^2
%%%%%%%%%%%%%%%% R^2 y,x (s,t)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R2y1y2=(sWdy1y2.*conj(sWdy1y2))./(smdY1.*smdY2); %(Ry,x)^2
%%%%%%%%%%%%%%%%%%% 1-R^2 y,Z(s,t)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sWdy1x1_c=conj(sWdy1x1);
sWdy1x1_sWdxx1_sWdy1x1=bsxfun(@times,sWdy1x1_sWdxx1,sWdy1x1_c);
R2y1x=sWdy1x1_sWdxx1_sWdy1x1./smdY1; %(Ry,Z)^2
one_minus_R2y1x=(1-R2y1x);

%%%%%%%%%%%%%%%%%%%%%  1-R^2 x,Z(s,t) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sWdy2x1_sWdxx1=bsxfun(@times,sWdy2x1,sWdxx1_i);
sWdy2x1_sWdxx1_sWdy2x1=bsxfun(@times,sWdy2x1_sWdxx1,sWdy2x1_c);
R2y2x=sWdy2x1_sWdxx1_sWdy2x1./smdY2; %(Rx,Z)^2
one_minus_R2y2x=(1-R2y2x);
%%%%%%%%%%%%%%%%%%  R^2 y,x.?Z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rsq_num=roundn(bsxfun(@times,abs_one_minus_R2y1y2x,R2y1y2),-16); %numerator part of the equation for squared pwc, if the value is less than -10^16, 0 is assigned
Rsq_den=real(bsxfun(@times,one_minus_R2y1x,one_minus_R2y2x));  %denominator part of the equation for squared pwc
Rsq= bsxfun(@rdivide,Rsq_num,Rsq_den); %squared partial wavelet coherence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else

% --------------- Partial wavelet coherence for two or more excluding variables ---------------------------------

%%%%%%%%%%%%%%%%%%%%% |1-R^2 y,x?Z (s,t)|^2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sWdy1x1_p=permute(sWdy1x1,[3,1,2]); %re-structure sWy1x1
sWdy1x1_rp=reshape(sWdy1x1_p,[1,nofe,m,n]);%change 3D to 4D matrix
sWdxx1_p=permute(sWdxx1,[4,3,1,2]); % re-structure sWxx1 (smoothed auto-or cross-wavelet power spectra for excluding factors)
for jj=1:n;
    for ii=1:m;
        sWdxx1_ip(:,:,ii,jj)=inv(sWdxx1_p(:,:,ii,jj)); % inversed sWxx1
    end
end
sWdy1x1_sWdxx1=bsxfun(@times,sWdy1x1_rp,sWdxx1_ip); % sWy1x1_rp*sWxx1_ip
sWdy1x1_sWdxx1_s=sum(sWdy1x1_sWdxx1,2);
sWdy1x1_sWdxx1_ps=permute(sWdy1x1_sWdxx1_s,[2,1,3,4]);  % re-structure 
sWdy1x1_sWdxx1_sps=squeeze(sWdy1x1_sWdxx1_ps); % 4D to 3D

sWdy2x1_c=conj(sWdy2x1); %conjugate of smoothed cross-wavelet power spectra between independent variable and excluding variables
sWdy2x1_pc=permute(sWdy2x1_c,[3,1,2]); % re-structure

sWdy1x1_sWdxx1_sWdy2x1=bsxfun(@times,sWdy1x1_sWdxx1_sps,sWdy2x1_pc); % sWy1x1_rp*sWxx1_ip*sWy2x1_pc
sWdy1x1_sWdxx1_sWdy2x1_s=sum(sWdy1x1_sWdxx1_sWdy2x1,1);
sWdy1x1_sWdxx1_sWdy2x1_ss=squeeze(sWdy1x1_sWdxx1_sWdy2x1_s);

R2y1y2x=sWdy1x1_sWdxx1_sWdy2x1_ss./sWdy1y2; %(Ry,x.Z)^2
abs_one_minus_R2y1y2x=(abs(1-R2y1y2x)).^2; %squared absolute of 1-(Ry,x.Z)^2
%%%%%%%%%%%%%%%% R^2 y,x (s,t)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R2y1y2=(sWdy1y2.*conj(sWdy1y2))./(smdY1.*smdY2); %(Ry,x)^2
%%%%%%%%%%%%%%%%%%% 1-R^2 y,Z(s,t)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sWdy1x1_c=conj(sWdy1x1);
sWdy1x1_pc=permute(sWdy1x1_c,[3,1,2]);
sWdy1x1_sWdxx1_sWdy1x1=bsxfun(@times,sWdy1x1_sWdxx1_sps,sWdy1x1_pc);
sWdy1x1_sWdxx1_sWdy1x1_s=sum(sWdy1x1_sWdxx1_sWdy1x1,1);
sWdy1x1_sWdxx1_sWdy1x1_ss=squeeze(sWdy1x1_sWdxx1_sWdy1x1_s);
R2y1x=sWdy1x1_sWdxx1_sWdy1x1_ss./smdY1; %(Ry,Z)^2
one_minus_R2y1x=(1-R2y1x);

%%%%%%%%%%%%%%%%%%%%%  1-R^2 x,Z(s,t) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sWdy2x1_p=permute(sWdy2x1,[3,1,2]);
sWdy2x1_rp=reshape(sWdy2x1_p,[1,nofe,m,n]);
sWdy2x1_sWdxx1=bsxfun(@times,sWdy2x1_rp,sWdxx1_ip);
sWdy2x1_sWdxx1_s=sum(sWdy2x1_sWdxx1,2);
sWdy2x1_sWdxx1_ps=permute(sWdy2x1_sWdxx1_s,[2,1,3,4]);
sWdy2x1_sWdxx1_sps=squeeze(sWdy2x1_sWdxx1_ps);
sWdy2x1_sWdxx1_sWdy2x1=bsxfun(@times,sWdy2x1_sWdxx1_sps,sWdy2x1_pc);
sWdy2x1_sWdxx1_sWdy2x1_s=sum(sWdy2x1_sWdxx1_sWdy2x1,1);
sWdy2x1_sWdxx1_sWdy2x1_ss=squeeze(sWdy2x1_sWdxx1_sWdy2x1_s);
R2y2x=sWdy2x1_sWdxx1_sWdy2x1_ss./smdY2; %(Rx,Z)^2
one_minus_R2y2x=(1-R2y2x);
%%%%%%%%%%%%%%%%%%  R^2 y,x.?Z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rsq_num=roundn(bsxfun(@times,abs_one_minus_R2y1y2x,R2y1y2),-16); %numerator part of the equation for squared pwc, if the value is less than -10^16, 0 is assigned
Rsq_den=real(bsxfun(@times,one_minus_R2y1x,one_minus_R2y2x));  %denominator part of the equation for squared pwc
Rsq= bsxfun(@rdivide,Rsq_num,Rsq_den); %squared partial wavelet coherence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

 for s=1:maxscale
        cd=Rsq(s,find(outsidecoi(s,:)));
        cd=max(min(cd,1),0);
        cd=floor(cd*(nbins-1))+1;
        for jj=1:length(cd)
            wlc(s,cd(jj))=wlc(s,cd(jj))+1;
        end
    end
end
close(wbh);

for s=1:maxscale
    rsqy=((1:nbins)-.5)/nbins;
    ptile=wlc(s,:);
    idx=find(ptile~=0);
    ptile=ptile(idx);
    rsqy=rsqy(idx);
    ptile=cumsum(ptile);
    ptile=(ptile-.5)/ptile(end);
    sig95(s)=interp1(ptile,rsqy,.95);
end
pwcsig=[scale' sig95'];

if any(isnan(sig95))&(~warned)
    warning('Sig95 calculation failed. (Some NaNs)');
else
    try
        save(cachefilename,'mccount','checkvalues','pwcsig'); %save to a cache....
    catch
        warning(['Unable to write to cache file: ' cachefilename]);
    end
end
