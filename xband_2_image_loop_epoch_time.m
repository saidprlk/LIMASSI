% Stacking of radar images for a chosen area on radar domain
% - Loop is added for multiple stacking
% - Time axis is corrected for adding tidal correction 

%% Rading radar data from binary files
raw_names=dir(fullfile("F:/xbr_images/raw_radar_01_03_2023/",'*imm*'));
sample_num=63;
for jj=1:numel(raw_names)/sample_num
BB=nan(1025,1025,sample_num);
for ii=(1:sample_num) % time stacking of the raw data
   % ii
fin=fopen([raw_names(1).folder '\' raw_names(ii+sample_num*(jj-1)).name]);
BB(:,:,ii)=fread(fin,[1025,Inf],'uint16');
%BB(:,:,ii)=fread(fin,[1025^2 Inf]);
BB(:,:,ii)=flip(BB(:,:,ii))';
fclose(fin);
end
BB_avg=sum(BB,3)/sample_num;
load radar_mask.mat
radar_mask=flip(radar_mask)';
%BBx=reshape(BB,1027,[]);

% Phase Correction Flag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase_cor=0;
if phase_cor
    disp("Phase Correction Active!!")
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Rectangular grid extraction 
% teta=45.4765;
% xnum=321;ynum=441;
% len_y=3300;
% len_x=2400;%2.5*len_y;
% load radar_central.mat
% X1=radar_central_meters.Position(1);
% X2=X1;%X1-len_y*cosd(teta);%+2*len_x*cosd(teta);
% X3=X1+len_y*cosd(teta);
% xx=zeros(xnum,ynum)+(X2:((X3-X2)/(ynum-1)):X3);
% dxx=0:len_x*cosd(teta)/(xnum-1):len_x*cosd(teta);
% XX=xx+dxx'-X1/2;
% 
% Y1=radar_central_meters.Position(2);
% Y2=Y1;%Y1-len_y*sind(teta);
% Y3=Y1+len_y*sind(teta);
% yy=zeros(xnum,ynum)+(Y2:(Y3-Y2)/(ynum-1):Y3);
% dyy=0:len_x*sind(teta)/(xnum-1):len_x*sind(teta);
% YY=yy+dyy';
% YY=flip(YY,2);
%% Rectangular grid extraction 
teta=45.4765;
xnum=321;ynum=441;
len_y=3300/2;
len_x=2400;%2.5*len_y;
load radar_central.mat
X1=radar_central_meters.Position(1);
X2=X1-len_y*cosd(teta);%+2*len_x*cosd(teta);
X3=X1+len_y*cosd(teta);
xx=zeros(xnum,ynum)+(X2:((X3-X2)/(ynum-1)):X3);
dxx=0:len_x*cosd(teta)/(xnum-1):len_x*cosd(teta);
XX=xx+dxx';%-X1/2;

Y1=radar_central_meters.Position(2);
Y2=Y1-len_y*sind(teta);
Y3=Y1+len_y*sind(teta);
yy=zeros(xnum,ynum)+(Y2:(Y3-Y2)/(ynum-1):Y3);
dyy=0:len_x*sind(teta)/(xnum-1):len_x*sind(teta);
YY=yy+dyy';
YY=flip(YY,2);
%% Testing of interpolated grid 
unit_pix=4500/sqrt(330^2+393^2);
unit_pix=4611.48/sqrt((513-883)^2+(513-867)^2);
[dist_x,dist_y]=meshgrid((0:1024)*unit_pix,(0:1024)*unit_pix);
CC=BB(:,:,1);
figure,
colormap("gray")
%BB(BB==0)=nan;
pcolor(dist_x,dist_y,CC.*radar_mask),shading flat
%plot(dist_x(1:25:end,1:25:end).*radar_mask(1:25:end,1:25:end),dist_y(1:25:end,1:25:end).*radar_mask(1:25:end,1:25:end),'r.',MarkerSize=15)
hold on
pp=plot(XX(1:25:end,1:25:end),YY(1:25:end,1:25:end),'r.',MarkerSize=15);
%c=colorbar;
clim([max(max(CC))*.2 max(max(CC))])
axis equal 
axis off

%% Testing rectangular dimensinal coordinate system w/ interpolated radar data
% x_shr=(0:len_x/(xnum-1):len_x)'.*ones(xnum,ynum);
% %a_shr=(-len_y:2*len_y/(ynum-1):len_y).*ones(xnum,ynum);
% a_shr=(0:len_y/(ynum-1):len_y).*ones(xnum,ynum);
x_shr=(0:len_x/(xnum-1):len_x)'.*ones(xnum,ynum);
a_shr=(-len_y:2*len_y/(ynum-1):len_y).*ones(xnum,ynum);
aa=interp2(dist_x,dist_y,BB(:,:,1),XX,YY); % interp to rect grid 
aa=flip(aa,2)/(max(max(aa))/255); % flip and normalize to gray-scale(max:255)
%aa=flip(aa,2)/(max(max(aa)*.5)/255); % flip and normalize to gray-scale(max:255)
aa(aa>=200)=255;
aa(aa<=max(max(aa))*.2)=max(max(aa))*.2;
figure,
t=tiledlayout(1,1);
nexttile
colormap("gray")
%BB(BB==0)=nan;
pcolor(x_shr,a_shr,aa),shading flat
%pcolor(x_shr,a_shr,rect_radar_stc(:,:,1)),shading flat
%c=colorbar;
%rect_radar=rect_radar/(max(max(rect_radar))/255);
%clim([max(max(rect_radar))*.2 max(max(rect_radar))])
axis equal 
%axis off
t.Padding="tight";t.TileSpacing='none';
%addpath("altmany-export_fig-3.43.0.0\")
%export_fig 'radar_bw.png' -png -transparent -r300
%saveas(gcf,['xband_stack_' datestr(str2num(raw_names(1+sample_num*(jj-1)).name(5:end))/86400+datenum(1970,1,1),'ddmmyyyy_HHMM') '.png'])
%close all
%% extracting rectangular image normalised data from raw signals and stack them
rect_radar_stc=nan(xnum,ynum,sample_num);
for ii=1:sample_num % time stacking of the raw data
aa=interp2(dist_x,dist_y,BB(:,:,ii),XX,YY); % interp to rect grid 
aa=flip(aa,2); % flip and normalize to gray-scale(max:255)
%aa=flip(aa,2)/(max(max(aa)*.5)/255); % flip and normalize to gray-scale(max:255)
%aa(aa>=255)=255;
%aa(aa<=max(max(aa))*.2)=max(max(aa))*.2;
rect_radar_stc(:,:,ii)=aa;
% pcolor(x_shr,a_shr,aa),shading flat
% axis equal 
% axis off
% colormap("gray")
% pause()
end

if phase_cor
    sp1=2.05;
    sf1=1/sp1;
    L1=sample_num;
    t1=(0:L1-1)*sp1;
    f = sf1/L1*(0:(L1/2));
    radar_fft=fft(rect_radar_stc,sample_num,3);
    radar_fft_n=radar_fft./abs(max(radar_fft,[],3));
    radar_fft_nh=radar_fft_n(:,:,1:L1/2+1);
    w = abs(sf1/L1*(-L1/2:L1/2-1))*2*pi;
    w_lrg=nan(401,401,sample_num);
    for ii=1:sample_num
        w_lrg(:,:,ii)=w(ii).*ones(401,401);
    end
    wa=2*pi/sp1;
    bin_angle=atan2d(a_shr,x_shr)+90;
    radar_fft_corr=radar_fft_n.*exp(-1i*(w_lrg/wa).*(deg2rad(bin_angle).*ones(size(rect_radar_stc))));
    radar_corr=ifft(radar_fft_corr,sample_num,3,"symmetric");
    for zz=1:sample_num
        bb=radar_corr(:,:,zz);
        bb=((bb-min(min(bb)))/(max(max(bb))-min(min(bb))))*255; % flip and normalize to gray-scale(max:255)
        radar_corr(:,:,zz)=bb;
    end
end
XYZ(:,1)=reshape(x_shr,[],1);
XYZ(:,2)=reshape(a_shr,[],1);
XYZ(:,3)=0;
if phase_cor
    RAW=squeeze(reshape(radar_corr,[],1,sample_num))';
else
    RAW=squeeze(reshape(rect_radar_stc,[],1,sample_num))';
end

Tstart=str2num(raw_names(1+sample_num*(jj-1)).name(5:end)); % days from Year (00,00,0000) 
Tend=str2num(raw_names(sample_num+sample_num*(jj-1)).name(5:end));  % days from Year (00,00,0000) 
T=Tstart:(Tend-Tstart)/(sample_num-1):Tend;  % days from Year (00,00,0000) 
CAM=ones(size(RAW,2),1);
stack_name=['cBathy-Toolbox-master\stacks\xband_stack_' ...
    datestr(str2num(raw_names(1+sample_num*(jj-1)).name(5:end))/86400+datenum(1970,1,1),'ddmmyyyy_HHMM_tt')  '_v2.mat'];
save(stack_name,'XYZ','RAW','CAM','T')
disp(['The file is saved: ' stack_name])
clear XYZ RAW CAM T

end

