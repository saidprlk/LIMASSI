function S = assimilate_images(raw_names,sample_num,jj)
% processing of raw images to prepare them for analyses
[XX,YY,x_shr,a_shr,rect_radar_stc,T_dur,Tx]=xband_2_image_loop_epoch_time_jjfix(raw_names,sample_num,jj);
% application of radon transformation to capture wave refraction, interpolation of data to refracted
% rays 
[ex_imgs,ray_x,ray_y,ref_raw,ang_dist_fit,dir_switch,x_vec_adj] = refracted_images(x_shr,a_shr,rect_radar_stc);
%x_vec_adj = x_shr(:,1);
% application of CWT to evaluate wavenumbers
[k_max,k_filt,cele_pw,cele_fit,cele_tanh]=wavelet_k_c_f(x_vec_adj,ang_dist_fit,ex_imgs,T_dur);

grav = 9.807; 

% for power function
freq_pw = cele_pw.*k_filt;%
figure
per_hist=histogram(freq_pw);bin_avg=av2(per_hist.BinEdges')';
per_cdf = cumsum(per_hist.BinCounts);
per_cdf_sm = movmean(per_cdf,[3 3]);
per_pdf_sm = [0 diff(per_cdf_sm*(per_cdf(end)/per_cdf_sm(end)))];
per_pw=1./mean(bin_avg(per_pdf_sm==max(per_pdf_sm)));
close(gcf)

% for tanh function
freq_tanh = cele_tanh.*k_filt;%
figure
per_hist=histogram(freq_tanh);bin_avg=av2(per_hist.BinEdges')';
per_cdf = cumsum(per_hist.BinCounts);
per_cdf_sm = movmean(per_cdf,[3 3]);
per_pdf_sm = [0 diff(per_cdf_sm*(per_cdf(end)/per_cdf_sm(end)))];
per_tanh=1./mean(bin_avg(per_pdf_sm==max(per_pdf_sm)));
close(gcf)
% eliminate outlier predictions
cc_pw = rmse(cele_pw(:,1:end-1)',mean(cele_pw(:,1:end-1),1)');
pw_mask = find(cc_pw>.25);
cc_tanh = rmse(cele_tanh(:,1:end-1)',mean(cele_tanh(:,1:end-1),1)');
tanh_mask = find(cc_tanh>.85);
cele_pw(pw_mask,:) = nan;
cele_tanh(tanh_mask,:) = nan;
% minimizing index between two celerity 
divider_idx = find(mean(cele_pw-cele_tanh,1,'omitmissing')<0,1);
cele_merged = cele_tanh;
cele_merged(:,1:divider_idx) = cele_pw(:,1:divider_idx);

per_dom = per_tanh;
cele_deep = -grav*per_dom/(2*pi);
c0=abs((1./per_dom./(k_filt)));
k0_pw = abs(abs(1./per_dom)./cele_pw);%k0_pw(k0_pw<min(k_filt(:),[],'omitmissing'))=min(k_filt(:),[],'omitmissing');
k0_tanh = abs(abs(1./per_dom)./cele_tanh);
cele_shallow = sqrt(abs(grav*0.04*per_dom*cele_deep));
cdiff_mean = mean(abs(c0)-abs(cele_merged),1,'omitmissing');
% histogramming
figure
cdiff_hist=histogram(abs(c0)-abs(cele_merged),NumBins=5);bin_avg=av2(cdiff_hist.BinEdges')'; %abs(c0)-abs(cele_merged),NumBins=5
c_diff=mean(bin_avg(cdiff_hist.BinCounts==max(cdiff_hist.BinCounts)));
close(gcf)

h_ini_pw = abs(atanh((abs(cele_pw)+c_diff)/cele_deep)./(2*pi*k0_pw));
h_ini_tanh = abs(atanh((abs(cele_tanh)+c_diff)/cele_deep)./(2*pi*k0_tanh));
h_ini_merged = h_ini_tanh;
h_ini_merged(:,1:divider_idx) = h_ini_pw(:,1:divider_idx);
%h_cor1 = h_ini;
%h_cor1(abs(cele_mat)<cele_shallow) = cele_mat(abs(cele_mat)<cele_shallow).^2/grav;
load cBathy-Toolbox-master\ancona_tide_gauge_hourly.mat
t_radar =Tx(1)/86400+datenum(1970,1,1);
tidx = find(abs(t_tide-t_radar)==min(abs(t_tide-t_radar)));
tide_date=datestr(t_tide(tidx));
elev_tide = tide_elev_mvm(tidx);
h_cor=h_ini_merged-elev_tide;

k_merged = k0_tanh;
k_merged(:,1:divider_idx) = k0_pw(:,1:divider_idx);

pix_extend = 25;

mex_image = rect_radar_stc;%-mean(xbr_03_04_2023(jj).raw_image,3);
samp_img = squeeze(mex_image(round(end/2)-pix_extend:round(end/2)+pix_extend,round(end/2)-pix_extend:round(end/2)+pix_extend,:));
figure
samp_hist=histogram(samp_img);
up33idx = find(abs(cumsum(samp_hist.Values)/max(cumsum(samp_hist.Values))-.666)==min(abs(cumsum(samp_hist.Values)/max(cumsum(samp_hist.Values))-.666)));
bin_avg = av2(samp_hist.BinEdges')';
sig_max =sum(bin_avg(up33idx:end).*samp_hist.Values(up33idx:end))./sum(samp_hist.Values(up33idx:end));
close(gcf)

mex_image2 = k_merged;%-mean(xbr_03_04_2023(jj).raw_image,3);
samp_img2 = 1./mex_image2;
figure
samp_hist2=histogram(samp_img2);
up33idx2 = find(abs(cumsum(samp_hist2.Values)/max(cumsum(samp_hist2.Values))-.99)==min(abs(cumsum(samp_hist2.Values)/max(cumsum(samp_hist2.Values))-.99)));
bin_avg2 = av2(samp_hist2.BinEdges')';
wlen_99 =sum(bin_avg2(up33idx2:end).*samp_hist2.Values(up33idx2:end))./sum(samp_hist2.Values(up33idx2:end));
close(gcf)


S.x_shr = x_shr;
S.a_shr = a_shr;
S.T_dur = T_dur;
S.T_ref = Tx;
S.raw_image = rect_radar_stc;
S.ray_x= ray_x;
S.ray_y= ray_y;
S.ray_img = ex_imgs;
S.raw_dir = ref_raw;
S.dir_fit = ang_dist_fit;
S.dir_switch = dir_switch;
S.k_max = k_max;
S.k_filt = k_filt;
S.k_merged = k_merged;
S.cele_fit = cele_fit;
S.cele_pw = cele_pw;
S.cele_tanh = cele_tanh;
S.h_ini = h_ini_merged;
S.h_cor = h_cor;
S.tide_elev=elev_tide;
S.tide_date=tide_date;
S.dom_period = per_dom;
S.wlen99 = wlen_99;
S.sigma_u33 = sig_max;