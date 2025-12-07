function [k_mat,k_filt,cele_mat,cele_fit,cele_mat2]=wavelet_k_c_f(x_vec,ang_dist_fit,extracted_imgs,T_dur)


y_vec_slope = ang_dist_fit(x_vec+1); %exp1  --> Y = a*exp(b*x)
sloped_x = cumsum([0;diff(x_vec)]./cosd(y_vec_slope),'omitmissing');
%sloped_x_reg =(sloped_x(1):(sloped_x(end)-sloped_x(1))/(length(sloped_x)-1):sloped_x(end))';
k_wdw = [2 2];
k_mat = zeros(size(extracted_imgs,1),length(x_vec));
for i=k_wdw(1)+1:size(extracted_imgs,1)-k_wdw(2)
    eximg_wdw = reshape(permute(extracted_imgs(i-k_wdw(1):i+k_wdw(2),:,:),[2 3 1]),size(extracted_imgs,2),[]);
    for img_idx=1:size(eximg_wdw,2)
        sample_data = detrend(squeeze(eximg_wdw(:,img_idx)),'omitmissing');
        sample_data_filled = sample_data; sample_data_filled(isnan(sample_data_filled))=0;
        %sample_data_reg = interp1(sloped_x,sample_data,sloped_x_reg);
        [wlet_single(:,:,img_idx),k_wt]=cwt(sample_data_filled,1/mean(diff(sloped_x))); 
    end
    %sloped_x_reg_big = repmat(sloped_x_reg',size(k_wt,1),1);
    %sloped_x_big = repmat(sloped_x',size(k_wt,1),1);
    %k_wt_big = repmat(k_wt,1,size(sloped_x_reg,1));
    wlet_single_mean=mean(abs(wlet_single),3);
    %WLS = scatteredInterpolant(sloped_x_reg_big(:),k_wt_big(:),wlet_single_mean(:));
    %wlet_single_int = WLS(sloped_x_big,k_wt_big);
    k_wt_long = repmat(k_wt,1,length(x_vec));
    max_k_wt = k_wt_long(wlet_single_mean==max(wlet_single_mean,[],"includemissing"));
    k_mat(i,:) = max_k_wt;%movmean(max_k_wt,[10 10]);
    %figure,pcolor(sloped_x_big,k_wt_big,wlet_single_mean), shading flat
    %pause()
end
k_mat(1:k_wdw(1),:)=repmat(k_mat(k_wdw(1)+1,:),k_wdw(1),1);
k_mat(end-k_wdw(2)+1:end,:)=repmat(k_mat(end-k_wdw(2),:),k_wdw(2),1);
k_trial = k_mat;
k_trial(k_trial<0.0075)=nan;
figure
k_hist = histogram(k_trial-mean(k_trial,1,'omitmissing'));
bin_avg=av2(k_hist.BinEdges')';avg_con=numel(k_hist.Data)/k_hist.NumBins;
bin_mask = bin_avg(k_hist.Values<avg_con);bin_mask = bin_mask(bin_mask<0);
k_masked = k_trial; k_masked(k_masked<(mean(k_trial,1,'omitmissing')+bin_mask(end)))=nan;
close(gcf)
k_masked2 = k_masked;k_masked2(abs(k_masked2-mean(k_masked,1,'omitmissing'))>0.003)=nan;
k_filt = k_masked2;

dt_org = mean(diff(T_dur));
%T_enlarge = T_dur(1):dt_org/4:T_dur(end);
%fs = sample_num/T_dur(end);
rad_theta=(0.25:0.25:180)-0.25;
%cele_mat = nan(size(k_max_masked));
[~,m,n] = size(extracted_imgs);
sub_sector_thick = 40 ;%thickness of window, direction is onshore to offshore, [pixels];
sub_sector_slide = 5; % slide amount in [pixels]
sub_sector_number_x = floor((m-sub_sector_thick)/sub_sector_slide);
sub_sector_angle=nan(size(extracted_imgs,1),sub_sector_number_x);%size(extracted_imgs,1)
sub_sector_distance = nan(size(extracted_imgs,1),sub_sector_number_x);

for ray_idx=1:size(extracted_imgs,1)
    ray_timestack = squeeze(extracted_imgs(ray_idx,:,:));
    for ss=1:sub_sector_number_x
        sub_img = detrend(ray_timestack((1:sub_sector_thick)+(ss-1)*sub_sector_slide,:),'omitmissing');
        %sub_img(isnan(sub_img))=0;
        [rad_img,xp] = radon(sub_img,rad_theta);
        var_rad = var(rad_img)/max(var(rad_img));
        if isnan(var_rad)
            sub_sector_angle(ray_idx,ss)=nan;
        else
            sub_sector_angle(ray_idx,ss)=rad_theta(var_rad==max(var_rad));
        end
        sub_sector_distance(ray_idx,ss) = mean(sloped_x((1:sub_sector_thick)+(ss-1)*sub_sector_slide),'omitmissing');
    end
end
% get the error from the expected value
err_angle = sub_sector_angle - mean(sub_sector_angle,'all','omitmissing');
% create histogram for error distribution
figure
err_hist = histogram(err_angle(:));bin_avg=av2(err_hist.BinEdges');
% adjust the mean to eliminate outliers
mean_adj=mean(sub_sector_angle,'all','omitmissing') + bin_avg(err_hist.BinCounts==max(err_hist.BinCounts));
close gcf
sub_sector_angle_nan=sub_sector_angle;
sub_sector_angle_nan(abs(sub_sector_angle_nan-mean_adj)>10)=nan;
% create a fit for spatially varying wave celerity
cele_window = [2 2];
cele_mat = zeros(size(extracted_imgs,1),length(x_vec));
cele_mat2 = zeros(size(extracted_imgs,1),length(x_vec));
tanhModel = @(params, x) params(1)*tanh(params(2)*x) + params(3);

% Initial guesses [a, b, c, d]
beta0 = [1, 1, 0];
% Bounds
lb = [0, 0, 0, 0];
ub = [Inf, Inf, Inf, Inf];



for cc=cele_window(1)+1:size(extracted_imgs,1)-cele_window(2)
    sub_sector_distance_wdw = sub_sector_distance(cc-cele_window(1):cc+cele_window(2),:)*5/sloped_x(end);
    sub_sector_angle_nan_wdw = sub_sector_angle_nan(cc-cele_window(1):cc+cele_window(2),:);
    cele_fit = fit(sub_sector_distance_wdw(~isnan(sub_sector_angle_nan_wdw)),sub_sector_angle_nan_wdw(~isnan(sub_sector_angle_nan_wdw)),'power2');
    params = lsqcurvefit(tanhModel, beta0, sub_sector_distance_wdw(~isnan(sub_sector_angle_nan_wdw)), sub_sector_angle_nan_wdw(~isnan(sub_sector_angle_nan_wdw)),lb,ub);
    ray_int = cumsum(tand(cele_fit(sloped_x*5/sloped_x(end))),'omitmissing').*dt_org; %.*[(diff(sloped_x)/min(abs(diff(sloped_x))));nan]
    ray_int2 = cumsum(tand(tanhModel(params,sloped_x*5/sloped_x(end))),'omitmissing').*dt_org; %.*[(diff(sloped_x)/min(abs(diff(sloped_x))));nan]
    cele_xbr = [diff(sloped_x)./diff(ray_int);nan];
    cele_xbr2 = [diff(sloped_x)./diff(ray_int2);nan];
    cele_mat(cc,:) = cele_xbr'; %repmat(cele_xbr',sum(cele_window)+1,1)
    cele_mat2(cc,:) = cele_xbr2'; 
    % plot(sub_sector_distance_wdw,sub_sector_angle_nan_wdw,'kx')
    % hold on
    % plot(sloped_x*5/sloped_x(end),tanhModel(params,sloped_x*5/sloped_x(end)),'r-')
    % plot(sloped_x*5/sloped_x(end),cele_fit(sloped_x*5/sloped_x(end)),'b-')
    % hold off
    % ylim([145 160])
    % pause()
end
cele_mat(1:cele_window(1),:)=repmat(cele_mat(cele_window(1)+1,:),cele_window(1),1);
cele_mat(end-cele_window(2)+1:end,:)=repmat(cele_mat(end-cele_window(2),:),cele_window(2),1);
cele_mat2(1:cele_window(1),:)=repmat(cele_mat2(cele_window(1)+1,:),cele_window(1),1);
cele_mat2(end-cele_window(2)+1:end,:)=repmat(cele_mat2(end-cele_window(2),:),cele_window(2),1);
end