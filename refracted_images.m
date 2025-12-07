% interpolation of image to refraction rays based on nonlinear fitting. 
% script is modified to account for regularly spaced rays in space
% xband_2_image_loop_epoch_time_jjfix should be run first
function [ex_imgs,ray_x,ray_y,ref_raw,ang_dist_fit,dir_switch,x_vec_adj] = refracted_images(x_shr,a_shr,rect_radar_stc)

rad_theta=(0.25:0.25:180)-0.25;
mean_prop_angle_local=nan(size(rect_radar_stc,3),1);

[m,~] = size(rect_radar_stc,[1 2]);
sub_sector_thick = 40 ;%thickness of window, direction is onshore to offshore, [pixels];
sub_sector_slide = 20; % slide amount in [pixels]
sub_sector_number = floor((m-sub_sector_thick)/sub_sector_slide+1);
sub_sector_angle=nan(sub_sector_number,size(rect_radar_stc,3));
%sub_sector_distance = nan(sub_sector_number,size(rect_radar_stc,3));
sub_sector_distance = nan(sub_sector_number,size(rect_radar_stc,3));
for img_idx=1:size(rect_radar_stc,3)
%img_idx = 1;
%img_raw = rect_radar_stc(:,:,img_idx);
img_raw_meanex = double(rect_radar_stc(:,:,img_idx))-mean(rect_radar_stc,3);


for ss=1:sub_sector_number

[rad_img,xp] = radon(img_raw_meanex((1:sub_sector_thick)+(ss-1)*sub_sector_slide,:),rad_theta);
var_rad = var(rad_img)/max(var(rad_img));
sub_sector_angle(ss,img_idx)=rad_theta(var_rad==max(var_rad))-90;
%sub_sector_angle(ss,img_idx)=sum(var_rad(var_rad>var_lim).*rad_theta(var_rad>=var_lim))/sum(var_rad(var_rad>=var_lim))-90;
sub_sector_distance(ss,img_idx) = mean(x_shr((1:sub_sector_thick)+(ss-1)*sub_sector_slide,1));
end
end
err_angle = sub_sector_angle - mean(sub_sector_angle,'all','omitmissing');
% create histogram for error distribution
figure
err_hist = histogram(err_angle(:));bin_avg=av2(err_hist.BinEdges');
% adjust the mean to eliminate outliers
mean_adj=mean(sub_sector_angle,'all','omitmissing') + mean(bin_avg(err_hist.BinCounts==max(err_hist.BinCounts)));
close(gcf)
sub_sector_angle_nan=sub_sector_angle;
sub_sector_angle_nan(abs(sub_sector_angle_nan-mean_adj)>20)=nan;
% sub_sector_angle_nan=sub_sector_angle;
% sub_sector_angle_nan(isoutlier(sub_sector_angle')')=nan;
mean_sub_sector_angle = mean(sub_sector_angle_nan,2,'omitmissing');
mean_sub_sector_angle_nan =mean_sub_sector_angle; %mean_sub_sector_angle_nan(isoutlier(mean_sub_sector_angle_nan,"median"))=nan;
ray_number = 70;
x_vec = x_shr(:,1);
if mean(mean_sub_sector_angle_nan,'omitmissing')>0
    mean_sub_sector_angle_nan([nan diff(mean_sub_sector_angle_nan)']<0)=nan;
    ang_dist_fit=fit(sub_sector_distance(~isnan(mean_sub_sector_angle_nan)),mean_sub_sector_angle_nan(~isnan(mean_sub_sector_angle_nan)),...
        'power1','lower',[-Inf 0],'upper',[Inf Inf]);
    y_vec_slope = ang_dist_fit(x_vec+1); %exp1  --> Y = a*exp(b*x)
    y_vec_first = cumsum(tand(y_vec_slope),1,"omitmissing")*mean(diff(x_vec));
    dir_switch = 1;
else
    mean_sub_sector_angle_nan([nan diff(mean_sub_sector_angle_nan)']>0)=nan;
    ang_dist_fit=fit(sub_sector_distance(~isnan(mean_sub_sector_angle_nan)),-mean_sub_sector_angle_nan(~isnan(mean_sub_sector_angle_nan)),...
        'power1','lower',[-Inf 0],'upper',[Inf Inf]);
    y_vec_slope = ang_dist_fit(x_vec+1); %exp1  --> Y = a*exp(b*x)
    y_vec_first = -cumsum(tand(y_vec_slope),1,"omitmissing")*mean(diff(x_vec));
    dir_switch = -1;
end
%ang_dist_fit=fit(sub_sector_distance(~isnan(mean_sub_sector_angle_nan)),-mean_sub_sector_angle_nan(~isnan(mean_sub_sector_angle_nan)),'power2','lower',[-Inf 0],'upper',[Inf Inf]);
% if ang_dist_fit.b<0
%     error(' a*x^b --> b < 0 --> trend is captured wrong!')
% end
ref_raw = mean(sub_sector_angle_nan,2,'omitmissing');



% create rays based on local angle 

y_vec = a_shr(1,:);
y_step = y_vec(1):diff(y_vec([1 end]))/(ray_number-1):y_vec(end);
%y_step = diff(y_vec([1 end]))/(ray_number+2);
% arrangements for equally spaced refracrion curve
sloped_x = cumsum([0;diff(x_vec)]./cosd(y_vec_slope));
sloped_x_reg = sloped_x(1):sloped_x(end)/(length(sloped_x)-1):sloped_x(end);
y_vec_slope_adj = interp1(sloped_x,y_vec_slope,sloped_x_reg);
x_vec_adj = cumsum([0;diff(sloped_x_reg)'].*cosd(y_vec_slope_adj'));

if mean(mean_sub_sector_angle_nan,'omitmissing')>0
    y_step = y_step(1:end-1);
else
    y_step = y_step(2:end);
end
for i=1:length(y_step)
    ray_x(i,:) = x_vec_adj;
    ray_y(i,:) = y_vec_first+y_step(i);
end
ray_y(ray_y>y_vec(end))=nan;ray_y(ray_y<y_vec(1))=nan;
ray_x(isnan(ray_y))=nan;



% loop for time index
ex_imgs=nan([size(ray_x) size(rect_radar_stc,3)]);

for img_idx=1:size(rect_radar_stc,3)
img_raw_meanex = rect_radar_stc(:,:,img_idx)-mean(rect_radar_stc,3);
FF = scatteredInterpolant(x_shr(:),a_shr(:),img_raw_meanex(:));

ex_imgs(:,:,img_idx) = FF(ray_x,ray_y);

end
