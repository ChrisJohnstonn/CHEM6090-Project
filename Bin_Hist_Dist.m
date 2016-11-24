bins = 100;
result_histogram_matrix = zeros(bins,bins);
%find min/max xyz values for each carbon
c_xyz = vertcat(c1_xyz,c2_xyz);


c_max = max(abs(c_xyz));
c_min = min(abs(c_xyz));
c_diff = 1.2*c_max - 0.8*c_min;
c_steps = c_diff / bins;

x_axis = (0.8*c_min(1):c_steps(1):1.2*c_max(1));
z_axis = (0.8*c_min(3):c_steps(3):1.2*c_max(3));

c_xz_ke = zeros(length(c_xyz),1);

% xbin = zeros(length(c_xyz),1);
% zbin = zeros(length(c_xyz),1);
for i = 1:length(c_xyz)
%     xbin = ceil(abs((c_xyz(i,1))) - x_axis(1) / c_steps(1));
%     zbin = ceil(abs((c_xyz(i,3))) - z_axis(1) / c_steps(3));
    xbin = ceil( (abs(c_xyz(i,1)) - x_axis(1))  / c_steps(1));
    zbin = ceil( (abs(c_xyz(i,3)) - z_axis(1))  / c_steps(3));
%      disp(xbin);
%      disp(zbin);
%      disp(i);
    if i <= 480000
        result_histogram_matrix(zbin,xbin) = result_histogram_matrix(zbin,xbin) + result_intensity(ceil(i/30),1);
    else
        result_histogram_matrix(zbin,xbin) = result_histogram_matrix(zbin,xbin) + result_intensity(ceil((i-480000)/30),1); 
    end
end

for i = 1:length(c_xyz)
   c_xz_ke(i) = sqrt(c_xyz(i,1)^2 + c_xyz(i,3)^2); 
end


histfit(c_xz_ke,bins,'weibull');
figure
imagesc(x_axis,z_axis,result_histogram_matrix);


%histogram2(c1_xyz(:,1),c1_xyz(:,3),'FaceColor','flat');
%hold
%histogram2(c2_xyz(:,1),c2_xyz(:,3),'FaceColor','flat');