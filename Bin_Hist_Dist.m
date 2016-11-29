%specify number of bins.
bins = 400;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%create result matrix
result_histogram_matrix = zeros(bins,bins);

%put data of both carbons into one array
c_xyz = vertcat(c1_xyz,c2_xyz);

%find min/max xyz values for each carbon
c_max = max(c_xyz);
c_min = min(c_xyz);
c_diff = 1.2*c_max - 1.2*c_min; %add empty space on either side
c_steps = c_diff / bins;

%define axis limits
x_axis = (1.2*c_min(1):c_steps(1):1.2*c_max(1));
z_axis = (1.2*c_min(3):c_steps(3):1.2*c_max(3));


for i = 1:length(c_xyz)
    %determine which bin simulation should go into
    xbin = ceil( (c_xyz(i,1) - x_axis(1)) / c_steps(1));
    zbin = ceil( (c_xyz(i,3) - z_axis(1)) / c_steps(3));
    
    %add intensity to result matrix
    if i <= length(c_xyz)/2
        result_histogram_matrix(zbin,xbin) = result_histogram_matrix(zbin,xbin) + result_intensity(ceil(i/rotation),1);
    else
        result_histogram_matrix(zbin,xbin) = result_histogram_matrix(zbin,xbin) + result_intensity(ceil((i-(length(c_xyz)/2))/rotation),1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_xz_ke = zeros(length(c_xyz),1);

for i = 1:length(c_xyz)
    %determine velocity in x-z plane
    c_xz_ke(i) = sqrt(c_xyz(i,1)^2 + c_xyz(i,3)^2);
end
c_ke_max = max(c_xz_ke);
c_ke_min = min(c_xz_ke);
c_ke_diff = 1.2*c_ke_max - 0.8*c_ke_min;
c_ke_steps = c_ke_diff / bins;

x_ke_axis = (0.8*c_ke_min:c_ke_steps:1.2*c_ke_max);

result_dist_matrix = zeros(length(x_ke_axis),2);
result_dist_matrix(:,1) = x_ke_axis;

for i = 1:length(c_xz_ke)
    ke_x_bin = ceil( ( c_xz_ke(i) - x_ke_axis(1) ) / c_ke_steps );
%     
%     disp(ke_x_bin)
%     disp(i)
    if i <= length(c_xyz)/2
        result_dist_matrix(ke_x_bin,2) = result_dist_matrix(ke_x_bin,2) + result_intensity(ceil(i/rotation),1);
    else
        result_dist_matrix(ke_x_bin,2) = result_dist_matrix(ke_x_bin,2) + result_intensity(ceil((i-(length(c_xyz)/2))/rotation),1);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%plot graphs


imagesc(x_axis,z_axis,result_histogram_matrix);
xlabel('x velocity');
ylabel('z velocity');
set(gca,'YDir','normal')
figure
bar(result_dist_matrix(:,1),result_dist_matrix(:,2),'BarWidth',1);
xlabel('xz velocity');
ylabel('intensity');
% figure
% 
% probplot('weibull',result_dist_matrix([35:305],2));
% pd1 = fitdist(result_dist_matrix([35:305],2),'weibull');
% figure
% probplot('ev',result_dist_matrix([35:305],2));
% pd2 = fitdist(result_dist_matrix([35:305],2),'ev');
% figure
% probplot('normal',result_dist_matrix([35:305],2));
% pd3 = fitdist(result_dist_matrix([35:305],2),'normal');
% figure
% probplot('lognormal',result_dist_matrix([35:305],2));
% pd4 = fitdist(result_dist_matrix([35:305],2),'lognormal');
% 

