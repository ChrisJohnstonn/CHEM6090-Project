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

c_xz_vel = zeros(length(c_xyz),1);
c_xz_energy = zeros(length(c_xyz),1);
c_xz_angles = zeros(length(c_xyz),1);

for i = 1:length(c_xyz)
    %determine velocity in x-z plane
    c_xz_vel(i) = sqrt(c_xyz(i,1)^2 + c_xyz(i,3)^2);
    %convert to energy and convert to eV
    c_xz_energy(i) = 0.5 * m_c * c_xz_vel(i) * c_xz_vel(i) * 6.2415e18;
    %angle between velocity xz vector and x axis
    c_xz_angles(i) = atan2(c_xyz(i,1),c_xyz(i,3))*180/pi;
end
c_vel_max = max(c_xz_vel);
c_vel_min = min(c_xz_vel);
c_vel_diff = 1.2*c_vel_max - 0.8*c_vel_min;
c_vel_steps = c_vel_diff / bins;

x_vel_axis = (0.8*c_vel_min:c_vel_steps:1.2*c_vel_max);

result_dist_matrix = zeros(length(x_vel_axis),2);
result_dist_matrix(:,1) = x_vel_axis;



c_energy_min = min(c_xz_energy);
c_energy_max = max(c_xz_energy);
c_energy_diff = 1.2*c_energy_max - 0.8*c_energy_min;
c_energy_steps = c_energy_diff / bins;

x_energy_axis = (0.8*c_energy_min:c_energy_steps:1.2*c_energy_max);

result_dist_matrix_energy = zeros(length(x_energy_axis),2);
result_dist_matrix_energy(:,1) = x_energy_axis;



c_angles_min = -180;
c_angles_max = 180;
c_angles_diff = c_angles_max - c_angles_min;
c_angles_steps = c_angles_diff / bins;

x_angles_axis = (c_angles_min:c_angles_steps:c_angles_max);

result_dist_matrix_angles = zeros(length(x_angles_axis),2);
result_dist_matrix_angles(:,1) = x_angles_axis;

for i = 1:length(c_xz_vel)
    vel_x_bin = ceil( ( c_xz_vel(i) - x_vel_axis(1) ) / c_vel_steps );

    if i <= length(c_xyz)/2
        result_dist_matrix(vel_x_bin,2) = result_dist_matrix(vel_x_bin,2) + result_intensity(ceil(i/rotation),1);
    else
        result_dist_matrix(vel_x_bin,2) = result_dist_matrix(vel_x_bin,2) + result_intensity(ceil((i-(length(c_xyz)/2))/rotation),1);
    end
    
    energy_x_bin = ceil( ( c_xz_energy(i) - x_energy_axis(1) ) / c_energy_steps );
    
    if i <= length(c_xyz)/2
        result_dist_matrix_energy(energy_x_bin,2) = result_dist_matrix_energy(energy_x_bin,2) + result_intensity(ceil(i/rotation),1);
    else
        result_dist_matrix_energy(energy_x_bin,2) = result_dist_matrix_energy(energy_x_bin,2) + result_intensity(ceil((i-(length(c_xyz)/2))/rotation),1);
    end
    
    angles_x_bin = ceil( ( c_xz_angles(i) - x_angles_axis(1) ) / c_angles_steps );
    
    if i <= length(c_xyz)/2
        result_dist_matrix_angles(angles_x_bin,2) = result_dist_matrix_angles(angles_x_bin,2) + result_intensity(ceil(i/rotation),1);
    else
        result_dist_matrix_angles(angles_x_bin,2) = result_dist_matrix_angles(angles_x_bin,2) + result_intensity(ceil((i-(length(c_xyz)/2))/rotation),1);
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

figure

bar(result_dist_matrix_energy(:,1),result_dist_matrix_energy(:,2),'BarWidth',1);
xlabel('xz energy / eV');
ylabel('intensity');

figure

bar(result_dist_matrix_angles(:,1),result_dist_matrix_angles(:,2),'BarWidth',1);
xlabel('Angle between Velocity Vector and positive z axis.');
ylabel('intensity');
% figure
% 
% probplot('weibull',result_dist_matrix([35:305],2));
% pd1 = fitdist(result_dist_matrix([35:305],2),'weibull');
% figure
% probplot('ev',result_dist_matrix([35:305],2));
% pd2 = fitdist(result_dist_matrix([35:305],2),'ev');
% figure
% probplot('normal',result_dist_matrix_angles([149:352],2));
% pd3 = fitdist(result_dist_matrix([35:305],2),'normal');
% figure
% probplot('lognormal',result_dist_matrix([35:305],2));
% pd4 = fitdist(result_dist_matrix([35:305],2),'lognormal');
% 

