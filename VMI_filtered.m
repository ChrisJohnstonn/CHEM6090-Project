%To run after Bin_Hist_Dist_Filtered

Filtered_Simulations = length(result_filtered_simulations);

%create result matrix
result_histogram_matrix = zeros(bins,bins);

%find min/max xyz values for each carbon
c_max = max(c_xyz);
c_min = min(c_xyz);
c_diff = 1.2*c_max - 1.2*c_min; %add empty space on either side
c_steps = c_diff / bins;

%define axis limits
x_axis = (1.2*c_min(1):c_steps(1):1.2*c_max(1));
z_axis = (1.2*c_min(2):c_steps(2):1.2*c_max(2));

for i = 1:Filtered_Simulations
    
    Simulation_Number = result_filtered_simulations(i,1) - 1;
    Rotation = result_filtered_simulations(i,2) + 1;
    %determine which bin simulation should go into
    xbin = (c1_xyz(Simulation_Number*Rotation_Steps + Rotation,1) - x_axis(1)) / c_steps(1);
    zbin = (c1_xyz(Simulation_Number*Rotation_Steps + Rotation,2) - z_axis(1)) / c_steps(2);
    
    xbin2 = (c2_xyz((Simulation_Number*Rotation_Steps + Rotation),1) - x_axis(1)) / c_steps(1);
    zbin2 = (c2_xyz((Simulation_Number*Rotation_Steps + Rotation),2) - z_axis(1)) / c_steps(2);
    
    lx = floor(xbin + 1);
    ux = ceil(xbin + 1);
    b1x = x_axis(lx);
    b2x = x_axis(ux);
    lx_p = ((b2x)-(c1_xyz(Simulation_Number*Rotation_Steps + Rotation,1)))/(b2x-b1x);
    ux_p = 1 - lx_p;
    
    lz = floor(zbin + 1);
    uz = ceil(zbin + 1);
    b1z = z_axis(lz);
    b2z = z_axis(uz);
    lz_p = ((b2z)-(c1_xyz(Simulation_Number*Rotation_Steps + Rotation,2)))/(b2z-b1z);
    uz_p = 1 - lz_p;

    lx_2 = floor(xbin2 + 1);
    ux_2 = ceil(xbin2 + 1);
    b1x_2 = x_axis(lx_2);
    b2x_2 = x_axis(ux_2);
    lx_p_2 = ((b2x_2)-(c2_xyz((Simulation_Number*Rotation_Steps + Rotation),1)))/(b2x_2-b1x_2);
    ux_p_2 = 1 - lx_p_2;
    
    lz_2 = floor(zbin2 + 1);
    uz_2 = ceil(zbin2 + 1);
    b1z_2 = z_axis(lz_2);
    b2z_2 = z_axis(uz_2);
    lz_p_2 = ((b2z_2)-(c2_xyz((Simulation_Number*Rotation_Steps + Rotation),2)))/(b2z_2-b1z_2);
    uz_p_2 = 1 - lz_p_2;    
    %add intensity to result matrix
        result_histogram_matrix(lz,lx) = result_histogram_matrix(lz,lx) + lz_p*lx_p*result_intensity(Simulation_Number+1,1);
        result_histogram_matrix(uz,lx) = result_histogram_matrix(uz,lx) + uz_p*lx_p*result_intensity(Simulation_Number+1,1);
        result_histogram_matrix(lz,ux) = result_histogram_matrix(lz,ux) + lz_p*ux_p*result_intensity(Simulation_Number+1,1);
        result_histogram_matrix(uz,ux) = result_histogram_matrix(uz,ux) + uz_p*ux_p*result_intensity(Simulation_Number+1,1);
        
        result_histogram_matrix(lz_2,lx_2) = result_histogram_matrix(lz_2,lx_2) + lz_p_2*lx_p_2*result_intensity(Simulation_Number+1,1);
        result_histogram_matrix(uz_2,lx_2) = result_histogram_matrix(uz_2,lx_2) + uz_p_2*lx_p_2*result_intensity(Simulation_Number+1,1);
        result_histogram_matrix(lz_2,ux_2) = result_histogram_matrix(lz_2,ux_2) + lz_p_2*ux_p_2*result_intensity(Simulation_Number+1,1);
        result_histogram_matrix(uz_2,ux_2) = result_histogram_matrix(uz_2,ux_2) + uz_p_2*ux_p_2*result_intensity(Simulation_Number+1,1);

end

imagesc(x_axis,z_axis,result_histogram_matrix);