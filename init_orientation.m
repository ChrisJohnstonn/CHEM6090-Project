
%calculate each carbon velocity

result_c_velocity = zeros(length(result_filtered_simulations),2);
result_flip = zeros(length(result_filtered_simulations),1);
result_c_start_all = zeros(Simulations_Amount,6);

for i = 1:length(result_filtered_simulations)
   
   simulation = result_filtered_simulations(i,1)*Rotation_Steps;
   result_c_velocity(i,1) = sqrt(result_c1_xyz(simulation,1)^2 + result_c1_xyz(simulation,2)^2 + result_c1_xyz(simulation,3)^2);
   result_c_velocity(i,2) = sqrt(result_c2_xyz(simulation,1)^2 + result_c2_xyz(simulation,2)^2 + result_c2_xyz(simulation,3)^2);
   
end

%calculate carbon initial position in new axis

c1_start_all = result_initial_positions(:,4:6);
c2_start_all = result_initial_positions(:,7:9);

for i = 1:Simulations_Amount
   
    %calculate new axis
    midpoint = [(result_initial_positions(i,4)+result_initial_positions(i,7))/2 ...
              (result_initial_positions(i,5)+result_initial_positions(i,8))/2 ...
              (result_initial_positions(i,6)+result_initial_positions(i,9))/2];
  
    init_direction_vector = [(result_initial_positions(i,4)-result_initial_positions(i,7)) ...
                            (result_initial_positions(i,5)-result_initial_positions(i,8)) ...
                            (result_initial_positions(i,6)-result_initial_positions(i,9))];
                  
    norm_direction_vector = init_direction_vector/norm(init_direction_vector);
  
    mid_he_direction_vector = [(midpoint(1)-result_initial_positions(i,13)) ...
                              (midpoint(2)-result_initial_positions(i,14)) ...
                              (midpoint(3)-result_initial_positions(i,15))];
  
    c = dot(mid_he_direction_vector,norm_direction_vector);
    c2 = init_direction_vector * (c/norm(init_direction_vector));
 
    perp_vector =  mid_he_direction_vector - c2; 
    norm_perp_vector = perp_vector/norm(perp_vector);
    
    z_axis = cross(init_direction_vector,perp_vector);
    norm_z_axis = z_axis/norm(z_axis);
    
    %resolve midpoint -> carbon vector in new axis. can ignore rotation as
    %acetylene molecule is not moved in the rotation steps.
    
    mid_c1_direction_vector = [(midpoint(1)-result_initial_positions(i,4)) ...
                              (midpoint(2)-result_initial_positions(i,5)) ...
                              (midpoint(3)-result_initial_positions(i,6))];
    mid_c2_direction_vector = [(midpoint(1)-result_initial_positions(i,7)) ...
                              (midpoint(2)-result_initial_positions(i,8)) ...
                              (midpoint(3)-result_initial_positions(i,9))];
                          
    result_c_start_all(i,1:3) = [dot(mid_c1_direction_vector,norm_direction_vector) ...
                                 dot(mid_c1_direction_vector,norm_perp_vector) ...
                                 dot(mid_c1_direction_vector,norm_z_axis)];

    result_c_start_all(i,4:6) = [dot(mid_c2_direction_vector,norm_direction_vector) ...
                                 dot(mid_c2_direction_vector,norm_perp_vector) ...
                                 dot(mid_c2_direction_vector,norm_z_axis)];                             
end



%determine if need to flip xbin

for i = 1:length(result_filtered_simulations)
    
   if result_c_velocity(i,1) > result_c_velocity(i,2)
      result_flip(i) = 1; 
   else
      result_flip(i) = 0;
   end
end

result_initial_histogram_matrix = zeros(bins, bins);


for i = 1:length(result_filtered_simulations)
    
    midpoint_x_bin = ceil( (0 - mid_he_x_axis(1)) / mid_he_steps(1));
    %determine which bin simulation should go into
    
    %in plane
%    xbin = ceil( (result_mid_he_xyz(i * He_Atoms - (He_Atoms - 1),1) - mid_he_x_axis(1)) / mid_he_steps(1));
%    if result_flip(i) == 1
%        if xbin > (midpoint_x_bin)
%            xbin = xbin - 2*(xbin-(midpoint_x_bin));
%        else
%            xbin = xbin + 2*((midpoint_x_bin)-xbin);
%        end              
%     end
%    ybin = ceil( (result_mid_he_xyz(i * He_Atoms - (He_Atoms - 1),2) - mid_he_y_axis(1)) / mid_he_steps(2));
%    result_initial_histogram_matrix(ybin,xbin) = result_initial_histogram_matrix(ybin,xbin) + result_intensity(result_filtered_simulations(i,1));
 
    %out of plane
    if He_Atoms > 1
        for j = 1:(He_Atoms - 1)
            xbin = ceil( (result_mid_he_xyz(i * He_Atoms - (He_Atoms - (j+1)),1) - mid_he_x_axis(1)) / mid_he_steps(1));
            ybin = ceil( (result_mid_he_xyz(i * He_Atoms - (He_Atoms - (j+1)),2) - mid_he_y_axis(1)) / mid_he_steps(2));
            
            if eq(result_flip(i), 1)
               if xbin > (midpoint_x_bin)
                   xbin = xbin - 2*(xbin-(midpoint_x_bin));
                   %fprintf(string(' xbin: ') + xbin + string(' i: ') + i);
               else
                   xbin = xbin + 2*((midpoint_x_bin)-xbin);
                   %fprintf(string(' xbin: ') + xbin + string(' i: ') + i);
               end              
            end
            result_initial_histogram_matrix(ybin,xbin) = result_initial_histogram_matrix(ybin,xbin) + result_intensity(result_filtered_simulations(i,1));
        end
    end
    
    %add intensity to result matrix
    %result_initial_histogram_matrix(ybin,xbin) = result_initial_histogram_matrix(ybin,xbin) + result_intensity(result_filtered_simulations(i,1));

end

imagesc(mid_he_x_axis,mid_he_y_axis,result_initial_histogram_matrix);
xlabel('x position');
ylabel('y position');
