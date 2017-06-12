Rotation_Steps = 300;
Rotation_Amount = 2*pi;

result_c1_xyz = zeros(Simulations_Amount*Rotation_Steps,3);
result_c2_xyz = zeros(Simulations_Amount*Rotation_Steps,3);



for i = 1:Simulations_Amount
    
    if mod(i,100) == 0
          disp(fprintf(string('Resolving Velocity of Simulation ') + i + string(' of ') + Simulations_Amount));
    end
    
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
    z_axis = cross(init_direction_vector,perp_vector);
    
    c1_final_vel = result_final_velocity_vector(i,1:3);
    c2_final_vel = result_final_velocity_vector(i,4:6);
    
    v1 = init_direction_vector;
    v2_1 = c1_final_vel;
    v2_2 = c2_final_vel;
    
    for j = 1:Rotation_Steps
            
            %To simulate smearing we use Rodrigues' Rotation formula. Let v
            %be our vector which we want to rotate, let k be a unit vector
            %of which we want to rotate around and let theta be the angle
            %in radians in which we want to rotate. We then have:
            %
            %v_rot = v*cos(theta) + (k x v)*sin(theta) +
            %k(k.v)(1-cos(theta)
            %
            %Note that x represent cross multiplying and . is a dot
            %product.
            
            k = v1/norm(v1); %Normalise initial vector to be rotated around.
            theta = (j-1) * (Rotation_Amount)/(Rotation_Steps);
            v2_1_rot = v2_1*cos(theta) + cross(k,v2_1)*sin(theta) + k*dot(k,v2_1)*(1-cos(theta));
            v2_2_rot = v2_2*cos(theta) + cross(k,v2_2)*sin(theta) + k*dot(k,v2_2)*(1-cos(theta));
            
            result_c1_xyz((i-1)*Rotation_Steps + j,:,:) = [dot(v2_1_rot,norm_direction_vector) ...
                                                           dot(v2_1_rot,perp_vector/norm(perp_vector)) ...
                                                           dot(v2_1_rot,z_axis/norm(z_axis))];
                              
            result_c2_xyz((i-1)*Rotation_Steps + j,:,:) = [dot(v2_2_rot,norm_direction_vector) ...
                                                           dot(v2_2_rot,perp_vector/norm(perp_vector)) ...
                                                           dot(v2_2_rot,z_axis/norm(z_axis))];
            
        
    end
                      
                     
end

