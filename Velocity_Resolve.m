function [result_c1_xyz,result_c2_xyz,Rotation_Steps] = Velocity_Resolve(init_vector,final_vel_vector)

Simulations_Amount = length(final_vel_vector);

Rotation_Amount = 2*pi;
Rotation_Steps = 300;


result_c1_xyz = zeros(Simulations_Amount*Rotation_Steps,3);
result_c2_xyz = zeros(Simulations_Amount*Rotation_Steps,3);



    for i = 1:Simulations_Amount
        
        if mod(i,100) == 0
            disp(fprintf(string('Resolving Velocity of Simulation ') + i + string(' of ') + Simulations_Amount));
        end
        v1 = init_vector(i,:); %Initial Direction Vector
        v2_1 = final_vel_vector(i,1:3); %End Velocity Vector for c1
        v2_2 = final_vel_vector(i,4:6); %End Velocity Vector for c2
        
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
            
            
            %Calculate component of v2 in direction of v1.
            %Vector Projection
            %a = u.v where v is normalised and . is dot product
            x1 = dot(v2_1_rot,v1/norm(v1)); 

            %Find 1st perpendicular vector.

            v3 = cross(v1,v2_1);

            %calculate component along v3.

            y1 = dot(v2_1_rot,v3/norm(v3));

            %final perpendicular vector to v1 and v3.

            v4 = cross(v1,v3);

            %calculate final component - for completeness.

            z1 = dot(v2_1_rot,v4/norm(v4));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            v2_2_rot = v2_2*cos(theta) + cross(k,v2_2)*sin(theta) + k*dot(k,v2_2)*(1-cos(theta));
            %Calculate component of v2 in direction of v1.
            %Vector Projection
            %a = u.v where v is normalised and . is dot product
            x2 = dot(v2_2_rot,v1/norm(v1)); 

            %Find 1st perpendicular vector.

            v5 = cross(v1,v2_2);

            %calculate component along v3.

            y2 = dot(v2_2_rot,v5/norm(v5));

            %final perpendicular vector to v1 and v3.

            v6 = cross(v1,v5);

            %calculate final component - for completeness.

            z2 = dot(v2_2_rot,v6/norm(v6));

            %save results 

            result_c1_xyz((i-1)*Rotation_Steps + j,:) = [x1 y1 z1];
            result_c2_xyz((i-1)*Rotation_Steps + j,:) = [x2 y2 z2];
        end
    end
    
    %3D
    
%     scatter3(result_c1_xyz(:,1),result_c1_xyz(:,2),result_c1_xyz(:,3),1);
%     hold all
%     scatter3(result_c2_xyz(:,1),result_c2_xyz(:,2),result_c2_xyz(:,3),1);
    
    %2D - in plane?
    
%     scatter(result_c1_xyz(:,1),result_c1_xyz(:,3),1);
%     hold all
%     scatter(result_c2_xyz(:,1),result_c2_xyz(:,3),1);
    
%     xlabel('x vel')
%     ylabel('y vel')
%     zlabel('z vel')
end

