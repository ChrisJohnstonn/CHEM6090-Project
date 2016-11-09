function Velocity_Resolve(init_vector,final_vel_vector)

Simulations_Amount = length(final_vel_vector);

result_c1_xy = zeros(Simulations_Amount,2);
result_c2_xy = zeros(Simulations_Amount,2);

    for i = 1:Simulations_Amount
        
        v1 = init_vector(i,:); %Initial Direction Vector
        v2_1 = final_vel_vector(i,1:3); %End Velocity Vector for c1
        v2_2 = final_vel_vector(i,4:6); %End Velocity Vector for c2

        %Calculate component of v2 in direction of v1.
        %Vector Projection
        %a = u.v where v is normalised and . is dot product
        x1 = dot(v2_1,v1/norm(v1)); 

        %Find 1st perpendicular vector.

        v3 = cross(v1,v2_1);

        %calculate component along v3.

        y1 = dot(v2_1,v3/norm(v3));

        %final perpendicular vector to v1 and v3.

        v4 = cross(v1,v3);

        %calculate final component - for completeness.

        z1 = dot(v2_1,v4/norm(v4));
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        %Calculate component of v2 in direction of v1.
        %Vector Projection
        %a = u.v where v is normalised and . is dot product
        x2 = dot(v2_2,v1/norm(v1)); 

        %Find 1st perpendicular vector.

        v5 = cross(v1,v2_2);

        %calculate component along v3.

        y2 = dot(v2_2,v5/norm(v5));

        %final perpendicular vector to v1 and v3.

        v6 = cross(v1,v5);

        %calculate final component - for completeness.

        z2 = dot(v2_2,v6/norm(v6));
        
        result_c1_xy(i,:) = [x1 y1];
        result_c2_xy(i,:) = [x2 y2];

    end
    
    scatter(result_c1_xy(:,1),result_c1_xy(:,2),1);
    hold
    scatter(result_c2_xy(:,1),result_c2_xy(:,2),1);
    xlabel('x vel')
    ylabel('y vel')
end

