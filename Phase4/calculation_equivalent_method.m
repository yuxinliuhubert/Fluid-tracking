for i=1:round(NOS/N)
    alpha=-theta*(proj_used_index-1);%alpha is for tracking the degree rotated from the 1st shot
    values_this_round=proj2r0_acc(xz_proj(proj_used_index:(proj_used_index+N-1),:),theta,SRD,RDD,delta_T);
    [x0, y0, z0 ,u ,v ,w ,a_x, a_y ,a_z]=deal(values_this_round(1),values_this_round(2),values_this_round(3),values_this_round(4),values_this_round(5),values_this_round(6),values_this_round(7),values_this_round(8),values_this_round(9));
   new_positions=[x0,y0,z0];
    for j =1:N-1
        time=delta_T*j;
        new_positions=[new_positions; x0+u*time+0.5*a_x*time^2, y0+v*time+0.5*a_y*time^2, z0+w*time+0.5*a_z*time^2 ];
    end
    for k =1:N
        new_positions_rotated=T(new_positions(k,:)',alpha)';
        positions_predicted=[positions_predicted;new_positions_rotated];
    end

    proj_used_index=proj_used_index+N;
end