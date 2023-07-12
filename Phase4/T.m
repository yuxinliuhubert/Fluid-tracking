%% function for rotating coordinate axes by a given angle alpha
    function [r2]=T(r1,alpha) %r1=[x1;y1;z1]and r2=[x2;y2;z2] describe the same physical pt but in 2 coordinate systems with
        % theta degress rotation CCW
        r2=[cos(-alpha) -sin(-alpha) 0; sin(-alpha) cos(-alpha) 0;0 0 1]*r1;
        % note that we are rotating the coordinate axes, not the location vector of the particle
    end