function [Al_thickness,PMMA_thickness, angles_out] = DeDecompose(proj_l, angles_l, proj_h, angles_h)
%DEDECOMPSE Maps CBCT projections taken at 80 and 140 kVp to equivalent
%thicknesses of Al and PMMA using precalibrated transfer functions.
%   Coefficients C and D are the coefficients for the third-degree transfer
%   functions (may need to be modified for different systems)
%
%   Proj_l and angles_l are the projection and angle data for the 80 kVp
%   scan
%
%   Proj_h and angles_h are the projection and angle data for the 140 kVp
%   scan

% calibrated coefficients for transfer functions
c = [80, -87.778, 11.875, -27.992, 16.506, -1.3622, 4.1877, -4.1797, 1.3352];
d = [-244.29, 312.13, -12.957, 36.246, -24.015, 2.3154, -7.553, 7.973, -2.7003];

%% check that the projections are the same size
for n = 1:2
    if size(proj_h,n) ~= size(proj_l,n)
        fprintf('High- and low-energy projections are not the same size!\n')
        return;
    end
end

% initialize skip counter for unmatched projections
skip_angle = 0.3; % degrees
skipcount = 0;

for k = 1:size(proj_h,3)
    if (~mod(k,50))
        fprintf("Decomposing: projection %d/%d\n", k, size(proj_h,3));
    end

    % find the corresponding angle from proj_l
    [closest, idx_l] = angle_align(angles_h(k), angles_l);
    % skip if no corresponding low-energy projection
    if (abs(closest) > (skip_angle*pi/180))
        fprintf("Alert: Angular variation between high and low energy" + ...
            " projections (%d and %d) exceeds %d degrees (%d)." + ...
            "  Skipping.\n", k, idx_l, skip_angle, abs(closest)*180/pi)
        skipcount = skipcount + 1;
        continue;
    end
    % map pixels from H and L to equivalent thicknesses
    for j = 1:size(proj_h,2)
        for i = 1:size(proj_h,1)
            H = proj_h(i,j,k);  
            L = proj_l(i,j,idx_l);
            Al_thickness(i,j,k-skipcount) = c(1)*L + c(2)*H + c(3)*(L.^2) + c(4)*L.*H + c(5)*(H.^2) + c(6)*(L.^3) + c(7)*(H.*(L.^2)) + c(8)*(L.*(H.^2)) + c(9)*(H.^3);
            PMMA_thickness(i,j,k-skipcount) = d(1)*L + d(2)*H + d(3)*(L.^2) + d(4)*L.*H + d(5)*(H.^2) + d(6)*(L.^3) + d(7)*(H.*(L.^2)) + d(8)*(L.*(H.^2)) + d(9)*(H.^3);
            angles_out(k-skipcount) = angles_h(k);
        end
    end
end
end

function [closest, idx] = angle_align(target, angle_array)

angle_diff = abs(angle_array - target);
[closest, idx] = min(angle_diff);

end