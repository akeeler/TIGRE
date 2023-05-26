function [VM_proj, geo, angles] = MaterialDecomposition(datafolder_l, datafolder_h, energy)
%MATERIALDECOMPOSITION computes VM projections, geometry, and angles
%
%   datafolder_l and datafolder_h = '~/path/to/varian/2023-01-01-123456/';
%   creates a set of Virtual monoenergetic CBCT projections at ENERGY from
%   80 kVp and 140 kVp raw data.
%
%   Tested using data loaded and preprocessed with VarianDataLoader, but 
%   should work with other data loaders as well.
%
%   Date: 2023-05-15
%   Author: Andrew Keeler (akeeler@luc.edu)

[proj_l, ~, angles_l] = VarianDataLoader(datafolder_l, bh=false);
for idx = 1:size(proj_l, 3)
    proj_l(:,:,idx) = imgaussfilt(proj_l(:,:,idx), 1.5);
end

[proj_h, geo, angles_h] = VarianDataLoader(datafolder_h, bh=false);
for idx = 1:size(proj_h, 3)
    proj_h(:,:,idx) = imgaussfilt(proj_h(:,:,idx), 1.5);
end

[Al_thickness, PMMA_thickness, angles] = DeDecompose(proj_l, angles_l, proj_h, angles_h);

% LUT of basis material linear attenuations (mm^-1)
switch energy
    case 20
        Al_atten = 0.1 * 3.441 * 2.7;     %  20 keV
        PMMA_atten = 0.1 * 0.5714 * 1.185;
    case 30
        Al_atten = 0.1 * 1.128 * 2.7;     %  30 keV
        PMMA_atten = 0.1 * 0.3032 * 1.185;
    case 40
        Al_atten = 0.1 * 0.5685 * 2.7;    %  40 keV
        PMMA_atten  = 0.1 * 0.2350 * 1.185;
    case 50
        Al_atten = 0.1 * 0.3681 * 2.7;    %  50 keV
        PMMA_atten = 0.1 * 0.2074 * 1.185;
    case 60
        Al_atten = 0.1 * 0.2778 * 2.7;    %  60 keV
        PMMA_atten = 0.1 * 0.1924 * 1.185;
    case 80
        Al_atten = 0.1 * 0.2018 * 2.7;    %  80 keV
        PMMA_atten = 0.1 * 0.1751 * 1.185;
    case 100
        Al_atten = 0.1 * 0.1704 * 2.7;    % 100 keV
        PMMA_atten = 0.1 * 0.1641 * 1.185;
    case 150
        Al_atten = 0.1 * 0.1378 * 2.7;    % 150 keV
        PMMA_atten = 0.1 * 0.1456 * 1.185;
    otherwise
        fprintf("Requested VM energy does not exist in the lookup table");
        Al_atten = 0;
        PMMA_atten = 0;
end

% create VM attenuation image from equivalent thicknesses 
VM_proj = Al_thickness * Al_atten + PMMA_thickness * PMMA_atten;

VM_proj = EnforcePositive(VM_proj);
VM_proj = single(VM_proj);

end