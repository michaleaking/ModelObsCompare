%% example_select_roi.m
% Interactively define a region of interest (ROI) over a glacier and save it
% for reuse in the spatial-difference and flowline workflows.
%
% Two ways to make an ROI:
%   (A) draw it by hand on top of an observed velocity backdrop, or
%   (B) take the full extent of a glacier-specific observed netCDF.

cfg = setup_moc();   % put ISSM + toolbox on the path

%% (A) Draw an ROI on an observed velocity backdrop -------------------------
% Load a gridded observed velocity file and show its median-speed field.
O = moc_load_obs_netcdf('Sentinel_Subset_Upernavik.nc', 'bands',{'vv'});

% Show a mid-record frame as the backdrop, then draw a rectangle (or polygon).
%   - drag to draw, then double-click / Enter to accept
%   - passing 'name' saves rois/upernavik_box.mat and .exp
roi = moc_select_roi(O, 'shape','rectangle', 'field','vv', ...
                      'tindex', round(numel(O.time)/2), 'name','upernavik_box');

%% (B) Or just use the netCDF extent (optionally inset by a margin) ---------
roiFull = moc_roi_from_obs(O, 'margin',2000, 'name','upernavik_full');

disp('Saved ROIs:');
disp(dir(fullfile(moc_path('roi_dir'),'*.mat')));
