function outputs = run_single_track_processing(mtinsar_dir, building_geojson, optical_height_tif, output_dir)
%RUN_SINGLE_TRACK_PROCESSING Core single-track Opt3D-assisted MT-InSAR processing.
%
% This function contains the core processing steps used for one SAR track:
%   1) MT-InSAR product loading and initial geocoding
%   2) edge-based global co-registration
%   3) vertical datum correction
%   4) geometric and SAR-directional scatterer-to-building filtering
%   5) duplicate-point competition
%   6) LoD1 building-height correction
%   7) scatterer repositioning onto corrected LoD1 surfaces
%
% INPUTS
%   mtinsar_dir       Directory containing the MT-InSAR products:
%                     slc_prm.mat, intf_prm.mat, and levelALL/
%   building_geojson  Path to the input LoD1 building GeoJSON file
%   optical_height_tif
%                     Path to the optical building-height raster used
%                     for edge-based global co-registration
%   output_dir        Directory for generated results
%
% OUTPUT
%   outputs           Structure containing the main processing outputs
%
% NOTE
%   This code calls project-specific helper functions that should be
%   available on the MATLAB path. See the repository README for details.
%
% Tested with MATLAB R2025a.

    if nargin ~= 4
        error(['Usage: outputs = run_single_track_processing(', ...
               'mtinsar_dir, building_geojson, optical_height_tif, output_dir)']);
    end

    if ~exist(mtinsar_dir, 'dir')
        error('MT-InSAR directory does not exist: %s', mtinsar_dir);
    end
    if ~isfile(building_geojson)
        error('Building GeoJSON does not exist: %s', building_geojson);
    end
    if ~isfile(optical_height_tif)
        error('Optical height raster does not exist: %s', optical_height_tif);
    end
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    %% Configuration
    % Global co-registration
    dh_step = 1.0;          % initial height-search step [m]
    dh_min  = -20.0;        % initial search lower bound [m]
    dh_max  = 30.0;         % initial search upper bound [m]
    dh_stop = 0.1;          % stop when search step reaches this value [m]

    % Vertical datum correction
    bin_width = 0.5;               % histogram bin width [m]
    peak_search_range = [5, 25];   % search interval for dominant low-height peak [m]
    smooth_span = 1;               % 1 = no smoothing
    refine_half_width = 1.0;       % median refinement window around the peak [m]

    % Building-association tolerances in geographic coordinates.
    % Replace these values if a metric-coordinate implementation is used.
    geo_buffer_tolerance = 5e-5;
    sar_buffer_tolerance = 3e-6;

    %% 1. Load MT-InSAR products
    load(fullfile(mtinsar_dir, 'slc_prm.mat'));   % expected to provide PRM / Pt_Pars
    load(fullfile(mtinsar_dir, 'intf_prm.mat'));  % expected to provide intf_prm

    nifg = length(intf_prm.master);
    P0 = Pt_Pars(0, 1, nifg, 3);

    level_folder = fullfile(mtinsar_dir, 'levelALL');
    P0 = P0.Pt_read_files([level_folder, filesep]);


    %% 2. Edge-based global co-registration
    H = imread(optical_height_tif);
    L_info = geotiffinfo(optical_height_tif);

    L_c_lon = L_info.CornerCoords.Lon(1);
    L_c_lat = L_info.CornerCoords.Lat(1);
    L_dlon  = L_info.PixelScale(1);
    L_dlat  = -L_info.PixelScale(2);
    W_lon   = L_info.Width;
    H_lon   = L_info.Height;

    edge_map = double(edge(H, 'log'));

    while true
        dh_list = dh_min:dh_step:dh_max;
        if isempty(dh_list)
            error('Invalid height-search interval during global co-registration.');
        end

        hgt_corr = zeros(numel(dh_list), 1);

        for i = 1:numel(dh_list)
            pt_geo_test = PS_geocode_mex( ...
                P0.pt(:,1:2), ...
                P0.hgt(:,1) + dh_list(i), ...
                PRM(intf_prm.master(1)).par, ...
                1, 1);

            im_idx_w = round((pt_geo_test(:,1) - L_c_lon) / L_dlon);
            im_idx_h = round((pt_geo_test(:,2) - L_c_lat) / L_dlat);

            valid_mask = ...
                im_idx_w >= 2 & im_idx_w <= W_lon - 1 & ...
                im_idx_h >= 2 & im_idx_h <= H_lon - 1;

            if nnz(valid_mask) < 10
                hgt_corr(i) = 0;
                continue;
            end

            idx_linear = sub2ind( ...
                size(edge_map), ...
                im_idx_h(valid_mask), ...
                im_idx_w(valid_mask));

            edge_values = edge_map(idx_linear);
            C = corrcoef(P0.hgt(valid_mask,1) + dh_list(i), edge_values);

            if numel(C) > 1 && isfinite(C(1,2))
                hgt_corr(i) = C(1,2);
            else
                hgt_corr(i) = 0;
            end
        end

        hgt_corr(~isfinite(hgt_corr)) = -1;
        [max_corr, max_idx] = max(hgt_corr);

        if isempty(max_idx) || max_corr <= 0
            warning('No valid positive correlation found; using zero height compensation.');
            dh_off0 = 0;
            break;
        end

        dh_off0 = dh_list(max_idx);

        if dh_step <= dh_stop
            break;
        end

        dh_min = dh_off0 - 3 * dh_step;
        dh_max = dh_off0 + 3 * dh_step;
        dh_step = dh_step / 10;
    end

    pt_geo_corrected = PS_geocode_mex( ...
        P0.pt(:,1:2), ...
        P0.hgt(:,1) + dh_off0, ...
        PRM(intf_prm.master(1)).par, ...
        1, 1);

    %% 3. Build corrected PS point cloud and remove vertical datum offset
    original_indices = (1:size(P0.pt, 1))';
    pts_3d = [ ...
        pt_geo_corrected, ...
        P0.hgt(:,1) + dh_off0, ...
        original_indices];

    H_raw = pts_3d(:,3);
    valid_height_mask = isfinite(H_raw);
    H_valid = H_raw(valid_height_mask);

    if isempty(H_valid)
        error('No valid PS heights are available for datum correction.');
    end

    edge_min = floor(min(H_valid) / bin_width) * bin_width;
    edge_max = ceil(max(H_valid) / bin_width) * bin_width;
    hist_edges = edge_min:bin_width:edge_max;

    [counts, hist_edges] = histcounts(H_valid, hist_edges);
    centers = hist_edges(1:end-1) + bin_width / 2;

    if smooth_span > 1
        counts_for_peak = smoothdata(counts, 'gaussian', smooth_span);
    else
        counts_for_peak = counts;
    end

    search_mask = ...
        centers >= peak_search_range(1) & ...
        centers <= peak_search_range(2);

    if ~any(search_mask)
        error('No histogram bins fall inside peak_search_range.');
    end

    centers_search = centers(search_mask);
    counts_search = counts_for_peak(search_mask);

    [~, peak_idx] = max(counts_search);
    h_peak_bin_center = centers_search(peak_idx);

    near_peak_mask = ...
        H_valid >= h_peak_bin_center - refine_half_width & ...
        H_valid <= h_peak_bin_center + refine_half_width;

    if nnz(near_peak_mask) > 100
        h_peak_offset = median(H_valid(near_peak_mask));
    else
        h_peak_offset = h_peak_bin_center;
    end

    pts_3d_peakcorr = pts_3d;
    pts_3d_peakcorr(:,3) = pts_3d(:,3) - h_peak_offset;

    ps_file = fullfile(output_dir, 'ps_datum_corrected.mat');
    save(ps_file, 'pts_3d_peakcorr', 'dh_off0', 'h_peak_offset');

    %% 4. Restrict building reference to the corrected PS coverage
    covered_buildings_file = fullfile(output_dir, 'buildings_covered.geojson');

    subset_geojson_by_point_coverage( ...
        building_geojson, ...
        pt_geo_corrected, ...
        covered_buildings_file);

    %% 5. Geometric and SAR-directional scatterer-to-building filtering
    [idx_geo, ~] = buffer_points_by_buildings( ...
        covered_buildings_file, ...
        ps_file, ...
        geo_buffer_tolerance);

    points_data_circle = pts_3d_peakcorr(idx_geo, :);
    geo_points_file = fullfile(output_dir, 'points_geometric_buffer.mat');
    save(geo_points_file, 'points_data_circle');

    [building_data_geo, ~] = export_building_pointcloud_data( ...
        covered_buildings_file, ...
        geo_points_file, ...
        geo_buffer_tolerance, ...
        fullfile(output_dir, 'building_data_geometric.mat'));

    [idx_sar, ~] = filter_points_by_sar_buffer( ...
        covered_buildings_file, ...
        ps_file, ...
        sar_buffer_tolerance, ...
        PRM(1).par);

    points_data_sar = pts_3d_peakcorr(idx_sar, :);
    sar_points_file = fullfile(output_dir, 'points_sar_buffer.mat');
    save(sar_points_file, 'points_data_sar');

    [building_data_sar, ~] = export_building_data_sar( ...
        covered_buildings_file, ...
        sar_points_file, ...
        sar_buffer_tolerance, ...
        fullfile(output_dir, 'building_data_sar.mat'), ...
        PRM(1).par);

    %% 6. Intersect the two assignments and resolve duplicate associations
    intersection_file = fullfile(output_dir, 'building_data_intersection.mat');

    [building_data_inter, ~] = generate_building_intersection_data( ...
        building_data_geo, ...
        building_data_sar, ...
        intersection_file);

    dup_info = analyze_duplicates(building_data_inter);
    building_data_norepeat = remove_duplicates_by_competition( ...
        building_data_inter, ...
        dup_info);

    %% 7. Attach LOS deformation velocities using the original PS indices
    for k = 1:numel(building_data_norepeat)
        if building_data_norepeat(k).pointcloud_count > 0
            ids = building_data_norepeat(k).pointcloud_coordinates(:,4);
            building_data_norepeat(k).deformation_velocity = P0.def_par(ids,1);
            building_data_norepeat(k).pointcloud_indices = ids;
        else
            building_data_norepeat(k).deformation_velocity = [];
            building_data_norepeat(k).pointcloud_indices = [];
        end
    end

    assigned_file = fullfile(output_dir, 'building_data_assigned.mat');
    save(assigned_file, 'building_data_norepeat');

    %% 8. Correct LoD1 building heights using attributed PS heights
    corrected_file = fullfile(output_dir, 'building_data_height_corrected.mat');

    [building_data_corrected, stats] = correct_building_heights_by_ps( ...
        building_data_norepeat, ...
        corrected_file);

    %% 9. Reposition PS points onto the corrected LoD1 surfaces
    data_attached = project_points_to_building_revision( ...
        building_data_corrected, ...
        PRM(1).par);

    final_file = fullfile(output_dir, 'single_track_output.mat');
    save(final_file, 'data_attached', 'dh_off0', 'h_peak_offset');

    %% Outputs
    outputs = struct();
    outputs.global_height_compensation = dh_off0;
    outputs.vertical_datum_offset = h_peak_offset;
    outputs.covered_buildings_file = covered_buildings_file;
    outputs.assigned_buildings = building_data_norepeat;
    outputs.corrected_buildings = building_data_corrected;
    outputs.corrected_height_statistics = stats;
    outputs.final_data = data_attached;
    outputs.final_file = final_file;

    fprintf('Single-track processing completed.\n');
    fprintf('Results saved to: %s\n', output_dir);
end
