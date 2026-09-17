function outputs = run_dual_track_fusion(asc_dir, desc_dir, output_dir, varargin)
%RUN_DUAL_TRACK_FUSION Fuse ascending- and descending-track MT-InSAR products.
%
% Core operations:
%   1) identify the geographic overlap of the two corrected PS point clouds;
%   2) associate building records from the two tracks;
%   3) fuse single-track LoD1 height estimates;
%   4) align the additive LOS-velocity reference offset;
%   5) reproject/reposition both tracks using the fused building geometry.
%
% INPUTS
%   asc_dir     Directory containing ascending-track single-track outputs.
%   desc_dir    Directory containing descending-track single-track outputs.
%   output_dir  Directory for fusion outputs.
%
% NAME-VALUE OPTION
%   'MatchTolerance'  Maximum centroid distance used to match the same
%                     building across track-specific structures.
%                     Default: 3e-5 degrees.
%
% EXPECTED INPUT FILES
%   Point cloud (one of):
%     pts_3d_common_source.mat   -> variable pts_3d
%     ps_datum_corrected.mat     -> variable pts_3d_peakcorr
%
%   Building data (one of):
%     building_data_corrected_summary.mat
%     building_data_height_corrected.mat
%       -> variable building_data_corrected
%
%   Satellite parameters:
%     slc_prm.mat -> variable PRM
%
% OPTIONAL INPUT
%   unassigned_all_points.mat -> variable unassigned_pts
%   If available, unassigned points inside the overlap are saved for
%   subsequent PS-level anomaly analysis.
%
% OUTPUT
%   outputs.Overlap_Fusion_Struct
%   outputs.poly_overlap
%   outputs.datum_shift
%   outputs.output_file
%
% EXTERNAL DEPENDENCY
%   project_points_to_building_model.m
%
% Tested with MATLAB R2025a.

    p = inputParser;
    addParameter(p, 'MatchTolerance', 3e-5, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parse(p, varargin{:});
    match_tolerance = p.Results.MatchTolerance;

    validateattributes(asc_dir, {'char','string'}, {'nonempty'});
    validateattributes(desc_dir, {'char','string'}, {'nonempty'});
    validateattributes(output_dir, {'char','string'}, {'nonempty'});

    asc_dir = char(asc_dir);
    desc_dir = char(desc_dir);
    output_dir = char(output_dir);

    if ~exist(asc_dir, 'dir'),  error('Ascending-track directory not found: %s', asc_dir); end
    if ~exist(desc_dir, 'dir'), error('Descending-track directory not found: %s', desc_dir); end
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    %% 1. Load corrected single-track products and identify overlap
    pts_asc  = load_track_points(asc_dir);
    pts_desc = load_track_points(desc_dir);

    PRM_Asc  = load_satellite_parameters(asc_dir);
    PRM_Desc = load_satellite_parameters(desc_dir);

    k_a = convhull(pts_asc(:,1), pts_asc(:,2));
    poly_asc = polyshape(pts_asc(k_a,1), pts_asc(k_a,2));

    k_d = convhull(pts_desc(:,1), pts_desc(:,2));
    poly_desc = polyshape(pts_desc(k_d,1), pts_desc(k_d,2));

    poly_overlap = intersect(poly_asc, poly_desc);
    if area(poly_overlap) <= 0
        error('The ascending and descending point clouds have no geographic overlap.');
    end

    % Optional: retain unassigned PS points inside the overlap for later
    % PS-level deformation-model statistics.
    unassigned_asc = load_unassigned_points(asc_dir);
    unassigned_desc = load_unassigned_points(desc_dir);

    pts_inside_overlap_asc = filter_points_with_polygon(unassigned_asc, poly_overlap);
    pts_inside_overlap_desc = filter_points_with_polygon(unassigned_desc, poly_overlap);

    save(fullfile(output_dir, 'pts_unassigned_filter_asc.mat'), ...
        'pts_inside_overlap_asc');
    save(fullfile(output_dir, 'pts_unassigned_filter_desc.mat'), ...
        'pts_inside_overlap_desc');

    %% 2. Load single-track building results
    bldg_src_asc  = load_building_data(asc_dir);
    bldg_src_desc = load_building_data(desc_dir);

    % Retain ascending-reference buildings whose centroids fall in the
    % dual-track overlap. Both tracks are assumed to use the same optical
    % building reference.
    coords_asc_all = building_centroids(bldg_src_asc);
    in_overlap = isinterior(poly_overlap, coords_asc_all(:,1), coords_asc_all(:,2));
    filtered_indices = find(in_overlap);

    if isempty(filtered_indices)
        error('No building records fall inside the dual-track overlap.');
    end

    num_overlap = numel(filtered_indices);
    coords_asc_query = coords_asc_all(filtered_indices,:);
    coords_desc = building_centroids(bldg_src_desc);

    % Match building records from the two independently processed tracks.
    [idx_nearest, dist_nearest] = knnsearch(coords_desc, coords_asc_query, 'K', 1);

    %% 3. Initialize the common building-object structure
    empty_struct = struct( ...
        'Basic_Info', struct(), ...
        'Raw_Data_Pool', struct(), ...
        'Fused_Params', struct(), ...
        'Final_Geometry', struct());

    Overlap_Fusion_Struct = repmat(empty_struct, num_overlap, 1);

    for k = 1:num_overlap
        idx_A = filtered_indices(k);
        src_A = bldg_src_asc(idx_A);

        Overlap_Fusion_Struct(k).Basic_Info.Original_Index_Asc = idx_A;
        Overlap_Fusion_Struct(k).Basic_Info.Original_Index_Desc = NaN;
        Overlap_Fusion_Struct(k).Basic_Info.Polygon = src_A.original_polygon;
        Overlap_Fusion_Struct(k).Basic_Info.Properties = get_field_or_default(src_A, 'properties', struct());
        Overlap_Fusion_Struct(k).Basic_Info.Optical_Height = get_field_or_default(src_A, 'original_height', NaN);
        Overlap_Fusion_Struct(k).Basic_Info.Center = coords_asc_query(k,:);

        if isfield(src_A, 'buffer_polygon')
            Overlap_Fusion_Struct(k).Basic_Info.Asc_Buffer_Polygon = src_A.buffer_polygon;
        else
            Overlap_Fusion_Struct(k).Basic_Info.Asc_Buffer_Polygon = [];
        end

        asc_pts = get_field_or_default(src_A, 'pointcloud_coordinates', []);
        Overlap_Fusion_Struct(k).Raw_Data_Pool.Has_Asc = ~isempty(asc_pts);
        Overlap_Fusion_Struct(k).Raw_Data_Pool.Asc_Pts_Raw = asc_pts;
        Overlap_Fusion_Struct(k).Raw_Data_Pool.H_Asc_Single = ...
            get_field_or_default(src_A, 'corrected_height', NaN);

        asc_vel = get_field_or_default(src_A, 'deformation_velocity', []);
        if isempty(asc_vel)
            Overlap_Fusion_Struct(k).Raw_Data_Pool.Vel_Asc_Raw = NaN;
        else
            Overlap_Fusion_Struct(k).Raw_Data_Pool.Vel_Asc_Raw = median(asc_vel, 'omitnan');
        end

        % Descending-track record for the same optical building object.
        match_idx = idx_nearest(k);
        is_match = dist_nearest(k) <= match_tolerance;

        if is_match
            src_D = bldg_src_desc(match_idx);

            Overlap_Fusion_Struct(k).Basic_Info.Original_Index_Desc = match_idx;
            if isfield(src_D, 'buffer_polygon')
                Overlap_Fusion_Struct(k).Basic_Info.Desc_Buffer_Polygon = src_D.buffer_polygon;
            else
                Overlap_Fusion_Struct(k).Basic_Info.Desc_Buffer_Polygon = [];
            end

            desc_pts = get_field_or_default(src_D, 'pointcloud_coordinates', []);
            Overlap_Fusion_Struct(k).Raw_Data_Pool.Has_Desc = ~isempty(desc_pts);
            Overlap_Fusion_Struct(k).Raw_Data_Pool.Desc_Pts_Raw = desc_pts;
            Overlap_Fusion_Struct(k).Raw_Data_Pool.H_Desc_Single = ...
                get_field_or_default(src_D, 'corrected_height', NaN);

            desc_vel = get_field_or_default(src_D, 'deformation_velocity', []);
            if isempty(desc_vel)
                Overlap_Fusion_Struct(k).Raw_Data_Pool.Vel_Desc_Raw = NaN;
            else
                Overlap_Fusion_Struct(k).Raw_Data_Pool.Vel_Desc_Raw = median(desc_vel, 'omitnan');
            end
        else
            Overlap_Fusion_Struct(k).Basic_Info.Desc_Buffer_Polygon = [];
            Overlap_Fusion_Struct(k).Raw_Data_Pool.Has_Desc = false;
            Overlap_Fusion_Struct(k).Raw_Data_Pool.Desc_Pts_Raw = [];
            Overlap_Fusion_Struct(k).Raw_Data_Pool.H_Desc_Single = NaN;
            Overlap_Fusion_Struct(k).Raw_Data_Pool.Vel_Desc_Raw = NaN;
        end
    end

    %% 4. Fuse building height and align the relative velocity reference
    vel_asc = [];
    vel_desc = [];

    for k = 1:num_overlap
        raw = Overlap_Fusion_Struct(k).Raw_Data_Pool;

        candidates = [];
        if raw.Has_Asc && isfinite(raw.H_Asc_Single)
            candidates(end+1) = raw.H_Asc_Single; %#ok<AGROW>
        end
        if raw.Has_Desc && isfinite(raw.H_Desc_Single)
            candidates(end+1) = raw.H_Desc_Single; %#ok<AGROW>
        end

        if isempty(candidates)
            h_fused = Overlap_Fusion_Struct(k).Basic_Info.Optical_Height;
        else
            h_fused = max(candidates);
        end
        Overlap_Fusion_Struct(k).Fused_Params.Height_Unified = h_fused;

        if raw.Has_Asc && raw.Has_Desc && ...
                isfinite(raw.Vel_Asc_Raw) && isfinite(raw.Vel_Desc_Raw)
            vel_asc(end+1,1) = raw.Vel_Asc_Raw; %#ok<AGROW>
            vel_desc(end+1,1) = raw.Vel_Desc_Raw; %#ok<AGROW>
        end
    end

    % Current convention: align the descending-track velocity reference to
    % the ascending-track reference.
    if isempty(vel_asc)
        datum_shift = 0;
        warning('No valid dual-track velocity pairs were found; datum shift set to zero.');
    else
        datum_shift = median(vel_asc - vel_desc, 'omitnan');
    end

    for k = 1:num_overlap
        raw = Overlap_Fusion_Struct(k).Raw_Data_Pool;
        Overlap_Fusion_Struct(k).Fused_Params.Vel_Asc_Aligned = raw.Vel_Asc_Raw;

        if raw.Has_Desc && isfinite(raw.Vel_Desc_Raw)
            Overlap_Fusion_Struct(k).Fused_Params.Vel_Desc_Aligned = ...
                raw.Vel_Desc_Raw + datum_shift;
        else
            Overlap_Fusion_Struct(k).Fused_Params.Vel_Desc_Aligned = NaN;
        end
    end

    %% 5. Reposition both tracks using the fused building height
    Proxy_Asc = struct('pointcloud_coordinates', {}, ...
                       'corrected_height', {}, ...
                       'original_polygon', {});
    Proxy_Desc = Proxy_Asc;

    for k = 1:num_overlap
        poly = Overlap_Fusion_Struct(k).Basic_Info.Polygon;
        h_target = Overlap_Fusion_Struct(k).Fused_Params.Height_Unified;

        Proxy_Asc(k).original_polygon = poly;
        Proxy_Asc(k).corrected_height = h_target;
        Proxy_Asc(k).pointcloud_coordinates = ...
            Overlap_Fusion_Struct(k).Raw_Data_Pool.Asc_Pts_Raw;

        Proxy_Desc(k).original_polygon = poly;
        Proxy_Desc(k).corrected_height = h_target;
        Proxy_Desc(k).pointcloud_coordinates = ...
            Overlap_Fusion_Struct(k).Raw_Data_Pool.Desc_Pts_Raw;
    end

    Output_Asc = project_points_to_building_model(Proxy_Asc, PRM_Asc);
    Output_Desc = project_points_to_building_model(Proxy_Desc, PRM_Desc);

    for k = 1:num_overlap
        if ~isempty(Output_Asc(k).attached_points)
            Overlap_Fusion_Struct(k).Final_Geometry.Asc_Pts_Snapped = ...
                Output_Asc(k).attached_points;
        else
            Overlap_Fusion_Struct(k).Final_Geometry.Asc_Pts_Snapped = [];
        end

        if ~isempty(Output_Desc(k).attached_points)
            Overlap_Fusion_Struct(k).Final_Geometry.Desc_Pts_Snapped = ...
                Output_Desc(k).attached_points;
        else
            Overlap_Fusion_Struct(k).Final_Geometry.Desc_Pts_Snapped = [];
        end
    end

    %% 6. Save fusion output
    output_file = fullfile(output_dir, 'Overlap_Fusion_Struct_Final.mat');
    save(output_file, 'Overlap_Fusion_Struct', 'poly_overlap', 'datum_shift', '-v7.3');

    outputs = struct();
    outputs.Overlap_Fusion_Struct = Overlap_Fusion_Struct;
    outputs.poly_overlap = poly_overlap;
    outputs.datum_shift = datum_shift;
    outputs.output_file = output_file;

    fprintf('Dual-track fusion completed.\n');
    fprintf('Output: %s\n', output_file);
end


function pts = load_track_points(track_dir)
    candidates = { ...
        {'pts_3d_common_source.mat', 'pts_3d'}, ...
        {'ps_datum_corrected.mat', 'pts_3d_peakcorr'}};

    for i = 1:numel(candidates)
        f = fullfile(track_dir, candidates{i}{1});
        if isfile(f)
            S = load(f, candidates{i}{2});
            if isfield(S, candidates{i}{2})
                pts = S.(candidates{i}{2});
                return;
            end
        end
    end

    error('No supported corrected point-cloud file found in %s.', track_dir);
end


function prm = load_satellite_parameters(track_dir)
    f = fullfile(track_dir, 'slc_prm.mat');
    if ~isfile(f)
        error('Missing slc_prm.mat in %s.', track_dir);
    end
    S = load(f, 'PRM');
    prm = S.PRM(1).par;
end


function data = load_building_data(track_dir)
    candidates = {'building_data_corrected_summary.mat', ...
                  'building_data_height_corrected.mat'};

    for i = 1:numel(candidates)
        f = fullfile(track_dir, candidates{i});
        if isfile(f)
            S = load(f, 'building_data_corrected');
            if isfield(S, 'building_data_corrected')
                data = S.building_data_corrected;
                return;
            end
        end
    end

    error('No supported building-data file found in %s.', track_dir);
end


function pts = load_unassigned_points(track_dir)
    f = fullfile(track_dir, 'unassigned_all_points.mat');
    if ~isfile(f)
        pts = [];
        return;
    end

    S = load(f, 'unassigned_pts');
    if isfield(S, 'unassigned_pts')
        pts = S.unassigned_pts;
    else
        pts = [];
    end
end


function pts_in = filter_points_with_polygon(pts, poly)
    if isempty(pts)
        pts_in = zeros(0,4);
        return;
    end
    inside = isinterior(poly, pts(:,1), pts(:,2));
    pts_in = pts(inside,:);
end


function coords = building_centroids(buildings)
    coords = nan(numel(buildings), 2);
    for i = 1:numel(buildings)
        poly = buildings(i).original_polygon;
        coords(i,:) = mean(poly, 1, 'omitnan');
    end
end


function value = get_field_or_default(S, field_name, default_value)
    if isfield(S, field_name)
        value = S.(field_name);
    else
        value = default_value;
    end
end
