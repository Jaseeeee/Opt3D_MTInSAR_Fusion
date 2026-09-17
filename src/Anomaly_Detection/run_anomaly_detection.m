function outputs = run_anomaly_detection(fusion_file, asc_mtinsar_dir, desc_mtinsar_dir, output_dir, varargin)
%RUN_ANOMALY_DETECTION Perform PS-level MHT and object-level cross-track anomaly detection.
%
% This function starts from the dual-track fusion output and implements:
%   1) PS displacement-time-series classification using MHT;
%   2) building-level dominant non-noise deformation labels;
%   3) consolidation into four object-level deformation classes;
%   4) cross-track match/mismatch assessment;
%   5) anomaly screening using matched non-linear classes and persistent
%      mismatches under a minimum scatterer-count threshold.
%
% INPUTS
%   fusion_file       Path to Overlap_Fusion_Struct_Final.mat
%   asc_mtinsar_dir   Ascending-track MT-InSAR directory
%   desc_mtinsar_dir  Descending-track MT-InSAR directory
%   output_dir        Directory for anomaly-detection outputs
%
% NAME-VALUE OPTION
%   'MinScatterers'   Minimum number of attributed scatterers in each track
%                     used to define persistent cross-track mismatches.
%                     Default: 20.
%
% REQUIRED EXTERNAL FUNCTIONS
%   load_ts_using_pt_pars.m
%   mnt_test.m
%   parse_mht_results_v2.m
%
% IMPORTANT
%   mnt_test.m should implement the corrected MHT formulation in which
%   the alternative specification matrix contains only the additional
%   H/B/P components relative to the linear null model.
%
% PS-LEVEL TYPE CODES EXPECTED FROM parse_mht_results_v2
%   -1 Noise
%    0 Linear
%    1 H
%    2 B
%    3 P
%    4 H+P
%    5 B+P
%    6 H+B
%    7 H+B+P
%
% OBJECT-LEVEL CLASSES
%    1 Non-anomalous:    Linear, P
%    2 Step-like:        H, H+P
%    3 Velocity-change:  B, B+P
%    4 Complex:          H+B, H+B+P
%
% Tested with MATLAB R2025a.

    p = inputParser;
    addParameter(p, 'MinScatterers', 20, ...
        @(x) isnumeric(x) && isscalar(x) && x >= 1 && mod(x,1) == 0);
    parse(p, varargin{:});
    min_scatterers = p.Results.MinScatterers;

    if ~isfile(fusion_file)
        error('Fusion file not found: %s', fusion_file);
    end
    if ~exist(asc_mtinsar_dir, 'dir')
        error('Ascending MT-InSAR directory not found: %s', asc_mtinsar_dir);
    end
    if ~exist(desc_mtinsar_dir, 'dir')
        error('Descending MT-InSAR directory not found: %s', desc_mtinsar_dir);
    end
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    S = load(fusion_file, 'Overlap_Fusion_Struct');
    if ~isfield(S, 'Overlap_Fusion_Struct')
        error('The fusion file does not contain Overlap_Fusion_Struct.');
    end
    Overlap_Fusion_Struct = S.Overlap_Fusion_Struct;

    %% 1. Load acquisition dates and displacement time series
    asc_prm = load(fullfile(asc_mtinsar_dir, 'intf_prm.mat'), 'intf_prm');
    desc_prm = load(fullfile(desc_mtinsar_dir, 'intf_prm.mat'), 'intf_prm');

    ts_date_asc = asc_prm.intf_prm.mydate;
    ts_date_desc = desc_prm.intf_prm.mydate;

    [TS_Full_Asc, ~] = load_ts_using_pt_pars(asc_mtinsar_dir, ts_date_asc);
    [TS_Full_Desc, ~] = load_ts_using_pt_pars(desc_mtinsar_dir, ts_date_desc);

    %% 2. Collect attributed PS identifiers for both tracks
    all_ids_asc = [];
    map_building_asc = [];
    all_ids_desc = [];
    map_building_desc = [];

    for k = 1:numel(Overlap_Fusion_Struct)
        raw = Overlap_Fusion_Struct(k).Raw_Data_Pool;

        if raw.Has_Asc && ~isempty(raw.Asc_Pts_Raw)
            ids = raw.Asc_Pts_Raw(:,4);
            all_ids_asc = [all_ids_asc; ids]; %#ok<AGROW>
            map_building_asc = [map_building_asc; repmat(k, numel(ids), 1)]; %#ok<AGROW>
        end

        if raw.Has_Desc && ~isempty(raw.Desc_Pts_Raw)
            ids = raw.Desc_Pts_Raw(:,4);
            all_ids_desc = [all_ids_desc; ids]; %#ok<AGROW>
            map_building_desc = [map_building_desc; repmat(k, numel(ids), 1)]; %#ok<AGROW>
        end
    end

    validate_point_ids(all_ids_asc, size(TS_Full_Asc,1), 'ascending');
    validate_point_ids(all_ids_desc, size(TS_Full_Desc,1), 'descending');

    %% 3. PS-level MHT classification for attributed scatterers
    types_asc = [];
    desc_asc = {};
    if ~isempty(all_ids_asc)
        ts_subset = TS_Full_Asc(all_ids_asc,:);
        [mht, ~, Ha] = mnt_test(ts_subset, ts_date_asc);
        [types_asc, desc_asc] = parse_mht_results_v2(mht, Ha);
    end

    types_desc = [];
    desc_desc = {};
    if ~isempty(all_ids_desc)
        ts_subset = TS_Full_Desc(all_ids_desc,:);
        [mht, ~, Ha] = mnt_test(ts_subset, ts_date_desc);
        [types_desc, desc_desc] = parse_mht_results_v2(mht, Ha);
    end

    %% 4. Store PS-level labels and obtain dominant non-noise label per building
    for k = 1:numel(Overlap_Fusion_Struct)
        Overlap_Fusion_Struct(k).MHT_Analysis.Asc_Type = [];
        Overlap_Fusion_Struct(k).MHT_Analysis.Asc_Desc = {};
        Overlap_Fusion_Struct(k).MHT_Analysis.Desc_Type = [];
        Overlap_Fusion_Struct(k).MHT_Analysis.Desc_Desc = {};
        Overlap_Fusion_Struct(k).MHT_Analysis.Asc_Dominant = NaN;
        Overlap_Fusion_Struct(k).MHT_Analysis.Desc_Dominant = NaN;
        Overlap_Fusion_Struct(k).MHT_Analysis.Asc_Object_Class = NaN;
        Overlap_Fusion_Struct(k).MHT_Analysis.Desc_Object_Class = NaN;
    end

    if ~isempty(types_asc)
        buildings = unique(map_building_asc);
        for i = 1:numel(buildings)
            k = buildings(i);
            idx = map_building_asc == k;
            Overlap_Fusion_Struct(k).MHT_Analysis.Asc_Type = types_asc(idx);
            Overlap_Fusion_Struct(k).MHT_Analysis.Asc_Desc = desc_asc(idx);
        end
    end

    if ~isempty(types_desc)
        buildings = unique(map_building_desc);
        for i = 1:numel(buildings)
            k = buildings(i);
            idx = map_building_desc == k;
            Overlap_Fusion_Struct(k).MHT_Analysis.Desc_Type = types_desc(idx);
            Overlap_Fusion_Struct(k).MHT_Analysis.Desc_Desc = desc_desc(idx);
        end
    end

    for k = 1:numel(Overlap_Fusion_Struct)
        asc_labels = Overlap_Fusion_Struct(k).MHT_Analysis.Asc_Type;
        desc_labels = Overlap_Fusion_Struct(k).MHT_Analysis.Desc_Type;

        asc_nonnoise = asc_labels(asc_labels >= 0 & isfinite(asc_labels));
        desc_nonnoise = desc_labels(desc_labels >= 0 & isfinite(desc_labels));

        if ~isempty(asc_nonnoise)
            asc_dom = mode(asc_nonnoise);
            Overlap_Fusion_Struct(k).MHT_Analysis.Asc_Dominant = asc_dom;
            Overlap_Fusion_Struct(k).MHT_Analysis.Asc_Object_Class = ...
                map_to_four_classes(asc_dom);
        end

        if ~isempty(desc_nonnoise)
            desc_dom = mode(desc_nonnoise);
            Overlap_Fusion_Struct(k).MHT_Analysis.Desc_Dominant = desc_dom;
            Overlap_Fusion_Struct(k).MHT_Analysis.Desc_Object_Class = ...
                map_to_four_classes(desc_dom);
        end
    end

    %% 5. Cross-track object-level consistency and anomaly screening
    n = numel(Overlap_Fusion_Struct);
    valid_dual = false(n,1);
    is_match = false(n,1);
    is_matched_nonlinear = false(n,1);
    is_persistent_mismatch = false(n,1);
    n_ps_asc = zeros(n,1);
    n_ps_desc = zeros(n,1);
    class_asc = nan(n,1);
    class_desc = nan(n,1);

    for k = 1:n
        raw = Overlap_Fusion_Struct(k).Raw_Data_Pool;
        mht = Overlap_Fusion_Struct(k).MHT_Analysis;

        n_ps_asc(k) = size(raw.Asc_Pts_Raw, 1);
        n_ps_desc(k) = size(raw.Desc_Pts_Raw, 1);
        class_asc(k) = mht.Asc_Object_Class;
        class_desc(k) = mht.Desc_Object_Class;

        valid_dual(k) = raw.Has_Asc && raw.Has_Desc && ...
            isfinite(class_asc(k)) && isfinite(class_desc(k));

        if valid_dual(k)
            is_match(k) = class_asc(k) == class_desc(k);

            % Matched non-linear classes: step-like, velocity-change, complex.
            is_matched_nonlinear(k) = is_match(k) && class_asc(k) ~= 1;

            % Persistent mismatch under adequate attributed-PS sampling.
            is_persistent_mismatch(k) = ~is_match(k) && ...
                n_ps_asc(k) >= min_scatterers && ...
                n_ps_desc(k) >= min_scatterers;
        end

        Overlap_Fusion_Struct(k).MHT_Analysis.Cross_Track_Match = ...
            valid_dual(k) && is_match(k);
        Overlap_Fusion_Struct(k).MHT_Analysis.Is_Anomaly = ...
            is_matched_nonlinear(k) || is_persistent_mismatch(k);
    end

    final_anomaly = is_matched_nonlinear | is_persistent_mismatch;

    %% 6. Agreement as a function of minimum attributed-PS count
    M_values = [1, 5, 10, 20];
    N_valid = zeros(size(M_values));
    N_match = zeros(size(M_values));
    agreement = nan(size(M_values));

    for i = 1:numel(M_values)
        M = M_values(i);
        eligible = valid_dual & n_ps_asc >= M & n_ps_desc >= M;
        N_valid(i) = nnz(eligible);
        N_match(i) = nnz(eligible & is_match);
        if N_valid(i) > 0
            agreement(i) = N_match(i) / N_valid(i);
        end
    end

    %% 7. Optional PS-level MHT for unassigned points in the overlap
    % This block is retained for reporting the deformation-model
    % distribution of all overlap-region PS points. It is not used in the
    % building-level anomaly decision.
    fusion_dir = fileparts(fusion_file);

    [types_asc_unassigned, coords_asc_unassigned] = classify_unassigned( ...
        fullfile(fusion_dir, 'pts_unassigned_filter_asc.mat'), ...
        'pts_inside_overlap_asc', TS_Full_Asc, ts_date_asc);

    [types_desc_unassigned, coords_desc_unassigned] = classify_unassigned( ...
        fullfile(fusion_dir, 'pts_unassigned_filter_desc.mat'), ...
        'pts_inside_overlap_desc', TS_Full_Desc, ts_date_desc);

    %% 8. Save results
    Anomaly_Summary = struct();
    Anomaly_Summary.MinScatterers = min_scatterers;
    Anomaly_Summary.ValidDualTrackCount = nnz(valid_dual);
    Anomaly_Summary.MatchCount = nnz(valid_dual & is_match);
    Anomaly_Summary.MismatchCount = nnz(valid_dual & ~is_match);
    Anomaly_Summary.MatchedNonlinearCount = nnz(is_matched_nonlinear);
    Anomaly_Summary.PersistentMismatchCount = nnz(is_persistent_mismatch);
    Anomaly_Summary.FinalAnomalyCount = nnz(final_anomaly);

    Anomaly_Summary.ValidDualTrackIndices = find(valid_dual);
    Anomaly_Summary.MatchedNonlinearIndices = find(is_matched_nonlinear);
    Anomaly_Summary.PersistentMismatchIndices = find(is_persistent_mismatch);
    Anomaly_Summary.FinalAnomalyIndices = find(final_anomaly);

    Anomaly_Summary.AgreementThresholds = M_values;
    Anomaly_Summary.AgreementValidCounts = N_valid;
    Anomaly_Summary.AgreementMatchCounts = N_match;
    Anomaly_Summary.AgreementRates = agreement;

    Anomaly_Summary.UnassignedAscTypes = types_asc_unassigned;
    Anomaly_Summary.UnassignedAscCoords = coords_asc_unassigned;
    Anomaly_Summary.UnassignedDescTypes = types_desc_unassigned;
    Anomaly_Summary.UnassignedDescCoords = coords_desc_unassigned;

    output_file = fullfile(output_dir, 'Overlap_Fusion_Struct_MHT.mat');
    save(output_file, 'Overlap_Fusion_Struct', 'Anomaly_Summary', '-v7.3');

    outputs = struct();
    outputs.Overlap_Fusion_Struct = Overlap_Fusion_Struct;
    outputs.Anomaly_Summary = Anomaly_Summary;
    outputs.output_file = output_file;

    fprintf('Anomaly detection completed.\n');
    fprintf('Valid dual-track buildings: %d\n', Anomaly_Summary.ValidDualTrackCount);
    fprintf('Matched non-linear buildings: %d\n', Anomaly_Summary.MatchedNonlinearCount);
    fprintf('Persistent mismatches (M >= %d): %d\n', ...
        min_scatterers, Anomaly_Summary.PersistentMismatchCount);
    fprintf('Final detected anomalies: %d\n', Anomaly_Summary.FinalAnomalyCount);
end


function validate_point_ids(ids, n_rows, track_name)
    if isempty(ids)
        return;
    end
    if any(~isfinite(ids)) || any(ids < 1) || any(mod(ids,1) ~= 0) || any(ids > n_rows)
        error('Invalid %s-track point identifiers for the loaded time-series matrix.', track_name);
    end
end


function object_class = map_to_four_classes(ps_type)
%MAP_TO_FOUR_CLASSES Map the nine PS-level labels to four object classes.
    if ismember(ps_type, [0, 3])
        object_class = 1;  % non-anomalous: Linear, P
    elseif ismember(ps_type, [1, 4])
        object_class = 2;  % step-like: H, H+P
    elseif ismember(ps_type, [2, 5])
        object_class = 3;  % velocity-change: B, B+P
    elseif ismember(ps_type, [6, 7])
        object_class = 4;  % complex: H+B, H+B+P
    else
        object_class = NaN;
    end
end


function [types, coords] = classify_unassigned(file_name, variable_name, TS_Full, ts_date)
    types = [];
    coords = [];

    if ~isfile(file_name)
        return;
    end

    S = load(file_name);
    if ~isfield(S, variable_name)
        return;
    end

    pts = S.(variable_name);
    if isempty(pts)
        return;
    end

    ids = pts(:,4);
    validate_point_ids(ids, size(TS_Full,1), 'unassigned');
    coords = pts(:,1:2);

    ts_subset = TS_Full(ids,:);
    [mht, ~, Ha] = mnt_test(ts_subset, ts_date);
    [types, ~] = parse_mht_results_v2(mht, Ha);
end
