function [building_data, duplicate_points] = export_building_data_sar(geojson_file, pt_matched_file, tolerance, output_mat_file, sar_params)
% EXPORT_BUILDING_DATA_SAR (适配 Nx4 索引版)
% 

    if nargin < 4 || isempty(output_mat_file), output_mat_file = 'building_data_sar.mat'; end
    if nargin < 3, tolerance = 1e-5; end
    
    % 解析 SAR 参数
    if isfield(sar_params, 'par'), real_p = sar_params.par; else, real_p = sar_params; end
    sat_heading = real_p.heading;
    sat_inc_rad = real_p.inc;
    layover_azi = sat_heading - 90;
    
    fprintf('  [Export SAR] 导出模式: 定向投影 (Heading=%.2f, Inc=%.2f)\n', sat_heading, rad2deg(sat_inc_rad));

    %% 读取
    data = jsondecode(fileread(geojson_file));
    data_pts = load(pt_matched_file);
    vars = fieldnames(data_pts);
    % 智能读取 (优先适配主脚本里的变量名)
    if ismember('pts_3d_common_source', vars), pts = data_pts.pts_3d_common_source;
    elseif ismember('points_3d_sar', vars), pts = data_pts.points_3d_sar;
    elseif ismember('pts_3d', vars), pts = data_pts.pts_3d;
    else, pts = data_pts.(vars{1}); end
    
    if size(pts,2) >= 2
        lon = pts(:,1);
        lat = pts(:,2);
    else
        error('点云格式错误');
    end

    %% 初始化
    building_data = struct();
    all_id_collection = []; % 改名：用于收集ID或坐标
    building_indices_for_points = [];
    valid_building_count = 0;
    features = data.features;
    num_bldgs = numel(features);
    
    warning('off', 'MATLAB:polyshape:repairedBySimplify');
    fprintf('  [Export SAR] 正在打包 %d 个建筑物数据...\n', num_bldgs);
    
    %% 循环
    for k = 1:num_bldgs
        f = features(k);
        if ~isfield(f,'geometry') || isempty(f.geometry), continue; end
        rings = extractRings(f.geometry);
        if isempty(rings), continue; end
        
        valid_building_count = valid_building_count + 1;
        outer = rings{1};
        
        % === 1. 基础信息 ===
        building_data(valid_building_count).building_index = k;
        building_data(valid_building_count).original_index = k;
        
        if isfield(f, 'properties'), building_data(valid_building_count).properties = f.properties;
        else, building_data(valid_building_count).properties = struct(); end
        
        building_data(valid_building_count).geometry_type = f.geometry.type;
        building_data(valid_building_count).original_polygon = outer;
        building_data(valid_building_count).polygon_area = polyarea(outer(:,1), outer(:,2));
        
        % 坐标范围
        building_data(valid_building_count).lon_range = [min(outer(:,1)), max(outer(:,1))];
        building_data(valid_building_count).lat_range = [min(outer(:,2)), max(outer(:,2))];
        
        % === 2. 生成 SAR 形状 ===
        if isfield(f, 'properties') && isfield(f.properties, 'height') && f.properties.height > 0
            h_val = f.properties.height;
        else
            h_val = 30; 
        end
        
        sar_poly = compute_sar_projection_export(outer, h_val, layover_azi, sat_inc_rad);
        
        try
            ps = polyshape(sar_poly(:,1), sar_poly(:,2));
            if ps.NumRegions > 0
                ps_buf = polybuffer(ps, tolerance);
                buffer_poly = ps_buf.Vertices;
                buffer_area_val = area(ps_buf);
            else
                buffer_poly = sar_poly;
                buffer_area_val = polyarea(sar_poly(:,1), sar_poly(:,2));
            end
        catch
            buffer_poly = sar_poly;
            buffer_area_val = polyarea(sar_poly(:,1), sar_poly(:,2));
        end
        
        building_data(valid_building_count).buffer_polygon = buffer_poly;
        building_data(valid_building_count).buffer_area = buffer_area_val;
        building_data(valid_building_count).buffer_method = 'sar_directional';
        
        % === 3. 点云归属 (核心修改) ===
        if ~isempty(buffer_poly)
            [in_buf, on_buf] = inpolygon(lon, lat, buffer_poly(:,1), buffer_poly(:,2));
            mask = in_buf | on_buf;
        else
            mask = false(size(lon));
            in_buf = mask; on_buf = mask;
        end
        
        pts_extracted = pts(mask, :);
        
        % A. 存储坐标 (Nx4)
        building_data(valid_building_count).pointcloud_coordinates = pts_extracted;
        
        % B. 存储索引 (优先用第4列 ID)
        if size(pts_extracted, 2) >= 4
            building_data(valid_building_count).pointcloud_indices = pts_extracted(:, 4);
        else
            building_data(valid_building_count).pointcloud_indices = find(mask);
        end
        
        building_data(valid_building_count).pointcloud_count = sum(mask);
        
        % 详细点位统计
        building_data(valid_building_count).points_inside = sum(in_buf);
        building_data(valid_building_count).points_on_boundary = sum(on_buf);
        
        % 密度计算
        if buffer_area_val > 0
            building_data(valid_building_count).pointcloud_density = sum(mask) / buffer_area_val;
        else
            building_data(valid_building_count).pointcloud_density = 0;
        end
        
        % C. 收集数据用于判重
        if sum(mask) > 0
            % 如果有第4列，收集ID (Nx1)；否则收集XYZ (Nx3)
            if size(pts_extracted, 2) >= 4
                all_id_collection = [all_id_collection; pts_extracted(:, 4)];
            else
                all_id_collection = [all_id_collection; pts_extracted(:, 1:3)];
            end
            
            building_indices_for_points = [building_indices_for_points; repmat(k, sum(mask), 1)];
            
            % 高程统计
            z_vals = pts_extracted(:,3);
            building_data(valid_building_count).elevation_stats = struct(...
                'min', min(z_vals), 'max', max(z_vals), 'mean', mean(z_vals), 'std', std(z_vals));
        else
            building_data(valid_building_count).elevation_stats = struct('min',NaN,'max',NaN,'mean',NaN,'std',NaN);
        end
        
        if mod(k, 2000) == 0, fprintf('    已处理 %d...\n', k); end
    end
    warning('on', 'MATLAB:polyshape:repairedBySimplify');
    
    %% 重复点检测 (适配 ID 检测)
    duplicate_points = detect_duplicate_points_full(all_id_collection, building_indices_for_points);
    
    %% 保存
    summary_info.total_buildings = numel(features);
    summary_info.valid_buildings = valid_building_count;
    summary_info.method = 'sar_directional';
    save(output_mat_file, 'building_data', 'summary_info', 'duplicate_points');
    fprintf('  [Export SAR] 导出完成 -> %s\n', output_mat_file);
end

%% === 内部子函数 ===
function proj_coords = compute_sar_projection_export(base_poly, h, azi_deg, inc_rad)
    center_pt = mean(base_poly(1:end-1,:), 1);
    ref_x = center_pt(1); ref_y = center_pt(2);
    scale_lat = 111320;
    scale_lon = 111320 * cosd(ref_y); 
    enu_x = (base_poly(:,1) - ref_x) * scale_lon;
    enu_y = (base_poly(:,2) - ref_y) * scale_lat;
    offset_dist = h / tan(inc_rad);
    off_x = offset_dist * sind(azi_deg);
    off_y = offset_dist * cosd(azi_deg);
    roof_x = enu_x + off_x; roof_y = enu_y + off_y;
    all_x = [enu_x; roof_x]; all_y = [enu_y; roof_y];
    try, k = convhull(all_x, all_y); hx = all_x(k); hy = all_y(k);
    catch, hx = enu_x; hy = enu_y; end
    proj_coords = [(hx / scale_lon) + ref_x, (hy / scale_lat) + ref_y]; 
end

function duplicate_info = detect_duplicate_points_full(data_col, building_indices)
    % 适配：data_col 可能是 Nx1 的 ID (整数)，也可能是 Nx3 的坐标 (浮点)
    if isempty(data_col)
        duplicate_info = struct('coordinates',[],'occurrence_count',[], 'building_indices',{}, 'occurrence_locations',{}); return;
    end
    
    % 分支处理
    if size(data_col, 2) == 1
        % ID 模式：速度极快
        [unique_vals, ~, ic] = unique(data_col);
    else
        % 坐标模式：兼容旧数据
        tolerance = 1e-8;
        try
            [unique_vals, ~, ic] = uniquetol(data_col, tolerance, 'ByRows', true, 'DataScale', 1);
        catch
            rounded = round(data_col/tolerance)*tolerance;
            [unique_vals, ~, ic] = unique(rounded, 'rows');
        end
    end
    
    unique_ic = unique(ic);
    % 使用 histcounts 加速统计
    counts = histcounts(ic, 1:(length(unique_ic)+1))';
    
    is_dup = counts > 1;
    dup_indices = unique_ic(is_dup);
    dup_counts = counts(is_dup);
    
    duplicate_info.coordinates = unique_vals(dup_indices, :);
    duplicate_info.occurrence_count = dup_counts;
    duplicate_info.building_indices = cell(length(dup_indices), 1);
    duplicate_info.occurrence_locations = cell(length(dup_indices), 1);
    
    for i = 1:length(dup_indices)
        idx = dup_indices(i);
        locs = find(ic == idx);
        duplicate_info.building_indices{i} = building_indices(locs);
        duplicate_info.occurrence_locations{i} = locs;
    end
end

function rings = extractRings(geometry)
    rings = {};
    if ~isfield(geometry,'coordinates') || isempty(geometry.coordinates), return; end
    coords = geometry.coordinates;
    if strcmpi(geometry.type,'Polygon')
        if iscell(coords), for k=1:numel(coords), rings{end+1}=toNx2(coords{k}); end
        else, rings{1}=toNx2(coords); end
    elseif strcmpi(geometry.type,'MultiPolygon')
        if iscell(coords)
            for p=1:numel(coords)
                poly=coords{p};
                if iscell(poly), for k=1:numel(poly), rings{end+1}=toNx2(poly{k}); end
                else, rings{end+1}=toNx2(poly); end
            end
        end
    end
end

function r2 = toNx2(r)
    if isempty(r), r2=[]; return; end
    if iscell(r)
        try r2=toNx2(cell2mat(r)); catch, r2=[]; end; return;
    end
    if isnumeric(r)
        if ndims(r)==3
            if size(r,1)==1, r2=squeeze(r(1,:,:));
            else, r2=reshape(r,[],size(r,3)); end
        else, r2=r(:,1:2); end
    else, r2=[]; end
end