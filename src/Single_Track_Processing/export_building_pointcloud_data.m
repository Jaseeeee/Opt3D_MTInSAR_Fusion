function [building_data, duplicate_points] = export_building_pointcloud_data(geojson_file, pt_matched_file, buffer_radius, output_mat_file)
% EXPORT_BUILDING_POINTCLOUD_DATA 

    if nargin < 4, output_mat_file = 'building_pointcloud_data.mat'; end
    if nargin < 3, buffer_radius = 0.00005; end

    %% 读取
    fprintf('加载数据文件...\n');
    data = jsondecode(fileread(geojson_file));
    data_pts = load(pt_matched_file);
    vars = fieldnames(data_pts);
    % 智能读取
    if ismember('pts_3d_common_source', vars), pts = data_pts.pts_3d_common_source; % 优先适配主脚本里的新名字
    elseif ismember('points_3d_circle', vars), pts = data_pts.points_3d_circle;
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
    fprintf('处理建筑物数据 (四周形缓冲区)...\n');
    
    for k = 1:num_bldgs
        f = features(k);
        if ~isfield(f,'geometry') || isempty(f.geometry), continue; end
        rings = extractRings(f.geometry);
        if isempty(rings), continue; end
        
        valid_building_count = valid_building_count + 1;
        outer = rings{1};
        
        % 基础信息填充 (保持不变)
        building_data(valid_building_count).building_index = k;
        building_data(valid_building_count).original_index = k;
        if isfield(f, 'properties'), building_data(valid_building_count).properties = f.properties;
        else, building_data(valid_building_count).properties = struct(); end
        building_data(valid_building_count).geometry_type = f.geometry.type;
        building_data(valid_building_count).original_polygon = outer;
        building_data(valid_building_count).polygon_area = polyarea(outer(:,1), outer(:,2));
        building_data(valid_building_count).lon_range = [min(outer(:,1)), max(outer(:,1))];
        building_data(valid_building_count).lat_range = [min(outer(:,2)), max(outer(:,2))];
        
        % 缓冲区计算 (保持不变)
        try
            polyshp = polyshape(outer(:,1), outer(:,2));
            bufshp = polybuffer(polyshp, buffer_radius);
            buffer_poly = bufshp.Vertices;
            buffer_area_val = area(bufshp);
            method = 'exact_round';
        catch
            lonmin = min(outer(:,1))-buffer_radius; lonmax = max(outer(:,1))+buffer_radius;
            latmin = min(outer(:,2))-buffer_radius; latmax = max(outer(:,2))+buffer_radius;
            buffer_poly = [lonmin, latmin; lonmax, latmin; lonmax, latmax; lonmin, latmax; lonmin, latmin];
            buffer_area_val = (lonmax-lonmin)*(latmax-latmin);
            method = 'bounding_box';
        end
        
        building_data(valid_building_count).buffer_polygon = buffer_poly;
        building_data(valid_building_count).buffer_area = buffer_area_val;
        building_data(valid_building_count).buffer_method = method;
        
        % === 点云归属 (核心修改) ===
        [in_buf, on_buf] = inpolygon(lon, lat, buffer_poly(:,1), buffer_poly(:,2));
        mask = in_buf | on_buf;
        
        pts_extracted = pts(mask, :);
        
        % 1. 存储坐标 (Nx4)
        building_data(valid_building_count).pointcloud_coordinates = pts_extracted;
        
        % 2. 存储索引 (优先用第4列)
        if size(pts_extracted, 2) >= 4
            building_data(valid_building_count).pointcloud_indices = pts_extracted(:, 4);
        else
            building_data(valid_building_count).pointcloud_indices = find(mask);
        end
        
        building_data(valid_building_count).pointcloud_count = sum(mask);
        building_data(valid_building_count).points_inside = sum(in_buf);
        building_data(valid_building_count).points_on_boundary = sum(on_buf);
        
        % 密度与统计
        if buffer_area_val > 0
            building_data(valid_building_count).pointcloud_density = sum(mask) / buffer_area_val;
        else
            building_data(valid_building_count).pointcloud_density = 0;
        end
        
        if sum(mask) > 0
            % 3. 收集数据用于查重 (优先用ID)
            if size(pts_extracted, 2) >= 4
                all_id_collection = [all_id_collection; pts_extracted(:, 4)]; % Nx1 (IDs)
            else
                all_id_collection = [all_id_collection; pts_extracted(:, 1:3)]; % Nx3 (XYZ)
            end
            building_indices_for_points = [building_indices_for_points; repmat(k, sum(mask), 1)];
            
            % 高程统计 (用第3列)
            z_vals = pts_extracted(:,3);
            building_data(valid_building_count).elevation_stats = struct(...
                'min', min(z_vals), 'max', max(z_vals), 'mean', mean(z_vals), 'std', std(z_vals));
        else
            building_data(valid_building_count).elevation_stats = struct('min',NaN,'max',NaN,'mean',NaN,'std',NaN);
        end
        
        if mod(k, 2000) == 0, fprintf('  已处理 %d...\n', k); end
    end
    warning('on', 'MATLAB:polyshape:repairedBySimplify');
    
    %% 重复点检测 (适配 ID 或 坐标)
    duplicate_points = detect_duplicate_points_full(all_id_collection, building_indices_for_points);
    
    %% 保存
    summary_info.total_buildings = numel(features);
    summary_info.valid_buildings = valid_building_count;
    summary_info.buffer_radius = buffer_radius;
    summary_info.method = 'circular';
    save(output_mat_file, 'building_data', 'summary_info', 'duplicate_points', 'buffer_radius');
    fprintf('  [Export Circle] 导出完成 -> %s\n', output_mat_file);
end

function duplicate_info = detect_duplicate_points_full(data_col, building_indices)
    % data_col 可能是 Nx1 的 ID 列，也可能是 Nx3 的坐标列
    if isempty(data_col)
        duplicate_info = struct('coordinates',[],'occurrence_count',[], 'building_indices',{}, 'occurrence_locations',{}); return;
    end
    
    % 如果是 ID (整数)，直接用 unique 很快
    if size(data_col, 2) == 1
        [unique_vals, ~, ic] = unique(data_col);
    else
        % 如果是 坐标 (浮点)，用 uniquetol
        tolerance = 1e-8;
        try
            [unique_vals, ~, ic] = uniquetol(data_col, tolerance, 'ByRows', true, 'DataScale', 1);
        catch
            rounded = round(data_col/tolerance)*tolerance;
            [unique_vals, ~, ic] = unique(rounded, 'rows');
        end
    end
    
    unique_ic = unique(ic);
    counts = zeros(length(unique_ic), 1);
    
    % 统计出现次数 (这里可以优化速度，但保持逻辑简单先)
    % histcounts 比循环快得多
    counts = histcounts(ic, 1:(length(unique_ic)+1))';
    
    is_dup = counts > 1;
    dup_indices = unique_ic(is_dup);
    dup_counts = counts(is_dup);
    
    duplicate_info.coordinates = unique_vals(dup_indices, :); % 这里可能是ID也可能是坐标
    duplicate_info.occurrence_count = dup_counts;
    duplicate_info.building_indices = cell(length(dup_indices), 1);
    duplicate_info.occurrence_locations = cell(length(dup_indices), 1);
    
    % 这里依然需要循环来收集索引
    for i = 1:length(dup_indices)
        idx = dup_indices(i);
        locs = find(ic == idx);
        duplicate_info.building_indices{i} = building_indices(locs);
        duplicate_info.occurrence_locations{i} = locs;
    end
end
% 辅助函数保持不变 (extractRings, toNx2) ...
function rings = extractRings(geometry)
    rings = {};
    if ~isfield(geometry,'coordinates') || isempty(geometry.coordinates)
        return;
    end
    coords = geometry.coordinates;
    gtype = geometry.type;
    if strcmpi(gtype,'Polygon')
        if iscell(coords)
            for kk = 1:numel(coords)
                rings{end+1} = toNx2(coords{kk});
            end
        else
            rings{1} = toNx2(coords);
        end
    elseif strcmpi(gtype,'MultiPolygon')
        if iscell(coords)
            for p = 1:numel(coords)
                poly = coords{p};
                if iscell(poly)
                    for kk = 1:numel(poly)
                        rings{end+1} = toNx2(poly{kk});
                    end
                else
                    rings{end+1} = toNx2(poly);
                end
            end
        end
    end
end

function r2 = toNx2(r)
    if isempty(r), r2 = []; return; end
    if iscell(r)
        try
            rnum = cell2mat(r);
            r2 = toNx2(rnum);
            return;
        catch
            r2 = [];
            return;
        end
    end
    if isnumeric(r)
        if ndims(r) == 3
            if size(r,1) == 1
                r2 = squeeze(r(1,:,:)); % N x 2
            else
                r2 = reshape(r, [], size(r,3));
            end
        else
            r2 = r(:,1:2);
        end
    else
        r2 = [];
    end
end