function [in_building_idx, unassigned_points] = buffer_points_by_buildings(geojson_file, pt_matched_file, buffer_radius)
% BUFFER_POINTS_BY_BUILDINGS 

    if nargin < 3
        buffer_radius = 0.00005; % 默认约5米
    end

    %% === 1. 读取 GeoJSON 文件 ===
    fprintf('加载 GeoJSON 文件...\n');
    data = jsondecode(fileread(geojson_file));

    %% === 2. 读取点云文件 ===
    fprintf('加载点云文件: %s\n', pt_matched_file);
    data_pts = load(pt_matched_file);
    vars = fieldnames(data_pts);
    
    % 优先读取 pts_3d (适配新流程)，否则盲读第一个变量 (适配旧习惯)
    if ismember('pts_3d', vars)
        pts = data_pts.pts_3d;
        fprintf('  -> 读取变量: "pts_3d"\n');
    else
        pts = data_pts.(vars{1});
        fprintf('  -> 读取默认变量: "%s"\n', vars{1});
    end
    
    lon = pts(:,1);
    lat = pts(:,2);

    %% === 3. 提取所有建筑的多边形 ===
    fprintf('提取建筑多边形...\n');
    building_polys = {};
    for k = 1:numel(data.features)
        f = data.features(k);
        if ~isfield(f,'geometry') || isempty(f.geometry), continue; end
        rings = extractRings(f.geometry);
        if isempty(rings), continue; end
        outer = rings{1}; % 仅取外环
        building_polys{end+1} = outer;
    end
    fprintf('共提取 %d 个建筑多边形。\n', numel(building_polys));

    %% === 4. 为每个建筑生成缓冲区 ===
    fprintf('为每个建筑生成 %.1f m 缓冲区...\n', buffer_radius*111320); 
    in_building_idx = false(size(lon));

    % 关闭 polyshape 警告
    warning('off', 'MATLAB:polyshape:repairedBySimplify');

    for i = 1:numel(building_polys)
        poly = building_polys{i};
        if isempty(poly), continue; end

        px = poly(:,1);
        py = poly(:,2);

        % === 经典逻辑：包含包围盒降级 ===
        try
            polyshp = polyshape(px, py);
            bufshp = polybuffer(polyshp, buffer_radius);
            [in, on] = inpolygon(lon, lat, bufshp.Vertices(:,1), bufshp.Vertices(:,2));
            in_building_idx = in_building_idx | in | on;
        catch
            lonmin = min(px)-buffer_radius; lonmax = max(px)+buffer_radius;
            latmin = min(py)-buffer_radius; latmax = max(py)+buffer_radius;
            in_box = lon>=lonmin & lon<=lonmax & lat>=latmin & lat<=latmax;
            in_building_idx = in_building_idx | in_box;
        end
        
        if mod(i, 2000) == 0, fprintf('  已处理 %d...\n', i); end
    end
    
    warning('on', 'MATLAB:polyshape:repairedBySimplify');

    %% === 5. 输出结果 ===
    unassigned_points = pts(~in_building_idx, :);
    
    fprintf('统计结果:\n');
    fprintf('  保留点数 (Buffer内): %d\n', sum(in_building_idx));
    fprintf('  额外点数 (Buffer外): %d\n', size(unassigned_points, 1));

    %% === 6. 可视化结果 (经典循环画法) ===
    % figure('Color','w', 'Name', 'Circular Buffer Result'); hold on;
    % 
    % % 红色：额外点
    % scatter(unassigned_points(:,1), unassigned_points(:,2), 2, [1 0 0], 'filled'); 
    % % 绿色：保留点
    % pts_in = pts(in_building_idx, :);
    % if ~isempty(pts_in)
    %     scatter(pts_in(:,1), pts_in(:,2), 2, [0 0.6 0], 'filled'); 
    % end
    % 
    % % === 恢复：使用循环绘制黑色轮廓 ===
    % fprintf('正在绘制建筑轮廓 ...\n');
    % for i = 1:numel(building_polys)
    %     poly = building_polys{i};
    %     if isempty(poly), continue; end
    %     plot(poly(:,1), poly(:,2), 'k-', 'LineWidth', 0.5);
    % end
    % 
    % legend({'额外点','建筑内点','建筑边界'}, 'Location', 'best');
    % xlabel('Longitude'); ylabel('Latitude');
    % title('点云与建筑缓冲区空间关系');
    % axis equal; grid on;
end

%% === 辅助函数 ===
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