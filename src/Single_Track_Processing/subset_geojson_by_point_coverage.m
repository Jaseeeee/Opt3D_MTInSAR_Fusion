function filtered_geojson_file = subset_geojson_by_point_coverage(input_geojson, pt_matched, output_geojson)
% SUBSET_GEOJSON_BY_POINT_COVERAGE
% 功能：根据点云的覆盖范围（凸包），筛选出范围内的建筑物，生成新的 GeoJSON。
%
% 输入：
%   input_geojson  : 原始 GeoJSON 路径 (13753个建筑)
%   pt_matched     : 配准后的点云坐标 (Nx2 或 Nx3)
%   output_geojson : 输出的新 GeoJSON 路径
%
% 输出：
%   filtered_geojson_file : 返回输出文件的路径

    fprintf('\n-------------------------------------------------\n');
    fprintf('>>> 正在根据点云覆盖范围筛选建筑物...\n');

    % 1. 读取原始 GeoJSON
    if ischar(input_geojson) || isstring(input_geojson)
        data = jsondecode(fileread(input_geojson));
    else
        error('Input GeoJSON must be a file path.');
    end
    
    features = data.features;
    total_bldgs = numel(features);
    fprintf('  原始建筑物数量: %d\n', total_bldgs);

    % 2. 计算点云的凸包 (Convex Hull)
    % 凸包是能包围所有点的最小多边形
    x = pt_matched(:,1);
    y = pt_matched(:,2);
    
    % 为了防止边缘的楼被切掉，稍微外扩一点点凸包 (Buffer) 是很难的
    % 所以我们直接用 convhull，但判断时只要建筑物有任意一点在里面就算
    k = convhull(x, y);
    hull_x = x(k);
    hull_y = y(k);
    
    % 可视化一下覆盖范围
    figure('Name', 'Coverage Filter'); hold on;
    plot(hull_x, hull_y, 'r-', 'LineWidth', 2);
    title('Point Cloud Convex Hull (Red Line)');
    axis equal; grid on;
    
    % 3. 筛选建筑物
    keep_indices = false(total_bldgs, 1);
    
    fprintf('  正在判定归属 (Convex Hull)...\n');
    
    for i = 1:total_bldgs
        feat = features(i);
        if ~isfield(feat,'geometry') || isempty(feat.geometry), continue; end
        
        % 获取建筑物坐标 (外环)
        rings = extractRings(feat.geometry);
        if isempty(rings), continue; end
        poly = rings{1};
        
        % 判定策略：只要建筑物的中心点在凸包内，就保留
        % (或者你可以改为：只要有任意顶点在凸包内)
        center = mean(poly(1:end-1, :), 1);
        
        if inpolygon(center(1), center(2), hull_x, hull_y)
            keep_indices(i) = true;
        end
        
        if mod(i, 2000) == 0, fprintf('    已检查 %d ...\n', i); end
    end
    
    % 4. 生成新的 Feature 列表
    data.features = features(keep_indices);
    valid_count = sum(keep_indices);
    
    fprintf('  筛选后建筑物数量: %d (剔除了 %d 个)\n', valid_count, total_bldgs - valid_count);
    
    % 5. 保存新的 GeoJSON
    json_str = jsonencode(data);
    fid = fopen(output_geojson, 'w');
    if fid == -1, error('无法写入输出文件'); end
    fwrite(fid, json_str, 'char');
    fclose(fid);
    
    filtered_geojson_file = output_geojson;
    fprintf('  已保存新的 GeoJSON: %s\n', output_geojson);
    
end

%% 辅助函数
function rings = extractRings(g), rings={}; if ~isfield(g,'coordinates')||isempty(g.coordinates),return;end; c=g.coordinates; if strcmpi(g.type,'Polygon'), if iscell(c), for k=1:numel(c), rings{end+1}=toNx2(c{k});end; else, rings{1}=toNx2(c);end; elseif strcmpi(g.type,'MultiPolygon'), if iscell(c), for p=1:numel(c), poly=c{p}; if iscell(poly), for k=1:numel(poly), rings{end+1}=toNx2(poly{k});end; else, rings{end+1}=toNx2(poly);end; end; end; end; end
function r2 = toNx2(r), if isempty(r),r2=[];return;end; if iscell(r), try r2=toNx2(cell2mat(r)); catch, r2=[];end; return;end; if isnumeric(r), if ndims(r)==3, if size(r,1)==1, r2=squeeze(r(1,:,:)); else, r2=reshape(r,[],size(r,3));end; else, r2=r(:,1:2);end; else, r2=[];end; end