function building_data_clean = remove_duplicates_by_competition(building_data, duplicate_info)
% REMOVE_DUPLICATES_BY_COMPETITION 基于距离竞争消除重复点 (适配 Nx4 ID & 索引同步版)
%
% 逻辑：
%   1. 读取 analyze_duplicates 发现的重复点。
%   2. 计算该点到所有"声称拥有它"的建筑物轮廓的距离。
%   3. 判定给距离最近的建筑物。
%   4. 从其他建筑物的列表中删除该点 (同步删除 Coordinates 和 Indices)。
%
% 输入：
%   building_data  : 包含重复点的原始结构体
%   duplicate_info : analyze_duplicates 的输出结果
%
% 输出：
%   building_data_clean : 去重后的干净结构体

    fprintf('\n>>> 开始基于距离竞争去除重复点 (Cleaning)...\n');
    building_data_clean = building_data;
    
    % 检查输入是否有数据
    if isempty(duplicate_info) || isempty(duplicate_info.coordinates)
        fprintf('  无重复点，无需处理。\n');
        return;
    end
    
    num_dups = size(duplicate_info.coordinates, 1);
    removed_count = 0;
    
    % 进度条
    fprintf('  正在仲裁 %d 个冲突点...\n', num_dups);
    
    for i = 1:num_dups
        % 1. 获取当前重复点
        pt = duplicate_info.coordinates(i, :); 
        pt_xy = pt(1:2);
        
        % 2. 获取涉及的建筑物索引
        bldg_indices = duplicate_info.building_indices{i}; 
        
        if isempty(bldg_indices), continue; end
        
        % 3. 计算到每个建筑物的距离
        dists = inf(length(bldg_indices), 1);
        
        for k = 1:length(bldg_indices)
            b_idx = bldg_indices(k);
            poly = building_data(b_idx).original_polygon;
            
            % 计算点到多边形(轮廓)的最短距离
            d = point_to_poly_distance_robust(pt_xy, poly);
            dists(k) = d;
        end
        
        % 4. 决出胜负 (Winner: 距离最近者)
        [~, min_idx] = min(dists);
        % winner_bldg_idx = bldg_indices(min_idx); % 赢家保留，不做操作
        
        % 5. 从输家 (Losers) 中删除该点
        for k = 1:length(bldg_indices)
            if k == min_idx, continue; end % 跳过赢家
            
            loser_bldg_idx = bldg_indices(k);
            
            % 获取输家的数据
            loser_pts = building_data_clean(loser_bldg_idx).pointcloud_coordinates;
            
            if isempty(loser_pts), continue; end
            
            % === 关键修改 A: 优先使用 ID 匹配 (第4列) ===
            match_mask = false(size(loser_pts, 1), 1);
            
            if size(pt, 2) >= 4 && size(loser_pts, 2) >= 4
                % ID 匹配 (绝对精准)
                match_mask = (loser_pts(:, 4) == pt(4));
            else
                % 降级：坐标匹配 (XYZ)
                match_mask = abs(loser_pts(:,1) - pt(1)) < 1e-8 & ...
                             abs(loser_pts(:,2) - pt(2)) < 1e-8;
                if size(pt,2) >= 3 && size(loser_pts,2) >= 3
                     match_mask = match_mask & (abs(loser_pts(:,3) - pt(3)) < 1e-8);
                end
            end
            
            if any(match_mask)
                % === 关键修改 B: 同步删除 Coordinates 和 Indices ===
                
                % 1. 删除坐标行
                loser_pts(match_mask, :) = [];
                building_data_clean(loser_bldg_idx).pointcloud_coordinates = loser_pts;
                
                % 2. 删除索引行 (如果存在)
                % 必须同步删除，否则 pointcloud_indices 和 coordinates 长度就不对应了
                if isfield(building_data_clean, 'pointcloud_indices') && ...
                   ~isempty(building_data_clean(loser_bldg_idx).pointcloud_indices)
               
                    current_indices = building_data_clean(loser_bldg_idx).pointcloud_indices;
                    % 确保长度匹配再删，防止报错
                    if length(current_indices) == length(match_mask) || size(current_indices,1) == length(match_mask)
                        current_indices(match_mask) = [];
                        building_data_clean(loser_bldg_idx).pointcloud_indices = current_indices;
                    end
                end
                
                % 3. 更新计数
                building_data_clean(loser_bldg_idx).pointcloud_count = size(loser_pts, 1);
                removed_count = removed_count + 1;
            end
        end
        
        if mod(i, 500) == 0, fprintf('    已处理 %d / %d ...\n', i, num_dups); end
    end
    
    fprintf('=== 去重完成 ===\n');
    fprintf('  共移除冗余归属: %d 次\n', removed_count);
    
    % 6. 重新计算 elevation_stats
    fprintf('  正在更新统计信息...\n');
    for k = 1:length(building_data_clean)
        if building_data_clean(k).pointcloud_count > 0
            z = building_data_clean(k).pointcloud_coordinates(:,3);
            building_data_clean(k).elevation_stats.min = min(z);
            building_data_clean(k).elevation_stats.max = max(z);
            building_data_clean(k).elevation_stats.mean = mean(z);
            building_data_clean(k).elevation_stats.std = std(z);
        else
            building_data_clean(k).elevation_stats = struct('min',nan,'max',nan,'mean',nan,'std',nan);
        end
    end

end % <--- 之前可能缺了这个 end，或者某个循环的 end 导致匹配错误

%% === 辅助函数：点到多边形距离 (绝对鲁棒版) ===
function min_d = point_to_poly_distance_robust(p, poly)
    % p: [1x2] double, poly: [Nx2] double
    
    % 1. 强制类型转换
    if iscell(poly), poly = cell2mat(poly); end
    p = double(p); poly = double(poly);
    
    % 2. 强制 p 为行向量
    if size(p, 1) > 1, p = p'; end
    
    % 3. 提取线段
    if size(poly, 1) < 2, min_d = 0; return; end
    
    p_start = poly(1:end-1, :);
    p_end = poly(2:end, :);
    
    % 4. 计算向量
    v = p_end - p_start;
    N = size(p_start, 1);
    
    % 手动扩充 p，避免使用隐式扩展
    p_rep = p(ones(N,1), :); % 等同于 repmat(p, N, 1)
    w = p_rep - p_start;
    
    % 5. 计算投影因子
    c1 = sum(w .* v, 2);
    c2 = sum(v .* v, 2);
    c2(c2 < 1e-12) = 1e-12; % 防除零
    
    b = c1 ./ c2;
    b = max(0, min(1, b)); % 限制范围
    
    % 6. 计算投影点
    b_rep = [b, b]; % 扩充为 [N, 2]
    proj = p_start + b_rep .* v;
    
    % 7. 计算距离
    dists_sq = sum((p_rep - proj).^2, 2);
    min_d = min(sqrt(dists_sq));
    
    % 8. 内部检查
    if inpolygon(p(1), p(2), poly(:,1), poly(:,2))
        min_d = 0;
    end
end