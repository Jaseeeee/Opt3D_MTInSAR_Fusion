function duplicate_info = analyze_duplicates(building_data)
% ANALYZE_DUPLICATES 分析建筑物数据中的点云重复情况 (适配 ID 查重版)
%
% 功能：
%   遍历 building_data 结构体，检测重复点云。优先使用 Global ID (第4列) 进行查重。
%
% 输入：
%   building_data : 包含 pointcloud_coordinates (Nx3 or Nx4) 的结构体数组
%
% 输出：
%   duplicate_info : 结构体
%       .coordinates          : 重复点的完整数据 (Nx3 or Nx4)
%       .occurrence_count     : 每个点出现的次数 (Nx1)
%       .building_indices     : (Cell) 每个重复点涉及的建筑物索引列表
%       .occurrence_locations : (Cell) 重复点在聚合列表中的位置索引

    fprintf('\n>>> 正在分析点云重复情况 (Duplicate Analysis)... \n');

    %% 1. 数据收集
    all_points = [];
    all_indices = [];
    
    num_bldgs = length(building_data);
    
    % 预分配内存 (估算)以加速，或者直接使用拼接
    % 这里为了代码简洁保持拼接，但逻辑上支持 Nx4
    for k = 1:num_bldgs
        pts = building_data(k).pointcloud_coordinates;
        if ~isempty(pts)
            all_points = [all_points; pts]; 
            all_indices = [all_indices; repmat(k, size(pts, 1), 1)];
        end
    end
    
    if isempty(all_points)
        fprintf('  没有点云数据，跳过分析。\n');
        duplicate_info = struct('coordinates',[], 'occurrence_count',[], ...
                                'building_indices',{}, 'occurrence_locations',{});
        return;
    end

    %% 2. 查重核心逻辑 (智能切换)
    
    % 判断是否包含 ID 列 (第4列)
    has_id = size(all_points, 2) >= 4;
    
    if has_id
        % === 方案 A: 基于 ID 查重 (极快且准确) ===
        % 只取第 4 列进行比对
        check_data = all_points(:, 4);
        [unique_vals, ~, ic] = unique(check_data);
        
        % 对应的完整点云数据 (取 unique 对应的行)
        % 注意：unique 的第二个返回值是 index，指向 unique_vals 在原数组中第一次出现的位置
        % 但我们需要重构 unique_pts 包含坐标以便输出
        
        % 方法：我们已经有了 unique_vals (ID)，我们需要找到它对应的坐标
        % 因为 ID 相同的点坐标肯定相同，所以随便取一个即可。
        % 利用 unique 的 index 特性：
        [~, first_occurrence_idx, ic] = unique(check_data);
        unique_pts = all_points(first_occurrence_idx, :);
        
    else
        % === 方案 B: 基于坐标查重 (兼容旧数据) ===
        tolerance = 1e-8;
        try
            [unique_pts, ~, ic] = uniquetol(all_points, tolerance, 'ByRows', true, 'DataScale', 1);
        catch
            rounded = round(all_points / tolerance) * tolerance;
            [unique_pts, ~, ic] = unique(rounded, 'rows');
        end
    end
    
    % 统计出现次数 (使用 histcounts 加速)
    % ic 是 1 到 length(unique_pts) 的索引
    counts = histcounts(ic, 1:(size(unique_pts,1)+1))';
    
    % 找出重复的点 (次数 > 1)
    is_dup_mask = counts > 1;
    
    dup_indices = find(is_dup_mask); % 这里是 unique 列表中的索引
    dup_counts = counts(is_dup_mask); % 重复次数
    
    %% 3. 构建详细输出
    num_unique_dups = length(dup_indices);
    
    % 提取重复点的完整信息 (Nx4)
    duplicate_info.coordinates = unique_pts(dup_indices, :);
    duplicate_info.occurrence_count = dup_counts;
    duplicate_info.building_indices = cell(num_unique_dups, 1);
    duplicate_info.occurrence_locations = cell(num_unique_dups, 1);
    
    if num_unique_dups > 0
        fprintf('  正在追踪重复点的来源建筑物...\n');
        for i = 1:num_unique_dups
            % 在 unique 列表中的索引
            target_unique_idx = dup_indices(i);
            
            % 在原始聚合列表(all_indices)中的位置
            locs = find(ic == target_unique_idx);
            
            % 获取对应的建筑物索引
            bldg_idxs = all_indices(locs);
            
            duplicate_info.building_indices{i} = unique(bldg_idxs)';
            duplicate_info.occurrence_locations{i} = locs;
        end
        
        % 按重复次数降序排列
        [~, sort_idx] = sort(duplicate_info.occurrence_count, 'descend');
        duplicate_info.coordinates = duplicate_info.coordinates(sort_idx, :);
        duplicate_info.occurrence_count = duplicate_info.occurrence_count(sort_idx);
        duplicate_info.building_indices = duplicate_info.building_indices(sort_idx);
        duplicate_info.occurrence_locations = duplicate_info.occurrence_locations(sort_idx);
    end
    
    %% 4. 输出报告
    total_dup_occurrences = sum(dup_counts);
    
    fprintf('----------------------------------\n');
    fprintf('  总点数:              %d\n', size(all_points, 1));
    fprintf('  重复点(唯一ID数):     %d\n', num_unique_dups);
    fprintf('  重复总次数(累加):     %d\n', total_dup_occurrences);
    
    if num_unique_dups == 0
        fprintf('  ✅ 完美！数据中不存在重复点。\n');
    else
        ratio = total_dup_occurrences / size(all_points, 1) * 100;
        fprintf('  ⚠️ 存在重复，占比: %.4f%%\n', ratio);
    end
    fprintf('----------------------------------\n');
end