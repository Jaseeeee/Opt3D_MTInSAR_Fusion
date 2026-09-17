function [building_data_final, summary] = generate_building_intersection_data(data_circle, data_sar, output_file)
% GENERATE_BUILDING_INTERSECTION_DATA (适配 Nx4 索引修正版)
% 

    fprintf('\n=================================================\n');
    fprintf('>>> 开始计算点云与几何形状的交集 (Intersection Strategy)...\n');

    % 1. 基础校验
    if length(data_circle) ~= length(data_sar)
        error('两个输入数据的建筑物数量不一致！请检查是否使用了同一个 GeoJSON。');
    end
    
    num_bldgs = length(data_circle);
    building_data_final = data_sar; % 继承基础结构
    
    total_pts_circle = 0;
    total_pts_sar = 0;
    total_pts_final = 0;
    
    % 关闭几何警告
    warning('off', 'MATLAB:polyshape:repairedBySimplify');
    
    % 2. 遍历每个建筑物
    for k = 1:num_bldgs
        % === A. 点云交集 (使用 ID 优化) ===
        pts_A = data_circle(k).pointcloud_coordinates;
        pts_B = data_sar(k).pointcloud_coordinates;
        
        count_A = size(pts_A, 1);
        count_B = size(pts_B, 1);
        total_pts_circle = total_pts_circle + count_A;
        total_pts_sar = total_pts_sar + count_B;
        
        pts_final = [];
        
        if count_A > 0 && count_B > 0
            % --- 核心修改：优先使用 ID (第4列) 进行交集运算 ---
            % 只要两个来源都是 Nx4 且第4列是 ID，这样做最快且最安全
            if size(pts_A, 2) >= 4 && size(pts_B, 2) >= 4
                id_A = pts_A(:, 4);
                id_B = pts_B(:, 4);
                
                % 找共同 ID
                common_ids = intersect(id_A, id_B);
                
                % 提取这些 ID 对应的完整行 (从 pts_A 中提取即可，因为数据源一样)
                [~, idx_in_A] = ismember(common_ids, id_A);
                pts_final = pts_A(idx_in_A, :);
            else
                % 如果没有 ID，回退到基于坐标的整行匹配
                pts_final = intersect(pts_A, pts_B, 'rows');
            end
        end
        
        count_final = size(pts_final, 1);
        total_pts_final = total_pts_final + count_final;
        
        % --- 关键修复：同时更新坐标和索引 ---
        building_data_final(k).pointcloud_coordinates = pts_final;
        
        if count_final > 0 && size(pts_final, 2) >= 4
            % 必须更新 indices，否则它保留的是 data_sar 的旧索引
            building_data_final(k).pointcloud_indices = pts_final(:, 4);
        elseif count_final == 0
            building_data_final(k).pointcloud_indices = [];
        else
            % 如果是 Nx3 数据，没有 ID，就不存 indices 或者置空
            building_data_final(k).pointcloud_indices = []; 
        end
        
        building_data_final(k).pointcloud_count = count_final;
        
        % === B. [新增] 缓冲区几何形状交集 (保持不变) ===
        poly_A = data_circle(k).buffer_polygon; 
        poly_B = data_sar(k).buffer_polygon;    
        
        final_poly = [];
        final_area = 0;
        
        try
            if ~isempty(poly_A) && ~isempty(poly_B)
                ps_A = polyshape(poly_A(:,1), poly_A(:,2));
                ps_B = polyshape(poly_B(:,1), poly_B(:,2));
                ps_final = intersect(ps_A, ps_B);
                
                if ps_final.NumRegions > 0
                    final_poly = ps_final.Vertices;
                    final_area = area(ps_final);
                else
                    final_poly = [];
                end
            end
        catch
            final_poly = poly_B;
            if ~isempty(poly_B)
                final_area = polyarea(poly_B(:,1), poly_B(:,2));
            end
        end
        
        building_data_final(k).buffer_polygon = final_poly;
        building_data_final(k).buffer_area = final_area;
        building_data_final(k).buffer_method = 'geometry_intersection';
        
        % === C. 更新高程统计 ===
        if count_final > 0 && size(pts_final, 2) >= 3
            z_vals = pts_final(:,3);
            building_data_final(k).elevation_stats = struct(...
                'min', min(z_vals), 'max', max(z_vals), 'mean', mean(z_vals), 'std', std(z_vals));
        else
            building_data_final(k).elevation_stats = struct('min',NaN,'max',NaN,'mean',NaN,'std',NaN);
        end
        
        if mod(k, 2000) == 0
            fprintf('  已处理 %d / %d ...\n', k, num_bldgs);
        end
    end
    
    warning('on', 'MATLAB:polyshape:repairedBySimplify');
    
    % 3. 保存与汇总
    summary.total_buildings = num_bldgs;
    summary.points_circle = total_pts_circle;
    summary.points_sar = total_pts_sar;
    summary.points_intersection = total_pts_final;
    summary.reduction_rate = (total_pts_circle - total_pts_final) / total_pts_circle * 100;
    
    if nargin >= 3 && ~isempty(output_file)
        save(output_file, 'building_data_final', 'summary');
        fprintf('  结果已保存至: %s\n', output_file);
    end
    
    fprintf('=== 交集计算完成 ===\n');
    fprintf('  - 几何形状已更新为 [四周形 ∩ 定向]\n');
    fprintf('  - 最终点云总数: %d\n', total_pts_final);
end