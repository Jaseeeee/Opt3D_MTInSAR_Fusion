function [building_data_corrected, summary_stats] = correct_building_heights_by_ps(building_data, output_file)
% CORRECT_BUILDING_HEIGHTS_BY_PS 利用PS点云修正建筑物高度
%
% 核心逻辑：
%   1. 提取每个建筑物缓冲区内的点云高程 (Z)。
%   2. 过滤掉非物理的负值 (Z < 0)。
%   3. 计算 95% 分位数作为修正后的建筑物高度 (Robust Estimate)。
%   4. 如果某建筑物没有点云，则保留原始光学高度或标记为 NaN。
%
% 输入：
%   building_data : 包含 pointcloud_coordinates 的结构体数组
%   output_file   : (可选) 结果保存路径，如 'building_data_corrected.mat'
%
% 输出：
%   building_data_corrected : 更新后的结构体，新增了 corrected_height 等字段
%   summary_stats           : 统计信息结构体

    if nargin < 2, output_file = ''; end

    fprintf('\n=================================================\n');
    fprintf('>>> 开始利用 PS 点云修正建筑物高度 (95th Percentile)...\n');
    
    % 复制输入数据，以免修改原始变量
    building_data_corrected = building_data;
    
    num_bldgs = length(building_data);
    
    % 统计变量
    count_corrected = 0;
    count_skipped = 0;
    height_diffs = [];
    final_heights_vector = nan(num_bldgs, 1); % 保持索引对应的向量
    
    %% 核心处理循环
    for k = 1:num_bldgs
        b = building_data(k);
        
        % 1. 获取原始高度 (作为备份/对比)
        original_h = 0;
        if isfield(b, 'properties') && isfield(b.properties, 'height') && ~isempty(b.properties.height)
            original_h = double(b.properties.height);
        end
        building_data_corrected(k).original_height = original_h;
        
        % 2. 提取点云高度
        pts = b.pointcloud_coordinates;
        
        % 3. 判断是否有点云
        if ~isempty(pts) && size(pts, 1) > 0
            z_vals = pts(:, 3);
            
            % --- 数据清洗 ---
            % 去除高度小于0的点 (地表以下通常是噪声)
            z_vals = z_vals(z_vals > 0);
            
           if ~isempty(z_vals)
                % --- 核心算法：95分位数 ---
                ps_height = prctile(z_vals, 95);
                
                % 存入结构体
                building_data_corrected(k).corrected_height = ps_height;
                building_data_corrected(k).height_source = 'InSAR_PS_95';
                building_data_corrected(k).properties.height = ps_height;
                final_heights_vector(k) = ps_height;
                
                % 【核心修复】：只有当光学高度有效 (>0) 时，才计算统计差值！
                if original_h > 0
                    height_diffs(end+1) = ps_height - original_h; 
                end
                
                count_corrected = count_corrected + 1;
            else
                % 有点云但全是负值 -> 视为无效
                building_data_corrected(k).corrected_height = original_h;
                building_data_corrected(k).height_source = 'Optical_Original (Negative PS removed)';
                final_heights_vector(k) = original_h;
                count_skipped = count_skipped + 1;
            end
        else
            % 没有点云 -> 保留原值
            building_data_corrected(k).corrected_height = original_h;
            building_data_corrected(k).height_source = 'Optical_Original (No PS)';
            final_heights_vector(k) = original_h;
            count_skipped = count_skipped + 1;
        end
        
        if mod(k, 2000) == 0, fprintf('  已处理 %d / %d ...\n', k, num_bldgs); end
    end
    
    %% 汇总统计
    summary_stats.total_buildings = num_bldgs;
    summary_stats.corrected_count = count_corrected;
    summary_stats.skipped_count = count_skipped;
    summary_stats.heights_vector = final_heights_vector; % N x 1 向量
    
    if ~isempty(height_diffs)
        summary_stats.mean_diff = mean(height_diffs);
        summary_stats.std_diff = std(height_diffs);
        summary_stats.max_diff = max(abs(height_diffs));
    else
        summary_stats.mean_diff = 0;
    end
    
    %% 保存与输出
    if ~isempty(output_file)
        save(output_file, 'building_data_corrected', 'summary_stats');
        fprintf('  数据已保存至: %s\n', output_file);
    end
    
    fprintf('=== 高度修正完成 ===\n');
    fprintf('  总建筑物数:     %d\n', num_bldgs);
    fprintf('  成功修正数:     %d (%.1f%%)\n', count_corrected, count_corrected/num_bldgs*100);
    fprintf('  跳过/保留原值:  %d\n', count_skipped);
    if count_corrected > 0
        fprintf('  平均高度变化:   %.2f 米\n', mean(height_diffs));
    end
    fprintf('=================================================\n');
    
    % %% 简单的可视化对比 (可选)
    % if count_corrected > 0
    %     figure('Name', 'Height Correction Analysis', 'Color', 'w');
    %     histogram(height_diffs, 50);
    %     xlabel('Height Difference (InSAR - Optical) [m]');
    %     ylabel('Count');
    %     title('InSAR Correction Distribution');
    %     grid on;
    % end
end