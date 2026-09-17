function [ts_matrix, global_ids] = load_ts_using_pt_pars(data_path, date_vec)
    % 1. 基础路径检查
    if ~exist(data_path, 'dir')
        error('路径不存在: %s', data_path);
    end
    
    % 2. 确定干涉对数量 nifg
      nifg = length(date_vec) - 1; 
    
    % 3. 初始化 Pt_Pars 对象
    % 【关键修改】我们将 np 设为 0。
    % 你的第一个脚本就是这样做的 (P0=Pt_Pars(0,1,nifg,3))。
    % 这样 Pt_read_files 会读取 pt1 文件，自动计算出真正的点数 np。
    P_obj = Pt_Pars(0, 1, nifg, 3);
    
    % 4. 读取文件
    % 指向 levelALL 文件夹
    level_folder = fullfile(data_path, 'levelALL', filesep);
    
    if ~exist(level_folder, 'dir')
        % 尝试不带 filesep 的情况，或者检查是否就在 data_path 下
        if exist(fullfile(data_path, 'pt1'), 'file')
            level_folder = [data_path, filesep];
        else
            error('找不到 levelALL 文件夹或 pt1 文件，请检查路径结构: %s', data_path);
        end
    end
    
    try
        % 调用 Pt_Pars 的读取方法
        P_obj = P_obj.Pt_read_files(level_folder);
    catch ME
        fprintf('读取出错！可能是 Data_Read 函数未在路径中。\n');
        rethrow(ME);
    end
    
    % 5. 提取时间序列 (ts)并转置
    % def_disp_all 也就是你的 pdef_disp_all1
    % 原始维度通常是 [Time x Points] 或 [Points x Time]，取决于 Pt_Pars 实现
    % 根据你的 mnt_test 需求，我们需要 [Points x Time]
    
    raw_ts = P_obj.def_disp_all; 
    
    % 获取真正的点数 np (读取后会自动更新)
    np = size(P_obj.pt, 1);
    
    % 维度检查与转置
    [dim1, dim2] = size(raw_ts);
    N_dates = length(date_vec);
    
    if dim1 == np
        % 已经是 [Points x Time]
        ts_matrix = raw_ts;
    elseif dim2 == np
        % 是 [Time x Points]，需要转置
        ts_matrix = raw_ts';
    else
        warning('时间序列矩阵维度 (%dx%d) 与点数 (%d) 不匹配，请检查 Pt_Pars。', dim1, dim2, np);
        ts_matrix = raw_ts'; % 盲猜转置
    end
    
    % 6. 二次检查时间维度
    if size(ts_matrix, 2) ~= N_dates
        % 如果列数不对 (例如少1列)，补0
        if size(ts_matrix, 2) == N_dates - 1
            ts_matrix = [zeros(np, 1), ts_matrix];
        else
             warning('警告：时间序列列数 (%d) 与日期数 (%d) 仍不一致。', size(ts_matrix, 2), N_dates);
        end
    end
    
    % 7. 生成 Global ID
    global_ids = (1:np)';
    
    fprintf('  成功加载: %d 个点, %d 个历元\n', np, size(ts_matrix, 2));
end