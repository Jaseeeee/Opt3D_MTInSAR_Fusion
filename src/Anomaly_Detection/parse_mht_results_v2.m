function [pdef_type, pdef_desc] = parse_mht_results_v2(mht_vec, Ha_struct)

    N = length(mht_vec);
    pdef_type = zeros(N, 1);
    pdef_desc = cell(N, 1);
    
    for i = 1:N
        idx = mht_vec(i);
        
        if idx == 0
            pdef_type(i) = 0; 
            pdef_desc{i} = 'Linear';
            
        elseif idx == -1
            pdef_type(i) = -1;
            pdef_desc{i} = 'Noise/Undefined';
            
        else
            model_components = Ha_struct(idx).models; 
            
            has_H = any(contains(model_components, 'heaviside'));
            has_B = any(contains(model_components, 'breakpoint'));
            has_P = any(contains(model_components, 'periodic'));
            
            % --- 互斥的全组合判定 ---
            
            % 1. 单一成分
            if has_H && ~has_B && ~has_P
                pdef_type(i) = 1; 
                pdef_desc{i} = 'Heaviside Only';
                
            elseif ~has_H && has_B && ~has_P
                pdef_type(i) = 2; 
                pdef_desc{i} = 'Breakpoint Only';
                
            elseif ~has_H && ~has_B && has_P
                pdef_type(i) = 3; 
                pdef_desc{i} = 'Periodic Only';
                
            % 2. 两种成分组合
            elseif has_H && ~has_B && has_P
                pdef_type(i) = 4; % [新编号] 阶跃 + 周期
                pdef_desc{i} = 'Heaviside + Periodic';
                
            elseif ~has_H && has_B && has_P
                pdef_type(i) = 5; % [新编号] 断点 + 周期
                pdef_desc{i} = 'Breakpoint + Periodic';
                
            elseif has_H && has_B && ~has_P
                pdef_type(i) = 6; % [新编号] 阶跃 + 断点
                pdef_desc{i} = 'Heaviside + Breakpoint';
                
            % 3. 三种成分全有
            elseif has_H && has_B && has_P
                pdef_type(i) = 7; % [新编号] 复杂全组合
                pdef_desc{i} = 'Heaviside + Breakpoint + Periodic';
                
            else
                pdef_type(i) = 99;
                pdef_desc{i} = 'Other Ha';
            end
        end
    end
end