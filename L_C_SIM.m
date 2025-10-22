% LC振荡电路交互式仿真 - 理想无损耗
% 通过输入 'toggle' 控制单刀双掷开关在A和B之间切换

clear all;
close all;
clc;

%% 参数设置
V = 10;              % 电源电压 (V)
C = 1e-6;            % 电容 (F) - 1 μF
L = 1e-3;            % 电感 (H) - 1 mH

% 计算振荡频率
omega0 = 1/sqrt(L*C);        % 角频率 (rad/s)
f0 = omega0/(2*pi);          % 频率 (Hz)
T0 = 1/f0;                   % 周期 (s)

%% 初始化仿真参数
dt = T0/1000;                % 时间步长（每个周期1000个点）
t_current = 0;               % 当前时间
vc = 0;                      % 电容电压初始值
i = 0;                       % 电流初始值
switch_pos = 'A';            % 开关初始位置：A=充电, B=振荡

% 历史数据记录 - 预分配空间
max_history_size = 1000000;  % 最大历史记录点数
t_history = zeros(1, max_history_size);
vc_history = zeros(1, max_history_size);
i_history = zeros(1, max_history_size);
history_index = 0;           % 当前历史记录索引

%% 创建图形窗口
fig = figure('Position', [100, 100, 1400, 800]);

fprintf('========== LC振荡电路交互式仿真 ==========\n');
fprintf('电源电压 V = %.2f V\n', V);
fprintf('电容 C = %.2f μF\n', C*1e6);
fprintf('电感 L = %.2f mH\n', L*1e3);
fprintf('固有频率 f0 = %.2f kHz\n', f0/1e3);
fprintf('振荡周期 T0 = %.2f μs\n', T0*1e6);
fprintf('==========================================\n\n');
fprintf('命令说明：\n');
fprintf('  sim <时间>      - 慢动作模拟指定时间 (单位: ms)，例如: sim 0.5\n');
fprintf('  fastsim <时间>  - 快速模拟指定时间 (单位: ms)，例如: fastsim 10\n');
fprintf('  toggle          - 切换开关位置 (A <-> B)\n');
fprintf('  x <时间>        - 设置横轴显示范围 (单位: ms)，例如: x 1\n');
fprintf('  status / s      - 显示当前状态\n');
fprintf('  quit / q        - 退出仿真\n\n');
fprintf('当前开关位置: %s\n', switch_pos);
fprintf('  (位置A: 电源给电容充电)\n\n');
fprintf('输入命令后按回车键...\n\n');

%% 主仿真循环
running = true;
step_count = 0;
steps_to_run = 50;  % 默认每次检查输入前运行50步（约0.01ms）
waiting_for_input = true;  % 初始等待用户输入
fast_mode = false;  % 是否为快速模式
x_range = 0.02;  % 横轴显示范围（ms），默认0.02ms

while running
    % 更新物理状态
    if switch_pos == 'A'
        % 位置A: 充电过程（使用小时间常数快速充电）
        tau = 0.5 * T0;  % 充电时间常数
        vc = vc + (V - vc) * dt / tau;
        i = (V - vc) / (tau/C);  % 充电电流
    else
        % 位置B: LC振荡（理想无损耗）
        % 使用辛欧拉法（Symplectic Euler），保持能量守恒
        % 先更新电流，再更新电压（或反之）
        % di/dt = -vc/L
        % dvc/dt = i/C
        
        % 方法1：辛欧拉法（半隐式）
        i = i - (vc / L) * dt;      % 先用旧的 vc 更新 i
        vc = vc + (i / C) * dt;     % 再用新的 i 更新 vc
        
        % 注：辛欧拉法能保持能量守恒，适合长时间模拟
    end
    
    % 更新时间和历史记录
    t_current = t_current + dt;
    history_index = history_index + 1;
    
    % 如果历史记录满了，扩展数组或采样保留
    if history_index > length(t_history)
        % 对历史数据进行降采样，保留一半的点
        keep_indices = 1:2:history_index-1;
        t_history(1:length(keep_indices)) = t_history(keep_indices);
        vc_history(1:length(keep_indices)) = vc_history(keep_indices);
        i_history(1:length(keep_indices)) = i_history(keep_indices);
        history_index = length(keep_indices) + 1;
    end
    
    t_history(history_index) = t_current;
    vc_history(history_index) = vc;
    i_history(history_index) = i;
    
    step_count = step_count + 1;
    
    % 定期更新图形（慢动作效果，快速模式下只在完成时更新）
    update_needed = false;
    if fast_mode
        % 快速模式：只在完成时更新一次
        if step_count >= steps_to_run
            update_needed = true;
        end
    else
        % 慢动作模式：每10步更新一次
        if mod(step_count, 10) == 0
            update_needed = true;
        end
    end
    
    if update_needed
        % 获取有效的历史数据
        valid_t = t_history(1:history_index);
        valid_vc = vc_history(1:history_index);
        valid_i = i_history(1:history_index);
        
        % 计算横轴范围（显示最后 x_range ms 的数据）
        if ~isempty(valid_t)
            xlim_max = valid_t(end) * 1e3;
            xlim_min = xlim_max - x_range;
        else
            xlim_min = 0;
            xlim_max = x_range;
        end
        
        % 绘图
        subplot(3,1,1);
        plot(valid_t*1e3, valid_vc, 'b-', 'LineWidth', 2);
        grid on;
        xlabel('时间 (ms)');
        ylabel('电容电压 v_C (V)');
        title(sprintf('电容电压 - 当前开关位置: %s  (时间: %.2f ms)', switch_pos, t_current*1e3));
        ylim([-V*1.2, V*1.2]);
        xlim([xlim_min, xlim_max]);
        
        subplot(3,1,2);
        plot(valid_t*1e3, valid_i*1e3, 'r-', 'LineWidth', 2);
        grid on;
        xlabel('时间 (ms)');
        ylabel('电流 i (mA)');
        xlim([xlim_min, xlim_max]);
        if switch_pos == 'A'
            title('充电电流');
        else
            title('振荡电流');
        end
        
        subplot(3,1,3);
        % 能量分布
        E_C = 0.5 * C * valid_vc.^2 * 1e6;  % 电容能量 (μJ)
        E_L = 0.5 * L * valid_i.^2 * 1e6;   % 电感能量 (μJ)
        E_total = E_C + E_L;
        plot(valid_t*1e3, E_C, 'b-', 'LineWidth', 2); hold on;
        plot(valid_t*1e3, E_L, 'r-', 'LineWidth', 2);
        plot(valid_t*1e3, E_total, 'k--', 'LineWidth', 1.5); hold off;
        grid on;
        xlabel('时间 (ms)');
        ylabel('能量 (μJ)');
        title('能量分布 (电场能 ↔ 磁场能)');
        legend('电场能 E_C', '磁场能 E_L', '总能量');
        xlim([xlim_min, xlim_max]);
        if max(E_total) > 0
            ylim([0, max(E_total)*1.2]);
        end
        
        drawnow;
        if ~fast_mode
            pause(0.001);  % 短暂暂停，产生慢动作效果（快速模式不暂停）
        end
    end
    
    % 定期检查用户输入或等待输入
    if waiting_for_input || step_count >= steps_to_run
        fprintf('> ');
        cmd = input('', 's');
        
        % 解析命令
        cmd_parts = strsplit(strtrim(cmd));
        cmd_name = '';
        if ~isempty(cmd_parts)
            cmd_name = lower(cmd_parts{1});
        end
        
        if strcmpi(cmd_name, 'sim')
            % sim 命令：慢动作模拟指定时间
            if length(cmd_parts) >= 2
                sim_time_ms = str2double(cmd_parts{2});
                if ~isnan(sim_time_ms) && sim_time_ms > 0
                    sim_time_s = sim_time_ms * 1e-3;  % 转换为秒
                    steps_to_run = round(sim_time_s / dt);
                    step_count = 0;
                    waiting_for_input = false;
                    fast_mode = false;
                    fprintf('\n开始慢动作模拟 %.3f ms (共 %d 步)...\n', sim_time_ms, steps_to_run);
                else
                    fprintf('错误：时间必须是正数\n');
                    fprintf('用法: sim <时间>  例如: sim 0.5\n\n');
                end
            else
                fprintf('错误：请指定模拟时间\n');
                fprintf('用法: sim <时间>  例如: sim 0.5\n\n');
            end
            
        elseif strcmpi(cmd_name, 'fastsim')
            % fastsim 命令：快速模拟指定时间（不慢放）
            if length(cmd_parts) >= 2
                sim_time_ms = str2double(cmd_parts{2});
                if ~isnan(sim_time_ms) && sim_time_ms > 0
                    sim_time_s = sim_time_ms * 1e-3;  % 转换为秒
                    steps_to_run = round(sim_time_s / dt);
                    step_count = 0;
                    waiting_for_input = false;
                    fast_mode = true;
                    fprintf('\n开始快速模拟 %.3f ms (共 %d 步)...\n', sim_time_ms, steps_to_run);
                    tic;  % 开始计时
                else
                    fprintf('错误：时间必须是正数\n');
                    fprintf('用法: fastsim <时间>  例如: fastsim 10\n\n');
                end
            else
                fprintf('错误：请指定模拟时间\n');
                fprintf('用法: fastsim <时间>  例如: fastsim 10\n\n');
            end
            
        elseif strcmpi(cmd_name, 'toggle')
            if switch_pos == 'A'
                switch_pos = 'B';
                fprintf('\n>>> 开关切换到位置 B (LC振荡回路) <<<\n');
                fprintf('电容开始放电，能量在电场和磁场之间振荡...\n');
                fprintf('当前电容电压: %.2f V\n', vc);
                fprintf('当前电流: %.2f mA\n\n', i*1e3);
            else
                switch_pos = 'A';
                fprintf('\n>>> 开关切换到位置 A (电源充电) <<<\n');
                fprintf('电源开始给电容充电...\n');
                fprintf('当前电容电压: %.2f V\n', vc);
                fprintf('当前电流: %.2f mA\n\n', i*1e3);
            end
            waiting_for_input = true;
            
        elseif strcmpi(cmd_name, 'x')
            % x 命令：设置横轴显示范围
            if length(cmd_parts) >= 2
                new_x_range = str2double(cmd_parts{2});
                if ~isnan(new_x_range) && new_x_range > 0
                    x_range = new_x_range;
                    fprintf('\n横轴显示范围已设置为: %.3f ms\n\n', x_range);
                else
                    fprintf('错误：范围必须是正数\n');
                    fprintf('用法: x <时间>  例如: x 1\n\n');
                end
            else
                fprintf('错误：请指定横轴显示范围\n');
                fprintf('用法: x <时间>  例如: x 1 (表示显示最近1ms的数据)\n');
                fprintf('当前横轴范围: %.3f ms\n\n', x_range);
            end
            waiting_for_input = true;
            
        elseif strcmpi(cmd_name, 'quit') || strcmpi(cmd_name, 'q') || strcmpi(cmd_name, 'exit')
            fprintf('\n退出仿真\n');
            running = false;
            
        elseif strcmpi(cmd_name, 'status') || strcmpi(cmd_name, 's')
            fprintf('\n当前状态:\n');
            fprintf('  开关位置: %s\n', switch_pos);
            fprintf('  仿真时间: %.3f ms\n', t_current*1e3);
            fprintf('  电容电压: %.3f V\n', vc);
            fprintf('  电流: %.3f mA\n', i*1e3);
            fprintf('  电场能: %.3f μJ\n', 0.5*C*vc^2*1e6);
            fprintf('  磁场能: %.3f μJ\n', 0.5*L*i^2*1e6);
            fprintf('  总能量: %.3f μJ\n\n', (0.5*C*vc^2 + 0.5*L*i^2)*1e6);
            fprintf('  横轴范围: %.3f ms\n\n', x_range);
            waiting_for_input = true;
            
        elseif isempty(cmd_name)
            % 空命令，什么都不做
            waiting_for_input = true;
            
        else
            fprintf('未知命令: %s\n', cmd);
            fprintf('可用命令:\n');
            fprintf('  sim <时间>      - 慢动作模拟 (ms)\n');
            fprintf('  fastsim <时间>  - 快速模拟 (ms)\n');
            fprintf('  toggle          - 切换开关\n');
            fprintf('  x <时间>        - 设置横轴范围 (ms)\n');
            fprintf('  status          - 显示状态\n');
            fprintf('  quit            - 退出\n\n');
            waiting_for_input = true;
        end
        
        % 如果完成了一轮模拟，重新等待输入
        if step_count >= steps_to_run && ~waiting_for_input
            if fast_mode
                elapsed = toc;
                fprintf('快速模拟完成！当前时间: %.3f ms (耗时: %.3f 秒)\n', t_current*1e3, elapsed);
            else
                fprintf('模拟完成！当前时间: %.3f ms\n', t_current*1e3);
            end
            waiting_for_input = true;
            step_count = 0;
            fast_mode = false;
        end
    end
end

fprintf('\n========== 仿真结束 ==========\n');
fprintf('总仿真时间: %.3f ms\n', t_current*1e3);
fprintf('最终电容电压: %.3f V\n', vc);
fprintf('最终电流: %.3f mA\n', i*1e3);
fprintf('==============================\n');
