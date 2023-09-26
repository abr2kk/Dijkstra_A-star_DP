%% 这是华为杯数学建模2019年f题第一问A*求解
%% 清空工作区
clc ,clear
%% 读取并处理数据
regulate = xlsread('excel1.xlsx');
% 将数据按x轴排序
regulate_ord = sortrows(regulate, 2);

% 分离水平和垂直校正点
regulate_x = regulate_ord(regulate_ord(:, 5) == 0, :);
regulate_y = regulate_ord(regulate_ord(:, 5) ~= 0, :); 

%% 初步处理数据，构建无向图
reg_size = size(regulate, 1);  % 构建的无向图大小（行列数）
S = ones(reg_size, reg_size) * 10e8;  % 初始值为无穷
% 初级约束，即两点之间距离大于了最大可校正误差时就认为两点不可到达
alpha2 = 15;
beta1 = 20;
theta = 30;
% 构建无向图
for i = 1:reg_size
    for j = i + 1:reg_size
        distance = norm(regulate_ord(i, 2:4) - regulate_ord(j, 2:4));  % 两点之间距离
        delta_error = distance*0.001;  % 两点之间航迹偏差
        if(j ~= reg_size) % 正常校正点，也就是非B点
            if(regulate_ord(i, 5) == 0)  % 水平校正
                if(delta_error < beta1)
                    S(i, j) = distance;
                end
            else  % 垂直校正点
                if(delta_error < alpha2)
                    S(i, j) = distance;
                end
            end
        else  % B点
            if(delta_error < theta)
                S(i, j) = distance;
            end
        end
%         S(j, i) = S(i, j);  % 无向图，i到j和i到i距离相等
    end
end

figure(1)
scatter3(regulate_x(:, 2), regulate_x(:, 3), regulate_x(:, 4), '*b');
hold on
scatter3(regulate_y(:, 2), regulate_y(:, 3), regulate_y(:, 4), '*g');
AB =  [regulate_ord(1, 2:4); regulate_ord(reg_size(1), 2:4)];
hold on
scatter3(AB(:, 1), AB(:, 2), AB(:, 3), '*r');
pause(5)
%% A*算法寻优
start = 1;  % 起点
dest = reg_size;  % 终点
tic
[dist,path] = A_starfun(S,start,dest,regulate_ord);  % Dijkstra算法求解
toc
path
dist
%% 画出轨迹
hold on
coor_line = zeros(size(path, 2), 3);
for i = 1:size(path, 2)
    coor_line(i, :) = regulate_ord(path(i), 2:4);
end
plot3(coor_line(:, 1), coor_line(:, 2), coor_line(:, 3), 'r');
scatter3(coor_line(:, 1), coor_line(:, 2), coor_line(:, 3), 'r');

%% A*算法嵌入版
function [dist, path] = A_starfun(S, start, dest, reg)
    %% 初始化变量 
    points_num = size(S, 1); % 路径点个数
    U = ones(1, points_num); % 待访问的点，1为待访问，0为已经访问
    U(start) = 0;   % 将起点置为已经访问
    path_all = zeros(1, points_num); % 路径表，表示当前点的上一个路径是啥
    Gfun = ones(1, points_num)*inf;
    Allfun = ones(1, points_num)*inf;
    Ffun = ones(1, points_num)*inf;
    error = zeros(points_num,2); % 每一个点的误差
    now_point = start;
    Ffun(start) = 0;
    %% 约束条件变量
    alpha1 = 25;
    alpha2 = 15;
    beta1 = 20;
    beta2 = 25;
    theta = 30;
    
    k = 1;
    while sum(U)
        %% 画出寻找过程
        disp(now_point)
        disp(Allfun(now_point))
        hold on
        scatter3(reg(now_point, 2), reg(now_point, 3), reg(now_point, 4), 'oblack');
        pause(0.05)
        %% 更新A*算法代价
         for i = 1:points_num
             if U(i) == 1  % 该点未被遍历过
                now_error_y = error(now_point, 1) + S(now_point, i) * 0.001;
                now_error_z = error(now_point, 2) + S(now_point, i) * 0.001;
                if i ~= dest
                    if reg(i, 5) == 0
                        if now_error_y < beta2 && now_error_z < beta1
                            F = Ffun(now_point) + S(now_point, i);
                            G = norm(reg(i, 2:4) - reg(dest, 2:4));
                            All = k*G + F;
                            if All < Allfun(i)  % 该点代价更小，更新该点代价
                                 Allfun(i) = All;
                                 Ffun(i) = F;
                                 Gfun(i) = G;
                                 path_all(i) = now_point;
                                 error(i, 1) = 0;
                                 error(i, 2) = now_error_z;
                            end
                        end
                    else
                        if now_error_y < alpha2 && now_error_z < alpha1
                            F = Ffun(now_point) + S(now_point, i);
                            G = norm(reg(i, 2:4) - reg(dest, 2:4));
                            All = k*G + F;
                            if All < Allfun(i)  % 该点代价更小，更新该点代价
                                 Allfun(i) = All;
                                 Ffun(i) = F;
                                 Gfun(i) = G;
                                 path_all(i) = now_point;
                                 error(i, 1) = now_error_y;
                                 error(i, 2) = 0;
                            end
                        end
                    end
                else
                    if now_error_y < theta && now_error_z < theta
                        F = Ffun(now_point) + S(now_point, i);
                        G = norm(reg(i, 2:4) - reg(dest, 2:4));
                        All = k*G + F;
                        if All < Allfun(i)  % 该点代价更小，更新该点代价
                             Allfun(i) = All;
                             Ffun(i) = F;
                             Gfun(i) = G;
                             path_all(i) = now_point;
                        end
                    end
                end
             end
         end
        %% 选出最小代价的点
        used_All = Allfun; % copy一份D用于操作
        used_All(U == 0) = inf; % 将遍历过的点与起点距离设置为inf，这样就不会被选择
        [~, min_idx] = min(used_All); % 找出距离起点最近的点
        %% 更新当前点
        now_point = min_idx;
        U(now_point) = 0;
        if now_point == dest
            break;
        end
    end
    %% 结果处理
    dist = Ffun(dest);
    path = [dest];
    while start ~= dest
        dest = path_all(dest);
        path = [path, dest];
    end
    path = fliplr(path);
end


