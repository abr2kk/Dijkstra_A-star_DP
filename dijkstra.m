%% 这是华为杯2019年f题第一问Dijkstra求解
%% 清空工作区
clc,clear
%% 读取并处理数据
regulate = xlsread('excel1.xlsx');
% 将数据按x轴排序
regulate_ord = sortrows(regulate, 2);
% 分离水平和垂直校正点
regulate_x = regulate_ord(regulate_ord(:, 5) == 0, :);
regulate_y = regulate_ord(regulate_ord(:, 5) ~= 0, :); 

%% 初步处理数据，构建无向图
reg_size = size(regulate, 1);  % 构建的无向图大小（行列数）
S = ones(reg_size, reg_size) * inf;  % 初始值为无穷
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
        S(j, i) = S(i, j);  % 无向图，i到j和i到i距离相等
    end
end
%% Dijkstra算法寻优
start = 1;  % 起点
dest = reg_size;  % 终点
[dist,path,Distance] = dijkstrafun(S,start,dest,regulate_ord);  % Dijkstra算法求解
path
dist
%% 画路径图
figure(1)
scatter3(regulate_x(:, 2), regulate_x(:, 3), regulate_x(:, 4), '*b');
hold on
scatter3(regulate_y(:, 2), regulate_y(:, 3), regulate_y(:, 4), '*g');
AB =  [regulate_ord(1, 2:4); regulate_ord(reg_size(1), 2:4)];
hold on
scatter3(AB(:, 1), AB(:, 2), AB(:, 3), '*r');

coor_line = zeros(size(path, 2), 3);
for i = 1:size(path, 2)
    coor_line(i, :) = regulate_ord(path(i), 2:4);
end
plot3(coor_line(:, 1), coor_line(:, 2), coor_line(:, 3));
scatter3(coor_line(:, 1), coor_line(:, 2), coor_line(:, 3), 'r');

function [dist,path,D] = dijkstrafun(A,start,dest,reg)
    %% 算法变量
    points_num = size(A, 1); % 路径点个数
    U = ones(1, points_num); % 待访问的点，1为待访问，0为已经访问
    U(start) = 0;   % 将起点置为已经访问
    D = A(start, :);%*inf; % 每个点到起点的距离
%     D(start) = 0;  % 起点到起点的距离置0
    path_all = zeros(1, points_num); % 路径表，表示当前点的上一个路径是啥
    path_all(D ~= inf) = start; % 初始时所有距离非无穷的点的前一个节点为起点
    error = zeros(points_num,2); % 每一个点的误差
    % 约束条件变量
    alpha1 = 25;
    alpha2 = 15;
    beta1 = 20;
    beta2 = 25;
    theta = 30;
    % 误差初始化
    for j = 1:points_num
        if D(j) < 10e8
            if reg(j, 5) == 0
                error(j, 1) = 0;
                error(j, 2) = D(j) * 0.001;
            else
                error(j, 1) = D(j) * 0.001;
                error(j, 2) = 0;
            end
        end
    end
    chengfa = 0;
    min_idx = start;
    %% 遍历寻优
    while sum(U) % 遍历直到遍历完所有元素
        %% 取出当前路径中距离起点最近的点
        used_D = D; % copy一份D用于操作
        used_D(U == 0) = inf; % 将遍历过的点与起点距离设置为inf，这样就不会被选择
        [min_val, min_idx] = min(used_D); % 找出距离起点最近的点
        %% 更新与该点连接点到起点的距离
        for i = 1:points_num
            if U(i) == 1   % 确定该点未被遍历
                now_error_y = error(min_idx, 1) + A(min_idx, i) * 0.001;
                now_error_z = error(min_idx, 2) + A(min_idx, i) * 0.001;
                if i ~= dest
                    if reg(i, 5) == 0
                        if now_error_y < beta2 && now_error_z < beta1
                            if D(i)>(D(min_idx)+A(min_idx, i))
                                D(i)=(D(min_idx)+A(min_idx,i) + chengfa);
                                path_all(i)=min_idx;
                                error(i, 1) = 0;
                                error(i, 2) = now_error_z;
                            end
                        end
                    else
                        if now_error_y < alpha2 && now_error_z < alpha1
                            if D(i)>(D(min_idx)+A(min_idx, i))
                                D(i)=(D(min_idx)+A(min_idx,i) + chengfa); 
                                path_all(i)=min_idx;
                                error(i, 1) = now_error_y;
                                error(i, 2) = 0;
                            end
                        end
                    end
                else
                    if now_error_y < theta && now_error_z < theta
                        if D(i)>(D(min_idx)+A(min_idx, i))
                            D(i)=(D(min_idx)+A(min_idx,i) + chengfa); 
                            path_all(i)=min_idx;
                        end
                    end
                end
            end
        end
        U(min_idx) = 0;
    end
    %% 结果处理
    dist = D(start, dest);
    path = [dest];
    while start ~= dest
        dest = path_all(dest);
        path = [path, dest];
    end
    path = fliplr(path);

end