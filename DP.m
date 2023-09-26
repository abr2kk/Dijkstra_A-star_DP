%% 创建日期：2023.8.10
%% 创建人：yqh
%% 作用：动态规划求解A点到B点的最优路径,华为杯2019年f题第一题
%% 清空工作区
clc, clear;

%% 读取数据，并将数据排序分割
regulate = xlsread('excel1.xlsx');
% 将数据按x轴排序，然后分离水平和垂直校正点
regulate_ord = sortrows(regulate, 2);
regulate_x = regulate_ord(regulate_ord(:, 5) == 0, :);
regulate_y = regulate_ord(regulate_ord(:, 5) ~= 0, :);

%% 下面是动态规划的一些参数

reg_size = size(regulate);  % 校正点矩阵形状
delta = 0.001;  % 误差变化率
% 能够校正的误差阈值
alpha1 = 25;
alpha2 = 15;
beta1 = 20;
beta2 = 25;
theta = 30;
coor = zeros(reg_size(1), 3);   % 飞机经过校正点的实际位置，此处认为是校正点位置
coor(1, :) = [0, 50000, 5000];  % A点实际位置
error = zeros(reg_size(1), 2);  % 飞机经过校正点时的误差
all_distance = ones(1, reg_size(1)) * inf;  % A点到其他点的最短路径
all_distance(1) = 0;    % A点到A点最短路径为0
coor_map = ones(reg_size(1), reg_size(1)+1);   % 存储路径
chengfa = 100;    %惩罚项，调制很大时，校正点个数最少

%% 下面开始动态规划寻优
for i = 2:reg_size(1)-1   % 遍历每一个校正点
    for j = 1:i-1    % 计算出到第i个校正点的最短路径
        %% 第一问，最少校正点数的路径
        distance = norm(regulate_ord(i, 2:4) - coor(j, :));   % 计算第i个校正点和第j个实际点之间的距离，coor(j)为经过第j个点时飞机的实际坐标
        delta_y = distance * delta + error(j, 1);   % 计算第j个点到第i个点之间的水平误差，error(j, 1)为飞机经过第j个点校正后的水平误差
        delta_z = distance * delta + error(j, 2);   % 计算第j个点到第i个点之间的垂直误差，error(j, 2)为飞机经过第j个点校正后的垂直误差
        % 在满足有效性的前提下，与当前最优路径进行比较
        if(regulate_ord(i, 5) == 0) % 第i个点为水平校正点，y->水平
            if(delta_y < beta2 && delta_z < beta1)  % 第j个点可到达第i个点
                if(all_distance(i) > all_distance(j) + distance)    % 如果存在更小距离， 就重新赋值第i个节点的量
                    all_distance(i) = all_distance(j) + distance;   % 重新赋值最短距离
                    % coor(i, :) = [regulate_ord(i, 2), regulate_ord(i, 3) + delta_y, regulate_ord(i, 4) + delta_z]; % 重新赋值飞机在第i点的实际坐标
                    coor(i, :) = regulate_ord(i, 2:4);
                    error(i, 1) = 0;    % 水平矫正误差
                    error(i, 2) = delta_z;  % 垂直误差保留
                    coor_map(i, 1:coor_map(j, reg_size(1)+1)) = coor_map(j, 1:coor_map(j, reg_size(1)+1));
                    coor_map(i, reg_size(1)+1) = coor_map(j, reg_size(1)+1) + 1;
                    coor_map(i, coor_map(i, reg_size(1)+1)) = i;
                end
            end
        else  % 第i个点为垂直校正点
            if(delta_y < alpha2 && delta_z < alpha1)    % 第j个点可到达第i个点
                if(all_distance(i) > all_distance(j) + distance)    % 如果存在更小距离， 就重新赋值第i个节点的量
                    all_distance(i) = all_distance(j) + distance;   % 重新赋值最短距离
                    % coor(i, :) = [regulate_ord(i, 2), regulate_ord(i, 3) + delta_y, regulate_ord(i, 4) + delta_z]; % 重新赋值飞机在第i点的实际坐标
                    coor(i, :) = regulate_ord(i, 2:4);
                    error(i, 1) = delta_y;  % 水平误差保留
                    error(i, 2) = 0;    % 垂直矫正误差
                    coor_map(i, 1:coor_map(j, reg_size(1)+1)) = coor_map(j, 1:coor_map(j, reg_size(1)+1));
                    coor_map(i, reg_size(1)+1) = coor_map(j, reg_size(1)+1) + 1;
                    coor_map(i, coor_map(i, reg_size(1)+1)) = i;
                end
            end
        end
         %%  第三问，通过率100%情况下的路径
%         % 先判断该点的有效性，即是否可到达
%         distance = norm(regulate_ord(i, 2:4) - coor(j, :));   % 计算第i个校正点和第j个实际点之间的距离，coor(j)为经过第j个点时飞机的实际坐标
%         delta_y = distance * delta + error(j, 1);   % 计算第j个点到第i个点之间的水平误差，error(j, 1)为飞机经过第j个点校正后的水平误差
%         delta_z = distance * delta + error(j, 2);   % 计算第j个点到第i个点之间的垂直误差，error(j, 2)为飞机经过第j个点校正后的垂直误差
%         % 在满足有效性的前提下，与当前最优路径进行比较
%         if(regulate_ord(i, 5) == 0) % 第i个点为水平校正点，y->水平
%             if(delta_y < beta2 && delta_z < beta1)  % 第j个点可到达第i个点
%                 if(all_distance(i) > all_distance(j) + distance)    % 如果存在更小距离， 就重新赋值第i个节点的量
%                     all_distance(i) = all_distance(j) + distance;   % 重新赋值最短距离
%                     % coor(i, :) = [regulate_ord(i, 2), regulate_ord(i, 3) + delta_y, regulate_ord(i, 4) + delta_z]; % 重新赋值飞机在第i点的实际坐标
%                     coor(i, :) = regulate_ord(i, 2:4);
%                     if regulate_ord(i, 6) == 1
%                         error(i, 1) = min(5, delta_y);    % 水平矫正误差
%                         error(i, 2) = delta_z;  % 垂直误差保留
%                     else
%                         error(i, 1) = 0;    % 水平矫正误差
%                         error(i, 2) = delta_z;  % 垂直误差保留
%                     end
%                     coor_map(i, 1:coor_map(j, reg_size(1)+1)) = coor_map(j, 1:coor_map(j, reg_size(1)+1));
%                     coor_map(i, reg_size(1)+1) = coor_map(j, reg_size(1)+1) + 1;
%                     coor_map(i, coor_map(i, reg_size(1)+1)) = i;
%                 end
%             end
%         else  % 第i个点为垂直校正点
%             if(delta_y < alpha2 && delta_z < alpha1)    % 第j个点可到达第i个点
%                 if(all_distance(i) > all_distance(j) + distance)    % 如果存在更小距离， 就重新赋值第i个节点的量
%                     all_distance(i) = all_distance(j) + distance;   % 重新赋值最短距离
%                     % coor(i, :) = [regulate_ord(i, 2), regulate_ord(i, 3) + delta_y, regulate_ord(i, 4) + delta_z]; % 重新赋值飞机在第i点的实际坐标
%                     coor(i, :) = regulate_ord(i, 2:4);
%                     if regulate_ord(i, 6) == 1
%                         error(i, 1) = delta_y;    % 水平矫正误差
%                         error(i, 2) = min(5, delta_z);  % 垂直误差保留
%                     else
%                         error(i, 1) = delta_y;    % 水平矫正误差
%                         error(i, 2) = 0;  % 垂直误差保留
%                     end
%                     coor_map(i, 1:coor_map(j, reg_size(1)+1)) = coor_map(j, 1:coor_map(j, reg_size(1)+1));
%                     coor_map(i, reg_size(1)+1) = coor_map(j, reg_size(1)+1) + 1;
%                     coor_map(i, coor_map(i, reg_size(1)+1)) = i;
%                 end
%             end
%         end
%             % 比较结果若是最优，就把数据保存下来
    end
    all_distance(i) = all_distance(i) + chengfa; % 将第i个校正点路程数据的数据保存下来
end

%% 单独处理B点
for j = 1:reg_size(1)-1
    % 先判断该点的有效性，即是否可到达
    distance = norm(regulate_ord(reg_size(1), 2:4) - coor(j, :));   % 计算第i个校正点和第j个实际点之间的距离，coor(j)为经过第j个点时飞机的实际坐标
    delta_y = distance * delta + error(j, 1);   % 计算第j个点到第i个点之间的水平误差，error(j, 1)为飞机经过第j个点校正后的水平误差
    delta_z = distance * delta + error(j, 2);   % 计算第j个点到第i个点之间的垂直误差，error(j, 2)为飞机经过第j个点校正后的垂直误差
    % 在满足有效性的前提下，与当前最优路径进行比较
    if(delta_y < theta && delta_z < theta)  % 第j个点可到达第i个点
        if(all_distance(reg_size(1)) > all_distance(j) + distance)    % 如果存在更小距离， 就重新赋值第i个节点的量
            all_distance(reg_size(1)) = all_distance(j) + distance;   % 重新赋值最短距离
            coor_map(reg_size(1), 1:coor_map(j, reg_size(1)+1)) = coor_map(j, 1:coor_map(j, reg_size(1)+1));
            coor_map(reg_size(1), reg_size(1)+1) = coor_map(j, reg_size(1)+1) + 1;
            coor_map(reg_size(1), coor_map(reg_size(1), reg_size(1)+1)) = reg_size(1);
        end
    end
end

%% 显示最短路径数据
points_num = coor_map(reg_size(1), reg_size(1)+1)-2
length = all_distance(reg_size(1)) - points_num * chengfa
points = coor_map(reg_size(1), 1:coor_map(reg_size(1), reg_size(1)+1));

%% 画图
figure(2)
scatter3(regulate_x(:, 2), regulate_x(:, 3), regulate_x(:, 4), '*b'); % 画水平校正点
hold on
scatter3(regulate_y(:, 2), regulate_y(:, 3), regulate_y(:, 4), '*g'); % 画垂直校正点
AB =  [regulate_ord(1, 2:4); regulate_ord(reg_size(1), 2:4)];
hold on
scatter3(AB(:, 1), AB(:, 2), AB(:, 3), '*r'); % 画AB两点

coor_line = zeros(coor_map(reg_size(1), reg_size(1)+1), 3);
for i = 1:coor_map(reg_size(1), reg_size(1)+1)
    coor_line(i, :) = regulate_ord(coor_map(reg_size(1), i), 2:4);
end
plot3(coor_line(:, 1), coor_line(:, 2), coor_line(:, 3)); % 画出路径
scatter3(coor_line(:, 1), coor_line(:, 2), coor_line(:, 3), 'r'); % 画出路径点

%% 将与原坐标点对应起来
coor_org = zeros(1, coor_map(reg_size(1), reg_size(1)+1)-2);
for i = 1:coor_map(reg_size(1), reg_size(1)+1)-2
    for j = 2:reg_size(1)-1
        if(coor_line(i+1, :) == regulate(j, 2:4))
            coor_org(i) = regulate(j, 1);
            break;
        end
    end
end
coor_org
