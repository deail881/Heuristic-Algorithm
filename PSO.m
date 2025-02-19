% 粒子群优化算法（PSO）示例
function [BestF,BestX,HisBestF]=PSO(nPop,MaxIt,Low,Up,Dim,fitness)
% 参数设置
numParticles = nPop;      % 粒子数量
numDimensions =  Dim;      % 维度
maxIterations = MaxIt;    % 最大迭代次数
w = 0.5;                % 惯性权重
c1 = 1.5;               % 自我认知权重
c2 = 1.5;               % 社会认知权重

% 适应度函数（目标函数）
fitness_function = fitness;

% 初始化粒子位置和速度
positions = initialize(nPop,Dim,Low,Up);% 在[low, up]范围内初始化位置
% positions(:,1:2) = randi([1,99],[numParticles, 2]);
velocities = rand(numParticles, numDimensions) * 20 - 10; % 初始化速度

% 初始化个体和全局最佳
personalBestPositions = positions;                % 个体最佳位置
% p = pen_count(positions)%动态惩罚因子
for i=1:numParticles
    personalBestScores(i,:) = fitness_function(personalBestPositions(i,:)); % 个体最佳得分
end
[globalBestScore, bestIdx] = min(personalBestScores); % 全局最佳得分
globalBestPosition = personalBestPositions(bestIdx, :); % 全局最佳位置
HisBestF = [];

% 迭代过程
for iter = 1:maxIterations
    %     if  iter>3&&HisBestF(iter-1)-HisBestF(iter-2)<HisBestF(iter-2)-HisBestF(iter-3)
    %     p = pen_count(positions)%动态惩罚因子
    %     end
    for i = 1:numParticles
        % 更新速度
        r1 = rand(size(positions(i, :))); % 随机数
        r2 = rand(size(positions(i, :))); % 随机数
        velocities(i, :) = w * velocities(i, :) ...
            + c1 * r1 .* (personalBestPositions(i, :) - positions(i, :)) ...
            + c2 * r2 .* (globalBestPosition - positions(i, :));

        % 更新位置
        positions(i, :) = positions(i, :) + velocities(i, :);

        % 确保位置在一定范围内（边界处理）
        positions(i, :) = boundary( positions(i, :),Up,Low);
        %         positions(i, 1:2) = floor(positions(i, 1:2));  % 取整
        %         positions(i, 1:2) = max(positions(i, 1:2), 1); % 下界
        %         positions(i, 1:2) = min(positions(i, 1:2), 99);  % 上界


        % 评估当前粒子位置
        currentScore = fitness_function(positions(i, :));

        % 更新个体最佳
        if currentScore < personalBestScores(i)
            personalBestScores(i) = currentScore;
            personalBestPositions(i, :) = positions(i, :);
        end
    end
    % 更新全局最佳
    [currentGlobalBestScore, bestIdx] = min(personalBestScores);
    if currentGlobalBestScore < globalBestScore
        globalBestScore = currentGlobalBestScore;
        globalBestPosition = personalBestPositions(bestIdx, :);
    end
    BestX = globalBestPosition;
    BestF = globalBestScore;
    HisBestF = [HisBestF,globalBestScore];
    % 显示当前迭代的全局最佳得分
    %     disp(['Iteration ' num2str(iter) ': Global Best Score = ' num2str(globalBestScore)]);
end

% % 显示全局最佳位置
% disp(['Global Best Position: ' num2str(globalBestPosition)]);
% disp(['Global Best Score: ' num2str(globalBestScore)]);
end
% function x  = initialize(N,dim,lb,ub)
% le = size(ub,2);
% if le == 1
%     x = rand(N,dim)*(ub-lb)+lb;
% else
%     for i = 1:dim
%         high=ub(i);low=lb(i);
%         x(:,i)=rand(N,1).*(high-low)+low;
%     end
% end
% end
function x1  = boundary(x,ub,lb)
x1 = max(x,lb);
x1 = min(x1,ub);
end