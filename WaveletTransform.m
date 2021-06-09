clc;
clear;
close all;
img = imread("Label/Label2.JPG");
img = imgaussfilt(img, 3);
% 高斯滤波
J = 8;
% 尺度
[a, h, v, d] = haart2(img, J);
transformed = uint8(normalization(h{1}+v{1}+d{1}));
imwrite(transformed, "transformed.png");
% Haar小波变换
Mod = {};
% 模
Af = {};
% 幅角
Threshold = 4.5;
% 灰度阈值

for k = 1:J
    horizontal = [];
    horizontal = h{k};
    vertical = [];
    vertical = v{k};
    [height_h, width_h] = size(horizontal);
    [height_v, width_v] = size(vertical);
    for j = 1:min(height_h, height_v)
        for i = 1:min(width_h, width_v)
            Mod{k}(j, i) = sqrt(horizontal(j, i) .^ 2 + vertical(j, i) .^ 2);
            Af{k}(j, i) = atan(vertical(j, i) / horizontal(j, i));
        end
    end
    Mod{k}(Mod{k} <= Threshold) = 0;
end
% 计算Haar小波变换后每个尺度对应的模与幅角

M = [];
M = Mod{J};
Angel = [];
Angel = Af{J};
[height, width] = size(M);
T = zeros(height+2, width+2);
T(2:height+1, 2:width+1) = M;
for j = 2:height+1
    for i = 2:width+1
        % [-pi/8, pi/8)
        if(Angel(j-1, i-1) >= -pi/8 && Angel(j-1, i-1) < pi/8)
            if(T(j, i) >= T(j+1, i+1) && T(j, i) >= T(j, i+1) && T(j, i) >= T(j-1, i+1))
                M(j-1, i-1) = 1;
            else
                M(j-1, i-1) = 0;
            end
        % [pi/8, 3pi/8)
        elseif(Angel(j-1, i-1) >= pi/8 && Angel(j-1, i-1) < 3*pi/8)
            if(T(j, i) >= T(j, i+1) && T(j, i) >= T(j-1, i+1) && T(j, i) >= T(j-1, i))
                M(j-1, i-1) = 1;
            else
                M(j-1, i-1) = 0;
            end
        % [3pi/8, 5pi/8)    
        elseif(Angel(j-1, i-1) >= 3*pi/8 && Angel(j-1, i-1) < 5*pi/8)
            if(T(j, i) >= T(j-1, i+1) && T(j, i) >= T(j-1, i) && T(j, i) >= T(j-1, i-1))
                M(j-1, i-1) = 1;
            else
                M(j-1, i-1) = 0;
            end
        % [5pi/8, 7pi/8)    
        elseif(Angel(j-1, i-1) >= 5*pi/8 && Angel(j-1, i-1) < 7*pi/8)
            if(T(j, i) >= T(j-1, i+1) && T(j, i) >= T(j-1, i) && T(j, i) >= T(j-1, i-1))
                M(j-1, i-1) = 1;
            else
                M(j-1, i-1) = 0;
            end
        % [7pi/8, pi) [-pi, -7pi/8)
        elseif((Angel(j-1, i-1) >= 7*pi/8 && Angel(j-1, i-1) < pi)||(Angel(j-1, i-1) >= -pi && Angel(j-1, i-1) < -7*pi/8))
            if(T(j, i) >= T(j-1, i-1) && T(j, i) >= T(j, i-1) && T(j, i) >= T(j+1, i-1))
                M(j-1, i-1) = 1;
            else
                M(j-1, i-1) = 0;
            end
        % [-3pi/8, -pi/8)
        elseif(Angel(j-1, i-1) >= -3*pi/8 && Angel(j-1, i-1) < -pi/8)
            if(T(j, i) >= T(j, i+1) && T(j, i) >= T(j+1, i+1) && T(j, i) >= T(j+1, i))
                M(j-1, i-1) = 1;
            else
                M(j-1, i-1) = 0;
            end
         % [-5pi/8, -3pi/8)
         elseif(Angel(j-1, i-1) >= -5*pi/8 && Angel(j-1, i-1) < -3*pi/8)
            if(T(j, i) >= T(j+1, i-1) && T(j, i) >= T(j+1, i) && T(j, i) >= T(j+1, i+1))
                M(j-1, i-1) = 1;
            else
                M(j-1, i-1) = 0;
            end
        % [-7pi/8, -5pi/8]
        elseif(Angel(j-1, i-1) >= -7*pi/8 && Angel(j-1, i-1) < -5*pi/8)
            if(T(j, i) >= T(j, i-1) && T(j, i) >= T(j+1, i-1) && T(j, i) >= T(j+1, i))
                M(j-1, i-1) = 1;
            else
                M(j-1, i-1) = 0;
            end
        end
    end
end
last_layer = M;
%计算最后一个尺度的边缘点图像

for k = J-1:-1:1
    M = [];
    M = Mod{k};
    Angel = [];
    Angel = Af{k};
    [height, width] = size(M);
    T = zeros(height+2, width+2);
    T(2:height+1, 2:width+1) = M;
    % 填充一圈省了越界以及图像越来越小
    for j = 1:size(last_layer, 1)
        for i = 1:size(last_layer, 2)
            % 根据上一个尺度的边缘点图像的点出发，如果这个点是1，就将坐标×2，找到这个尺度下
            % 这个点的坐标
            if(last_layer(j ,i) == 1)
                for this_layer_x = 2*j-1:2*j+1
                    % 找到对应点后搜寻周围的3×3区域
                    if(this_layer_x > size(Angel, 1))
                        break;
                    end
                    % 如果越界了就退出
                    for this_layer_y = 2*i-1:2*i+1
                        if(this_layer_y > size(Angel, 2))
                            break;
                        end
                        % 越界了就退出
                        % [-pi/8, pi/8)
                        if(Angel(this_layer_x, this_layer_y) >= -pi/8 && Angel(this_layer_x, this_layer_y) < pi/8)
                            if(T(this_layer_x+1, this_layer_y+1) >= T(this_layer_x+2, this_layer_y+2) && T(this_layer_x+1, this_layer_y+1) >= T(this_layer_x+1, this_layer_y+2) && T(this_layer_x+1, this_layer_y+1) >= T(this_layer_x, this_layer_y+2))
                                M(this_layer_x, this_layer_y) = 1;
                            else
                                M(this_layer_x, this_layer_y) = 0;
                            end
                            % [pi/8, 3pi/8)
                        elseif(Angel(this_layer_x, this_layer_y) >= pi/8 && Angel(this_layer_x, this_layer_y) < 3*pi/8)
                            if(T(this_layer_x+1, this_layer_y+1) >= T(this_layer_x+1, this_layer_y+2) && T(this_layer_x+1, this_layer_y+1) >= T(this_layer_x, this_layer_y+2) && T(this_layer_x+1, this_layer_y+1) >= T(this_layer_x, this_layer_y+1))
                                M(this_layer_x, this_layer_y) = 1;
                            else
                                M(this_layer_x, this_layer_y) = 0;
                            end
                            % [3pi/8, 5pi/8)
                        elseif(Angel(this_layer_x, this_layer_y) >= 3*pi/8 && Angel(this_layer_x, this_layer_y) < 5*pi/8)
                            if(T(this_layer_x+1, this_layer_y+1) >= T(this_layer_x, this_layer_y+2) && T(this_layer_x+1, this_layer_y+1) >= T(this_layer_x, this_layer_y+1) && T(this_layer_x+1, this_layer_y+1) >= T(this_layer_x, this_layer_y))
                                M(this_layer_x, this_layer_y) = 1;
                            else
                                M(this_layer_x, this_layer_y) = 0;
                            end
                            % [5pi/8, 7pi/8)
                        elseif(Angel(this_layer_x, this_layer_y) >= 5*pi/8 && Angel(this_layer_x, this_layer_y) < 7*pi/8)
                            if(T(this_layer_+1, this_layer_y+1) >= T(this_layer_x, this_layer_y+2) && T(this_layer_+1, this_layer_y+1) >= T(this_layer_x, this_layer_y+1) && T(this_layer_+1, this_layer_y+1) >= T(this_layer_x, this_layer_y))
                                M(this_layer_x, this_layer_y) = 1;
                            else
                                M(this_layer_x, this_layer_y) = 0;
                            end
                            % [7pi/8, pi) [-pi, -7pi/8)
                        elseif((Angel(this_layer_x, this_layer_y) >= 7*pi/8 && Angel(this_layer_x, this_layer_y) < pi)||(Angel(this_layer_x, this_layer_y) >= -pi && Angel(this_layer_x, this_layer_y) < -7*pi/8))
                            if(T(this_layer_+1, this_layer_y+1) >= T(this_layer_x, this_layer_y) && T(this_layer_+1, this_layer_y+1) >= T(this_layer_x+1, this_layer_y) && T(this_layer_x+1, this_layer_y+1) >= T(this_layer_x+2, this_layer_y))
                                M(this_layer_x, this_layer_y) = 1;
                            else
                                M(this_layer_x, this_layer_y) = 0;
                            end
                            % [-3pi/8, -pi/8)
                        elseif(Angel(this_layer_x, this_layer_y) >= -3*pi/8 && Angel(this_layer_x, this_layer_y) < -pi/8)
                            if(T(this_layer_x+1, this_layer_y+1) >= T(this_layer_x+1, this_layer_y+2) && T(this_layer_x+1, this_layer_y+1) >= T(this_layer_x+2, this_layer_y+2) && T(this_layer_x+1, this_layer_y+1) >= T(this_layer_x+2, this_layer_y+1))
                                M(this_layer_x, this_layer_y) = 1;
                            else
                                M(this_layer_x, this_layer_y) = 0;
                            end
                            % [-5pi/8, -3pi/8)
                        elseif(Angel(this_layer_x, this_layer_y) >= -5*pi/8 && Angel(this_layer_x, this_layer_y) < -3*pi/8)
                            if(T(this_layer_x+1, this_layer_y+1) >= T(this_layer_x+2, this_layer_y) && T(this_layer_x+1, this_layer_y+1) >= T(this_layer_x+2, this_layer_y+1) && T(this_layer_x+1, this_layer_y+1) >= T(this_layer_x+2, this_layer_y+2))
                                M(this_layer_x, this_layer_y) = 1;
                            else
                                M(this_layer_x, this_layer_y) = 0;
                            end
                            % [-7pi/8, -5pi/8]
                        elseif(Angel(this_layer_x, this_layer_y) >= -7*pi/8 && Angel(this_layer_x, this_layer_y) < -5*pi/8)
                            if(T(this_layer_+1, this_layer_y+1) >= T(this_layer_x+1, this_layer_y) && T(this_layer_+1, this_layer_y+1) >= T(this_layer_x+2, this_layer_y) && T(this_layer_+1, this_layer_y+1) >= T(this_layer_x+2, this_layer_y+1))
                                M(this_layer_x, this_layer_y) = 1;
                            else
                                M(this_layer_x, this_layer_y) = 0;
                            end
                        end
                    end
                end
            end
        end
    end
    last_layer = M;
    % 经过处理的M成为前一个尺度的参考
end
%从后到前计算边缘点图像，后一个尺度的边缘点图像会成为前一个边缘点图像的参考图像

res = uint8(normalization(last_layer));
figure;
imshow(res);
imwrite(res, "res.png");
             
function res = normalization(img)
    max_value = max(max(img));
    min_value = min(min(img));
    res = img;
    res = 255 * (res - min_value) / (max_value - min_value);
end
% 图像归一化