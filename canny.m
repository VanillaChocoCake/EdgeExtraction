clc;
clear;
img = imread("Label/Label2.JPG");
img = imgaussfilt(img, 3);
Sobel_x = [-1,0,1;-2,0,2;-1,0,1];
Sobel_y = [-1,-2,-1;0,0,0;1,2,1];
[height, width] = size(img);
extended_img = zeros(height+2, width+2);
extended_img(2:height+1, 2:width+1) = img(:, :);
G_x = zeros(height, width);
G_y = zeros(height, width);
for i = 2:width + 1
    for j = 2 :height + 1
        mat = [extended_img(j-1, i-1) extended_img(j-1, i) extended_img(j-1, i+1);
            extended_img(j, i-1) extended_img(j, i) extended_img(j, i+1);
            extended_img(j+1, i-1) extended_img(j+1, i) extended_img(j+1, i+1)];
        G_x(j-1, i-1) = sum(sum(Sobel_x .* mat));
        G_y(j-1, i-1) = sum(sum(Sobel_y .* mat));
    end
end
G = sqrt(G_x .* G_x + G_y .* G_y);
G_MAX = max(max(G));
G_MIN = min(min(G));
G = (G - G_MIN) / (G_MAX - G_MIN);
imwrite(G, "grad.png");
grad_area = zeros(height, width);
% 区域1为(3Pi/4, Pi] && [-Pi/4 ,0)
% 区域2为(Pi/2, 3Pi/4] && [-Pi/2, -Pi/4)
% 区域3为(Pi/4, Pi/2] && [-3Pi/4, -Pi/2)
% 区域4位(0, Pi/2] && [-Pi, -3Pi/4)
% 不在任一一个区域的为区域0
for j = 1:height
    for i = 1:width
        grad_x = G_x(j, i);
        grad_y = G_y(j, i);
        if(grad_x == 0 && grad_y == 0)
            grad_area(j, i) = 0;
            continue;
        end
        if(grad_x <= 0 && grad_y >= 0 && grad_y < -grad_x || grad_x >=0 && grad_y < 0 && grad_y >= -grad_x)
            grad_area(j, i) = 1;
        elseif(grad_x <= 0 && grad_y >= 0 && grad_y >= -grad_x || grad_x >=0 && grad_y <= 0 && grad_y < -grad_x)
            grad_area(j, i) = 2;
        elseif(grad_x >= 0 && grad_y >= 0 && grad_y > grad_x || grad_x <=0 && grad_y <= 0 && grad_y <= grad_x)
            grad_area(j, i) = 3;
        elseif(grad_x >= 0 && grad_y >= 0 && grad_y <= grad_x || grad_x <=0 && grad_y <= 0 && grad_y > grad_x)
            grad_area(j, i) = 4;
        end
    end
end
grad_up = zeros(height, width);
grad_down = zeros(height, width);
filled_G = zeros(height+2, width+2);
filled_G(2:height+1, 2:width+1) = G;
% 计算梯度
for j = 2:height+1
    for i = 2:width+1
        t = abs(G_x(j-1, i-1) ./ G_y(j-1, i-1));
        if(grad_area(j-1, i-1) == 1)
            grad_up(j-1, i-1) = filled_G(j, i+1) .* (1-t) + filled_G(j-1, i+1) .* t;
            grad_down(j-1, i-1) = filled_G(j, i-1) .* (1-t) + filled_G(j+1, i-1) .* t;
        elseif(grad_area(j-1, i-1) == 2)
            grad_up(j-1, i-1) = filled_G(j-1, i) .* (1-t) + filled_G(j-1 ,i+1) .* t;
            grad_down(j-1, i-1) = filled_G(j+1, i) .* (1-t) + filled_G(j+1, i-1) .* t;
        elseif(grad_area(j-1, i-1) == 3)
            grad_up(j-1, i-1) = filled_G(j-1, i) .* (1-t) + filled_G(j-1 ,i-1) .* t;
            grad_down(j-1, i-1) = filled_G(j+1, i) .* (1-t) + filled_G(j+1, i+1) .* t;
        elseif(grad_area(j-1, i-1) == 4)
            grad_up(j-1, i-1) = filled_G(j, i-1) .* (1-t) + filled_G(j-1 ,i-1) .* t;
            grad_down(j-1, i-1) = filled_G(j, i+1) .* (1-t) + filled_G(j+1, i+1) .* t;
        end
    end
end
G_MAX = zeros(height, width);
% 判断是否为梯度方向极大值
for j = 1:height
    for i = 1:width
        if(G(j, i) >= grad_up(j, i) && G(j, i) >= grad_down(j, i))
            G_MAX(j, i) = G(j, i);
        end
    end
end
res = zeros(height, width);
Thershold_Strong = 8e-2;
Thershold_Weak = 7e-2;
% 强边缘与弱边缘分辨
for j = 2:height-1
    for i = 2:width-1
        if(G_MAX(j, i) >= Thershold_Strong)
            res(j, i) = 1;
        elseif(G_MAX(j, i) < Thershold_Weak)
            res(j, i) = 0;
        else
            num = 0;% 梯度极大值的相邻像素点个数
            if(G_MAX(j-1, i-1) ~= 0)
                num = num + 1;
            end
            if(G_MAX(j-1, i) ~= 0)
                num = num + 1;
            end
            if(G_MAX(j-1, i+1) ~= 0)
                num = num + 1;
            end
            if(G_MAX(j, i-1) ~= 0)
                num = num + 1;
            end
            if(G_MAX(j, i+1) ~= 0)
                num = num + 1;
            end
            if(G_MAX(j+1, i-1) ~= 0)
                num = num + 1;
            end
            if(G_MAX(j+1, i) ~= 0)
                num = num + 1;
            end
            if(G_MAX(j+1, i+1) ~= 0)
                num = num + 1;
            end
            if num >= 4
                res(j, i) = 1;
            end
        end
    end
end
imshow(res, []);
imwrite(res, "canny_Label3.png");

            


        
            

