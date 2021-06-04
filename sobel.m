clc;
clear;
img = imread("Label/Label3.JPG");
img = imgaussfilt(img, 3);
imwrite(img, "gaussian.png");
Sobel_Threshold = 60;
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
Sobel_res = G;
for j = 1:height
    for i = 1:width
        if(Sobel_res(j, i) >= Sobel_Threshold)
            Sobel_res(j, i) = 255;
        else
            Sobel_res(j, i) = 0;
        end
    end
end
imshow(Sobel_res, []);
imwrite(Sobel_res, "Label_3.png");