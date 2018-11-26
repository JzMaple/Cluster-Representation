img = imread('img_1229.png');
[m,n,k]=size(img);
img1 = img(:,1:n/2,:);
img2 = img(:,n/2+1:end,:);

h1_1 = zeros(256,1); h1_2 = zeros(256,1); h1_3 = zeros(256,1);
h2_1 = zeros(256,1); h2_2 = zeros(256,1); h2_3 = zeros(256,1);
ch1_1 = img1(:,:,1); ch1_2 = img1(:,:,2); ch1_3 = img1(:,:,3);
ch2_1 = img2(:,:,1); ch2_2 = img2(:,:,2); ch2_3 = img2(:,:,3);
tmp1_1 = ch1_1(:); tmp1_2 = ch1_2(:); tmp1_3 = ch1_3(:);
tmp2_1 = ch2_1(:); tmp2_2 = ch2_2(:); tmp2_3 = ch2_3(:);

for i=1:length(tmp1_1)
    h1_1(tmp1_1(i)+1) = h1_1(tmp1_1(i)+1)+1;
    h1_2(tmp1_2(i)+1) = h1_2(tmp1_2(i)+1)+1;
    h1_3(tmp1_3(i)+1) = h1_3(tmp1_3(i)+1)+1;
end

figure, hold on; title('Histogram of source image','FontSize',12);
h_r = plot([0:255], h1_1,'Color','r','LineWidth',1.5);
h_g = plot([0:255], h1_2,'Color','g','LineWidth',1.5);
h_b = plot([0:255], h1_3,'Color','b','LineWidth',1.5);

L=legend([h_r,h_g,h_b],'Red','Green','Blue');
L.FontSize = 12;
ylabel('\fontname{times new roman}Accumulated counts','FontSize',15);
xlabel('\fontname{times new roman}Pixel intensity','FontSize',15);
axis tight;
hold off;

for i=1:length(tmp2_1)
    h2_1(tmp2_1(i)+1) = h2_1(tmp2_1(i)+1)+1;
    h2_2(tmp2_2(i)+1) = h2_2(tmp2_2(i)+1)+1;
    h2_3(tmp2_3(i)+1) = h2_3(tmp2_3(i)+1)+1;
end

figure, hold on; title('Histogram of generated image','FontSize',12);
h_r = plot([0:255], h2_1,'Color','r','LineWidth',1.5);
h_g = plot([0:255], h2_2,'Color','g','LineWidth',1.5);
h_b = plot([0:255], h2_3,'Color','b','LineWidth',1.5);

L=legend([h_r,h_g,h_b],'Red','Green','Blue');
L.FontSize = 12;
ylabel('\fontname{times new roman}Accumulated counts','FontSize',15);
xlabel('\fontname{times new roman}Pixel intensity','FontSize',15);
axis tight;
hold off;

diff_image = double(img1)-double(img2);
ch_dif_1 = diff_image(:,:,1); ch_dif_2 = diff_image(:,:,2); ch_dif_3 = diff_image(:,:,3);
tmp_dif_1=ch_dif_1(:); tmp_dif_2=ch_dif_2(:); tmp_dif_3=ch_dif_3(:); 
edges = [-100:1:100];
h1 = zeros(201,1); h2 = zeros(201,1); h3 = zeros(201,1);

for i=1:length(tmp_dif_1)
    if tmp_dif_1(i)>=-100 && tmp_dif_1(i)<=100
        h1(tmp_dif_1(i)+101)=h1(tmp_dif_1(i)+101)+1;
    end
    if tmp_dif_2(i)>=-100 && tmp_dif_2(i)<=100
        h2(tmp_dif_2(i)+101)=h2(tmp_dif_2(i)+101)+1;
    end
    if tmp_dif_3(i)>=-100 && tmp_dif_3(i)<=100
        h3(tmp_dif_3(i)+101)=h3(tmp_dif_3(i)+101)+1;
    end
end
figure, hold on; title('Histogram of the difference image','FontSize',12);
h_r = plot([-100:1:100],h1,'Color','r','LineWidth',1.5);
h_g = plot([-100:1:100],h2,'Color','g','LineWidth',1.5);
h_b = plot([-100:1:100],h3,'Color','b','LineWidth',1.5);
L=legend([h_r,h_g,h_b],'Red','Green','Blue');
L.FontSize = 12;
ylabel('\fontname{times new roman}Accumulated counts','FontSize',15);
xlabel('\fontname{times new roman}Pixel difference intensity','FontSize',15);
axis tight;
hold off;

img_diff = double(img1)-double(img2);
figure;
imshow(uint8(abs(img_diff)));
