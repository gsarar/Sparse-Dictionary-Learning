% load('denoised_image_approximation.mat')
im=imread('im450.jpg');
im=rgb2gray(im);
sigma=30;
var=(sigma/255)^2;
%noised_image = imnoise(im,'gaussian',0,var);
noised_image = imnoise(im,'speckle',var);
lambda=[30/sigma,10,30,50];
%lambda=50;
for kk=1:length(lambda)
    %taken from denoised_image_approximation according to relevant
    %dictionary and sigma
    approximation=reshape(approximation,8,8,111969);
    %Creating weight and merged_data matrices-> merged_data is when all
    %patches are placed one by one
    merged_data=zeros(384,304);
    weight_mat=zeros(384,304);
    for i=1:384-7
        for j=1:304-7
            merged_data(i:i+7,j:j+7)=approximation(:,:,i*297-297+j)+merged_data(i:i+7,j:j+7);
            weight_mat(i:i+7,j:j+7)=ones(8,8)+weight_mat(i:i+7,j:j+7);
        end
    end
    result=(lambda(kk)*double(noised_image(1:384,1:end))+merged_data)./(lambda(kk)+weight_mat);
    [peaksnr(kk),snr(kk)] = psnr(uint8(result),im(1:384,1:end));
    [n_peaksnr(kk),n_snr(kk)] = psnr(noised_image,im(1:384,1:end));
end
%imshow(uint8(result))