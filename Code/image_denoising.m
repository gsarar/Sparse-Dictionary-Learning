clear
clc
im=imread('im450.jpg');
im=rgb2gray(im);
%sigma_values=[10 20 30];
sigma_values=30;
%[0,1]range std
for l=1:length(sigma_values)
    %noisy image creation
    sigma=sigma_values(l);
    var=(sigma/255)^2;
    noised_image = imnoise(im,'gaussian',0,var);
    %noised_image = imnoise(im,'speckle',var);
    % %imshow(noised_image)
    noised_data=zeros(8,8);
    for i=1:384-7
        for j=1:304-7
            noised_data=cat(3,noised_data,noised_image(i:i+7,j:j+7));
        end
    end
    noised_data=noised_data(:,:,2:end);
    noised_data=reshape(noised_data,1,64,111969);
    noised_data_matrix=zeros(64,111969);
    for i=1:64
        noised_data_matrix(i,:)=noised_data(1,i,:);
    end
    %DCT dictionary (just DCT dictionary is taken online)
    bb=8;
    eps=0.01;
    Pn=ceil(sqrt(441));
    DCT=zeros(bb,Pn);
    for kkk=0:1:Pn-1,
        V=cos([0:1:bb-1]'*kkk*pi/Pn);
        if kkk>0, V=V-mean(V); end;
        DCT(:,kkk+1)=V/norm(V);
    end;
    DCT=kron(DCT,DCT);
    dictionary=DCT;
    %
    %load('globalTrainedDictionary.mat');
    %dictionary=[ones(length(currDictionary(:,1)),1)/norm(ones(length(currDictionary(:,1)),1)) currDictionary];
    %load('learned_dictionary_oneimage_180.mat')
    %load('learned_dictionary_oneimage_80.mat')
    %load('good_working_learned_dictionary_oneimage_80.mat')
    
    representation_vectors=zeros(length(dictionary(1,:)),length(noised_data_matrix(1,:)));
    approximation=zeros(size(noised_data_matrix));
    error=1.15*sigma;
    for i=1:length(noised_data_matrix(1,:))
        if(mod(i,1000)==0)
            i
        end
        %[approximation_data,estimate_ideal_signal,index_set,residual,chosen_atoms]=omp(dictionary,training_data_matrix(:,i),sparsity_level);
        [approximation(:,i),estimate_ideal_signal,index_set,residual]=omp_error_bound(dictionary,noised_data_matrix(:,i),error);
        representation_vectors(:,i)=estimate_ideal_signal;
    end
    m = matfile('denoised_image_approximation.mat','Writable',true);
    eval(sprintf('m.approximation_dct%d=approximation;',sigma));
end