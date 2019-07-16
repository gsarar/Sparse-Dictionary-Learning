# Sparse-Dictionary-Learning
UCSD ECE269 Project

The two main algortihms OMP and K-SVD are implemented primarily in omp.m and ksvd.m. The other .m files, which has omp and ksvd in their names are slightly different versions of these two codes, so that I could use them for different parts. For example omp with an error bound or giving a sparsity level as smaller equal than T instead of equal T. That’s why I didn’t comment all these files but only these two main ones.

omp_main.m is the code where I check the performance of omp.m with the exact recovery condition in omp_check.m.

ksvd.m is the algorithm for synthetic data
image_ksvd.m is for dictionary learning from image, which is stored in im_data_matrix.mat
image_denoising is for patch_wise OMP implementation with omp_error_boud.m. The results of all dictionaries with different level of noising is in denoised_image_approximation.mat
patch_merging.m is for creating the single image from the sparse patches.

I excluded denoised_image_approximation.mat, due to the large file size.
