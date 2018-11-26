function [clusters,results, cluster_elements] = random_iterative_clustering_by_mapping (matrix,k,m,number_of_iter,nus,threshold_bound,iteration2_threshold, num_of_repeats,gold_clusters)

%matrix - the DPP L-ensemple matrix 
%k - number of poiints in a subset
%m - number of samples in each iteration
%number_of_iter - number of iters for this protocol
%nus - number of used samples from DPP for the clustering process
%iteration2_threshold - inner paramater for clustering_by_mapping. For each threshold value, how many samples the algorithm tries to fit to the threshold before increasing the threshold
%num_of_repeats - how many times one runs the random baseline


  results = zeros (num_of_repeats,number_of_iter,3);

  cluster_elements = zeros (num_of_repeats,number_of_iter);

for r=1:1:num_of_repeats


  out_indices = [];

  size_matrix = max(size(matrix));

  for i=1:1:number_of_iter

      J = setdiff ([1:size_matrix],out_indices);

      matrix_J = matrix(J,J);

      size_matrix_J = max (size(matrix_J));

     if (min (size(matrix_J)) == 0)
 
       break;
  
     end

     if (size_matrix_J <= k)

       clusters (i,J) = floor(rand(1,size_matrix_J)*k) + 1;
        
       out_indices = [out_indices,J];

     else

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This  is for the dpp code

%      L = decompose_kernel (matrix_J);

%      [samples,stat,scores] = sample_dpp_mtimes(L,k,m); 

%      [sort_val,sort_loc] = sort (scores);

%      top_samples = samples (sort_loc (m - nus + 1:m),:);

%      top_scores = scores (sort_loc (m - nus + 1:m));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%       top_samples =  floor(rand (nus,k)*size_matrix_J) + 1

        top_samples = zeros(nus,k);

        for r=1:1:nus

	  flag = 0;

          while (flag == 0)

	    top_samples(r,:) = floor (rand(1,k)*size_matrix_J) + 1;

            if(max(size(unique(top_samples(r,:)))) == k)

	      flag = 1;

            end

          end            

        end

	max(top_samples)

      top_scores = rand (nus);

      [iter_clusters] = clustering_by_mapping (matrix_J, top_samples, top_scores, threshold_bound,iteration2_threshold);

      K = find (iter_clusters ~= 0);

      clusters (i,:) = zeros (1,size_matrix);

      clusters (i,J) = iter_clusters;

      out_indices = [out_indices,J(K)];

    end

      conf_mat_rand = create_confusion_matrix_zeros (clusters(i,:), gold_clusters);

     results (r, i, 1) = calculate_greedy (conf_mat_rand);

     results (r, i, 2) = calculate_greedy (conf_mat_rand');
 
    r = r

    i  = i

%     results (r, i, 3) = calculate_hungarian_munkers (conf_mat_rand');

      cluster_elements (r,i) = max(size(find (iter_clusters ~= 0)));

  end



end
