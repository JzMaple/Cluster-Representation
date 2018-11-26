function clusters = random_iterative_clustering_by_mapping_updated (matrix, k, m, number_of_iter, nus, local_threshold_bound, iteration2_threshold)

%matrix - the DPP L-ensemple matrix 
%k - number of poiints in a subset
%m - number of samples in each iteration
%number_of_iter - number of iters for this protocol
%nus - number of used samples from DPP for the clustering process
%iteration2_threshold - inner paramater for clustering_by_mapping. For each threshold value, how many samples the algorithm tries to fit to the threshold before increasing the threshold

  out_indices = [];

  size_matrix = max(size(matrix));


  for i=1:1:number_of_iter

	  fprintf (2, '===================================iteration = %d========================\n', i);

      J = setdiff ([1:size_matrix],out_indices);

      size_J = max(size(J));

      if ( (size_J <= k) & (i > 1))

        clusters(i,:) = zeros(1,size_matrix);

        clusters (i,J) = [1:size_J];

	break;

      end


      matrix_J = matrix(J,J);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      L = decompose_kernel (matrix_J);

     size_L_M = size(L.M)

    [samples,stat,scores] = sample_dpp_mtimes(L,k,m); 

      [sort_val,sort_loc] = sort (scores);

      top_samples = samples (sort_loc (m - nus + 1:m),:);

      top_scores = scores (sort_loc (m - nus + 1:m));


%       top_samples = floor (rand (nus,k)*size_J) + 1; 

%       top_scores = randperm(nus);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     flag = 1;

     threshold_bound = local_threshold_bound;

     while (flag == 1)
     

      [iter_clusters] = clustering_by_mapping (matrix_J, top_samples, top_scores, threshold_bound,iteration2_threshold);

      K = find (iter_clusters ~= 0);

     if ( ( max(size(K)) > k) | (threshold_bound == local_threshold_bound + 3))

	  flag = 0;

      else
     
	threshold_bound = threshold_bound + 1;

      end

     end

	   max_size_k = max(size(K))

      clusters (i,:) = zeros (1,size_matrix);

      clusters (i,J) = iter_clusters;

      out_indices = [out_indices,J(K)];

  end



  print_matrix_int (clusters, 'clusters_from_last_experiment');
