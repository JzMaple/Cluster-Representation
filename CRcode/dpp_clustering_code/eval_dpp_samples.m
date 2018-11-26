function [X,Y,samples, quality,diversity]  = eval_dpp_samples (matrix,trials, gold, scores_sort_loc, samples_matrix, k,eval_indexes)


  dd = diag (matrix);

  for i=1:1:trials
 
  
	  size(samples_matrix)

   samples(i,:) = samples_matrix (scores_sort_loc(length(scores_sort_loc) -i + 1),:)

    sample = samples(i,:);

    quality (i) = mean(dd(sample));

   diversity (i,:) = sum(matrix(:,sample)) - diag (matrix (sample,sample))';

   [clusters,results,alpha,max_vals] = knn_mapping (matrix, sample', 2, gold, eval_indexes);

   X(i,:) = results(:,1)';

   Y(i,:) = results(:,2)';


   %conf_mat = create_confusion_matrix (clusters, gold);

   %X(i) = calculate_greedy (conf_mat');

   %Y(i) = calculate_greedy (conf_mat);

  end
