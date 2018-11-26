function [X,Y,samples]  = eval_random (matrix,trials, gold, k, norm_flag,eval_indexes)

  size1 = size(matrix,1);

  for i=1:1:trials
  
    sample = 0;

    while (max(size(unique(sample))) ~= k)
       sample = floor (rand(1,k)*size1 ) + 1;
    end

   samples (i,:) = sample;

   [clusters,results,alphs,max_vals] = knn_mapping (matrix, sample', norm_flag, gold,eval_indexes);

   X(i,:) = results(:,1)';

   Y(i,:) = results(:,2)';

   %conf_mat = create_confusion_matrix (clusters, gold);

   %X(i) = calculate_greedy (conf_mat');

   %Y(i) = calculate_greedy (conf_mat);

  end
