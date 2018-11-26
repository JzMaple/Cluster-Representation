function [std_out,out_scores] = plot_qstd_vs_score (matrix, samples,scores, gold_clusters,eval_indexes)


  size_samples = max(size(samples));
  
  [a,b] = sort (scores);

  dd_matrix = diag(matrix);


   for i = 1:1:100


      std_out(i) = std(dd_matrix(samples(b(size_samples-i+1),:)));


     cluster_center_indexes = samples (b(size_samples-i+1),:);

     norm_flag = 0;

     [clusters,results,out_alpha,max_vals] = knn_mapping (matrix, cluster_center_indexes',norm_flag, gold_clusters,eval_indexes);


      


       out_scores(i,:) = results;


   end
