function [cluster_schemes] = print_hierarchical_clustering_no_evaluation (cell_matrix,  experiment_name)

  number_of_cluster_schemes = max(size(cell_matrix));

  number_of_points = size(cell_matrix{1},2);

  cluster_schemes = zeros (number_of_cluster_schemes,number_of_points);

  
  for i=1:1:number_of_cluster_schemes

    size_scheme = size(cell_matrix{i},1);

    total_clusters = 0;
 
      for j=1:1:size_scheme
 
	I = find (cell_matrix{i} (j,:) ~= 0);

        cluster_schemes (i,I) = cell_matrix{i} (j,I) + total_clusters;

        total_clusters = total_clusters + max (cell_matrix{i} (j,I));

     end

  end


  for i=1:1:number_of_cluster_schemes

	  %       conf_mat = create_confusion_matrix (cluster_schemes(i,:), gold);

	  % results(i) = calculate_greedy (conf_mat');

       file_name = sprintf ('clusters_%s_scheme%d', experiment_name,i)

       
       print_matrix_int (cluster_schemes(i,:), file_name);


 end
