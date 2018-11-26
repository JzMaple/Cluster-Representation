function [dpp_clusters, cl_clusters] = run_iterative_and_hierarchical (matrix, k, m, number_of_iter, nus, local_threshold_bound, iteration2_threshold)


  dpp_clusters = iterative_clustering_by_mapping (matrix, k, m, number_of_iter, nus, local_threshold_bound, iteration2_threshold);

cl_clusters = hierarhical_clusterings_cl (dpp_clusters, matrix, k);


