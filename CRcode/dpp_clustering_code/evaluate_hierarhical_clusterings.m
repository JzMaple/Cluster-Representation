function [results,avg_results,number_of_elements]  = evaluate_hierarhical_clusterings (h_clustering,gold)

  size_h = max(size(h_clustering));

  results = cell (1,size_h);

  number_of_elements = cell (1,size_h);

  avg_results  = cell (1,size_h);

 for i=1:1:size_h


    num_of_schemes = size(h_clustering{i},1);

    results{i} = zeros (num_of_schemes,4);

    avg_results{i} = zeros (1,4);

    number_of_elements{i} = zeros (1,num_of_schemes);

    total_elements = 0;

    for j=1:1:num_of_schemes

       conf_mat = create_confusion_matrix_zeros (h_clustering{i}(j,:),gold);

       number_of_elements{i}(j) = sum(sum(conf_mat)); 
   
       total_elements = total_elements +  sum(sum(conf_mat)); 
  
       results{i}(j,1) = calculate_greedy (conf_mat);

       avg_results{i}(1) = avg_results{i}(1) + results{i}(j,1)*number_of_elements{i}(j);

       results{i}(j,2) = calculate_greedy (conf_mat'); %'

       avg_results{i}(2)  = avg_results{i}(2) + results{i}(j,2)*number_of_elements{i}(j);

       results{i}(j,3) = calculate_hungarian_munkers  (conf_mat);

       avg_results{i}(3) = avg_results{i}(3) + results{i}(j,3)*number_of_elements{i}(j);

       results{i}(j,4) = calculate_v_measure (conf_mat);

       avg_results{i}(4) =  avg_results{i}(4) + results{i}(j,4)*number_of_elements{i}(j);
       

    end

    avg_results{i}(1) = avg_results{i}(1)/total_elements;

    avg_results{i}(2) = avg_results{i}(2)/total_elements;

    avg_results{i}(3) = avg_results{i}(3)/total_elements;

    avg_results{i}(4) = avg_results{i}(4)/total_elements;

 end
