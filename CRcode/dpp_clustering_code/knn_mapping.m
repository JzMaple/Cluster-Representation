function [clusters,results,out_alpha,max_vals] = knn_mapping (matrix, cluster_center_indexes,norm_flag, gold_clusters,eval_indexes)

%%%%%%%%%%%%%%%%%%%%%%%%% computing q values


  
  dd = diag (matrix);

  size_dd = max (size (dd));

  q = zeros (1,size_dd);

  for i = 1:1:size_dd

    if (dd(i) < 0)

      dd(i)
      return;

    end

  q(i) = sqrt (dd(i));

  end


 if (norm_flag == 1)

  for i = 1:1:size_dd
   for j = 1:1:size_dd

    if (i < j)
      matrix(i,j) = matrix (i,j)/(q(i)*q(j));
    else 
     matrix(i,j) = matrix (i,j)/(q(j)*q(i));
    end

   end
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  size_matrix = max(size(matrix))
  number_of_clusters = max(size(cluster_center_indexes));

%  vals = zeros(size_matrix,number_of_clusters);

t = 0;

if (norm_flag <= 1)

  for i=1:1:size_matrix

     if (size(intersect(i,eval_indexes),1) == 0)
       continue;
     end

     t = t+1;

     if (size(intersect(i,cluster_center_indexes),1) == 1)

       i;
       
       for j=1:1:number_of_clusters

              %vals(i,j) = matrix(i,cluster_center_indexes(j));

	  if (cluster_center_indexes(j) == i)             
              clusters(t) = j;
        end

      end


      continue;

     end

     

       max_val = 0;
       index = 0;

       for j=1:1:number_of_clusters

            %vals(i,j) = matrix(i,cluster_center_indexes(j));

          if (matrix(i,cluster_center_indexes(j)) > max_val)        
	    max_val = matrix(i,cluster_center_indexes(j));
            index = j;
          end

       end

       clusters(t) = index;

       max_vals (t) = max_val;
     

  end % for i



   conf_mat = create_confusion_matrix (clusters, gold_clusters);
  
   results (1) = calculate_greedy (conf_mat); 

   results (2) = calculate_greedy (conf_mat'); 

   out_alpha = 1;

end %norm_flag <= 1


if (norm_flag == 2)



  k = 0;



for alpha = -2:0.25:5

  k = k + 1;

  t = 0;

  for i=1:1:size_matrix


     if (size(intersect(i,eval_indexes),1) == 0)
       continue;
     end

     t = t+1;



     if (size(intersect(i,cluster_center_indexes),1) == 1)

       i;
       
       for j=1:1:number_of_clusters

              %vals(i,j) = matrix(i,cluster_center_indexes(j));

	  if (cluster_center_indexes(j) == i)             
	    clusters(k,t) = j;
        end

      end


      continue;

     end

     

       max_val = 0;
       index = 0;

       for j=1:1:number_of_clusters

            %vals(i,j) = matrix(i,cluster_center_indexes(j));

          if (matrix(i,cluster_center_indexes(j))*(q(cluster_center_indexes(j))^alpha) > max_val)        
	    max_val = matrix(i,cluster_center_indexes(j))*(q(cluster_center_indexes(j))^alpha);
            index = j;
          end

       end

     if (index > 0)
       clusters(k,t) = index;

     elseif (index == 0)
        clusters(k,t) = floor (rand*number_of_clusters) + 1;
     
     end


       max_vals (k,t) = max_val;
     

    end %for i


  
   find (clusters(k,:) == 0)

   

   conf_mat = create_confusion_matrix (clusters(k,:), gold_clusters);

   
  

   results (k,1) = calculate_greedy (conf_mat); 

   results (k,2) = calculate_greedy (conf_mat'); 

   out_alpha (k) = alpha;

  end % for alpha

end % norm_flag == 2
