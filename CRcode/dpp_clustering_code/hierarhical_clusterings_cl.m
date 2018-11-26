function clusters  = hierarhical_clusterings_cl (dpp_output_clusters, matrix, number_of_clusters)


  number_of_output_sets = size(dpp_output_clusters,1);

  size_matrix = max(size(matrix));

  centers = zeros (number_of_output_sets, size_matrix);

  clusters = cell (number_of_output_sets - 1,1); %we start with dpp_output_clusters and every iteration there will be one less
  
  clusters {1} = dpp_output_clusters;


for j = 1:1:number_of_output_sets -1

  j = j

  centers = zeros (number_of_output_sets - j +1, size_matrix);

  current_number_of_output_sets = number_of_output_sets - j + 1;

% Mapping: First creating confusion matrix, then using it for mapping

 
  conf_mat = ones (current_number_of_output_sets,current_number_of_output_sets);

  for l=1:1:current_number_of_output_sets
    for k=1:1:current_number_of_output_sets

         

            if (l == k)
              conf_mat(l,k) = 0;
   	      continue;
            end

	    for m=1:1:size_matrix                  
	       for s=1:1:size_matrix
                  
		   if ( (clusters{j}(l,m) > 0) &  (clusters{j}(k,s) > 0) & (matrix(m,s) < conf_mat(l,k) ) & (matrix(m,s) > 0 ))


           	          conf_mat (l,k) = matrix(m,s);


                        

                 end

             end % end for s

          end % end for m

	  a = conf_mat (l,k);

        

    end % for l
 end % for k


 current_number_of_output_sets

 conf_mat
 
index1 = 0;

 index2 = 0;

 max_val = 0;

 for l=1:1:current_number_of_output_sets
    for k=1:1:current_number_of_output_sets

	if (conf_mat (l,k) > max_val)

	  max_val = conf_mat (l,k);
          index1 = l;
          index2 = k;

        end % for k

    end % for l
 end

 mappaing_element = index1

  mappaing_element = index2

%return

% creating the unified clsutering

% First: copy the clusters that are not unified

 clusters{j+1} = zeros (number_of_output_sets -j, size_matrix);

 iters = number_of_output_sets -j +1;
 
 inner_index = 1;

 for m=1:1:iters

    if ( (m ~= index1) & (m ~= index2))

      clusters{j+1}(inner_index,:) = clusters{j}(m,:);

       inner_index = inner_index + 1;

    end %if

 end % for m

%Second: unifiy the most similar clusters


 conf_mat = ones (number_of_clusters, number_of_clusters)*10000*(-1);

 for l=1:1:size_matrix
   for k=1:1:size_matrix
 
   if ( (clusters{j}(index1,l) > 0) & (clusters{j}(index2,m) > 0) )
 
          val1 = clusters{j}(index1,l);

          val2 = clusters{j}(index2,m);

          if (conf_mat (val1,val2) < (-1)*matrix(l,k))

	    conf_mat (val1,val2) = (-1)*matrix(l,k);

          end

   end

   end %for k
 end % for l

 sum(sum(conf_mat));

 sum1 = calculate_greedy (conf_mat);

 sum2 = calculate_greedy (conf_mat');%'

 if (sum1 > sum2)

   [val,I] = max(conf_mat);
 
 else

   [val,I] = max(conf_mat');%'

 end

 clusters_index1 =  clusters{j}(index1,:);

 clusters_index2 =  clusters{j}(index2,:);

 for m=1:1:size_matrix

     if ((clusters_index2(m) ~= 0) & (clusters_index1(m) ~= 0) )

       clusters_index1(m) = clusters_index1(m);

       clusters_index2(m) = clusters_index2(m);

       return;

     end

     if (sum1 > sum2)

       if (clusters_index2(m) == 0)
         continue;
       end

       clusters_index1(m) = I(clusters_index2(m));

     else

       if (clusters_index1(m) == 0)
         continue;
        end

       clusters_index2(m) = I(clusters_index1(m));

     end
  
 end % for m

   if (sum1 > sum2)

     clusters{j+1}(inner_index,:) = clusters_index1;

     else

     clusters{j+1}(inner_index,:) = clusters_index2;
  
   end

 end % for j

%end

  
