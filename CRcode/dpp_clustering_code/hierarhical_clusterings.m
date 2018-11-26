function clusters  = hierarhical_clusterings (dpp_output_clusters, matrix, number_of_clusters)


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

  for i = 1:1:current_number_of_output_sets

      T = find (clusters{j}(i,:) ~= 0);

    centers(i,:) = mean(matrix(T,:));

  end % for i


  conf_mat = zeros (current_number_of_output_sets,current_number_of_output_sets);

  for l=1:1:current_number_of_output_sets
    for k=1:1:current_number_of_output_sets

            if (l == k)
   	      continue;
            end

	    conf_mat (l,k) = sum( centers(l,:).*centers(k,:));

%/(norm(centers(l,:))*norm(centers(k,:)));

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

 cluster_center_index1 = zeros (number_of_clusters, size_matrix);

 cluster_center_index2 = zeros (number_of_clusters, size_matrix);

 for m=1:1:number_of_clusters

   temp = size(find (clusters{j}(index1,:)  == m),2);
   
   if (temp > 0)
     I1 = find (clusters{j}(index1,:)  == m);
     cluster_center_index1(m,:) = mean(matrix(I1,:));
   end


   temp = size(find (clusters{j}(index2,:)  == m),2);
   
   if (temp > 0)
     I2 = find (clusters{j}(index2,:)  == m);
     cluster_center_index2(m,:) = mean(matrix(I2,:));
   end

 end% for m

 conf_mat = zeros (number_of_clusters, number_of_clusters);

 for l=1:1:current_number_of_output_sets
   for k=1:1:current_number_of_output_sets     
 
	   conf_mat(l,k) = sum(cluster_center_index1(l,:).*cluster_center_index2(k,:))./(norm(cluster_center_index1(l,:))*norm(cluster_center_index2(k,:)));

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

       clusters_index1(m) = clusters_index1(m)

       clusters_index2(m) = clusters_index2(m)

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

  
