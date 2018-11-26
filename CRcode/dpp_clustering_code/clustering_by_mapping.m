function [clusters] = clustering_by_mapping (matrix,samples,scores,threshold_bound, iteration2_thresh)

%  threshold_bound = 2;

  [val,loc] = sort (scores);

  size_matrix = max (size(matrix));

  clusters = zeros (1,size_matrix);
  
  number_of_samples = size(samples,1);

  number_of_clusters = size(samples,2);

  centers = zeros (number_of_clusters,size_matrix);

  flag = 1;

  first_iteration = 1;

  iterations = 0;

  iteration2 = 0;

  threshold = 0;

  while  (flag == 1)
 

    if (size(find(clusters == 0),2) == 0)
      flag = 0;
      continue;
    end

   
    if (number_of_samples == 0)
      break;
    end 

    not_assigned = size(find(clusters == 0))

    if (max(size(find(clusters == 0))) == 1)

      iterations = iterations + 1;

     if (iterations == 50)

%       return;

      flag = 0;

     end;

    end


    points = samples (loc(number_of_samples),:);


    if (first_iteration == 1)

      first_iteration = 0;

      cluster_centers = matrix(points,:);

      size(cluster_centers)

      for j=1:1:number_of_clusters

	clusters (points(j)) = j;

     end

    number_of_samples = number_of_samples - 1;

    continue;

   end


   already_clusteres = find (clusters > 0);

  XX =  size(intersect(points,already_clusteres))

  threshold = threshold

  if (size(intersect(points,already_clusteres),2) >= threshold)

    if (threshold == number_of_clusters)

    number_of_samples = number_of_samples - 1;
    continue;
    
   end 
  
    iteration2 = iteration2  + 1

    threhsold = threshold
  
      threshold_bound = threshold_bound

   if ( (iteration2 >= iteration2_thresh) & (threshold < threshold_bound))

      threshold = threshold + 1
      iteration2 = 0;
  

%    if (iteration2 == 25000)
%      return;
%       flag = 1;
%         break;
    
    elseif (iteration2 < iteration2_thresh)
 
     number_of_samples = number_of_samples - 1;     
    % iteration2 = 0;

   elseif (threshold == threshold_bound)

    break;

   end


    continue;

 else

      iteration2 = 0;

 end




    number_of_samples = number_of_samples - 1;

   conf_mat = zeros(number_of_clusters,number_of_clusters);

   for j= 1:1:number_of_clusters
     for k = 1:1:number_of_clusters

     conf_mat(j,k) = cluster_centers(k,points(j));

    end
   end

   conf_mat

    size(conf_mat)


   [score,map,sum1] = calculate_hungarian_munkers(conf_mat);

   for j=1:1:number_of_clusters

      if (clusters(points(j)) ~= 0)
	continue;
      end


     clusters (points(j)) = map(j);

     t = map(j);
  
    cluster_centers (t,:) = cluster_centers (t,:) + matrix (points(j),:);

   end

end
   

%mapping the rest of the points


  finalize = 0

if (finalize == 1)

   I = find (clusters == 0);

   J = find (clusters ~= 0);

   matrix_greedy = matrix(I,J);
  
   [a,b] = max(matrix_greedy'); 

   clusters(I) = clusters(J(b));

end  

%=============================================

  finalize = 0

if (finalize == 1)

  for i=1:1:size_matrix


   if (clusters(i) ~= 0)
    continue;
   end

   distance = zeros (1,number_of_clusters);

   for j=1:1:number_of_clusters   
      
      J = find (clusters == j);

      distance(j) = mean(matrix(i,J));     

   end
   
   [temp,clusters(i)] = max(distance);

  end
end

%=============================================

