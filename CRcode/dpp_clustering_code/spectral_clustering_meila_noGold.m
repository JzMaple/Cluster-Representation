function [clusters] = spectral_clustering_meila_noGold (matrix,k,init_flag, center_indexes)


  disp ('Transfering matrix to a similarity matrix W');

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


%  for i = 1:1:size_dd
 %  for j = 1:1:size_dd

  %  if (i < j)
   %   matrix(i,j) = matrix (i,j)/(q(i)*q(j));
    %else 
    % matrix(i,j) = matrix (i,j)/(q(j)*q(i));
    %end

   %end
  %end

%Notat that Lin used skewed divegence while I am using cosine similarity

%  matrix  = exp (matrix);


 disp ('creating P matrix');

  D = zeros(size_dd,size_dd);

  for i= 1:1:size_dd
    for j=1:1:size_dd
	    D(i,i) = D(i,i) + matrix(i,j);
    end
  end


P = inv(D)*matrix;

disp ('creating Y');

[v,d] = eig (P);  

d_vector = diag(d);

[d_sorted_val,d_sorted_loc] =sort(d_vector);

size_d_vector = max(size(d_vector));

for i=1:1:k-1

   relevant_eigs(i) = d_sorted_loc(size_d_vector -1 - i);

end 

X = v(:,relevant_eigs);

size_rows_x = size(X,1);

size_columns_x = size(X,2);

  for i=1:1:size_rows_x
    for j=1:1:size_columns_x  
       if (norm(X(i,:)) ~= 0)
	 Y(i,j) = X(i,j)/norm(X(i,:));
       end
   end
 end


disp ('running K-means')

   if (init_flag == 0)

     for i=1:1:100

	     min(Y)

	          clusters = kmeans (Y,k, 'distance','cosine','replicates',1,'start','sample','emptyaction','singleton');

%   clusters = kmeans (Y,k, 'distance','sqEuclidean','replicates',1,'start','sample','emptyaction','singleton');

% conf_mat = create_confusion_matrix (clusters',gold);

%  results_x(i) = calculate_greedy (conf_mat);

%  results_y(i) = calculate_greedy (conf_mat');

   end

  elseif (init_flag == 1)

  max_sumd = 0;

  for i=1:1:100

    init_matrix = Y(center_indexes(i,:),:);

[temp,C,sumd] = kmeans (Y,k, 'distance','sqEuclidean','replicates',1,'start',init_matrix,'emptyaction','singleton');
 
% conf_mat = create_confusion_matrix (temp,gold);

%results_x(i) = calculate_greedy (conf_mat);

%results_y(i) = calculate_greedy (conf_mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clusters_iter(:,i) = temp;

sum_sumd(i) = sum(sumd) ;


if (sum(sumd) > max_sumd)

  max_sumd = sum(sumd);

   index = i

end

  max_sumd_vec(i) = max_sumd;

  clusters = clusters_iter (:,i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  end

  end
