function [centers, point_center,avg_my_center,avg_other_center] = my_kmeans_cosine (data,k, init_start_flag, init_start,fixed_points_flag, fixed_points_numbers, fixed_points_matrix)

%k-means is implemented according to Gal Chechik's book, page 91, the
%distance metric is L1 norm
%samples are rows of data - there are p samples of size n
%k is the number of clusters
%init_start is an array of points that intialized the clusters. its size
%must be equal to k.
  
  
n = size (data, 2)
p = size (data, 1)

fixed_points_flag

disp (sprintf ('I am here1'));

%this version implements the option to add a matrix of fixed points (fixed_points_matrix) that contribute to  a cluster centroid and their assignment cannot be changed
%during the run of the algorithm. fixed_points_numbers should include the number of points associated with its class. The first fixed_points_numbers(1) raws in fixed_p%oints_matrix are associated with class 1, the next fixed_points_numbers(2) with class 2 and so on). 
% The option is used if fixed_points_flag == 1. Note that the points in fixed_points_matrix should not be included in the data matrix

if (fixed_points_flag == 1)

disp (sprintf ('fixed_points_flag == 1'));

if (sum(fixed_points_numbers) ~= size(fixed_points_matrix,1))

size_fixed_points_matrix = size(fixed_points_matrix)

sum_fixed_points_numbers = sum(fixed_points_numers);
 
return;

end


if (max(size(fixed_points_numbers)) ~= k) 

size_fixed_points_numbers = size(fixed_points_numbers)

k = k

return;

end

end



disp (sprintf ('I am here1'));

%step 1: initialization:

if ( (k ~= length (init_start)) & (init_start_flag == 1) )
%   fprintf (1,'%s\n','error, number of initial centers is different than
%   'k\n');
%    disp(sprintf('error, number of initial centers is different than
%    k'));
    
   k = k

   
     init_start = init_start ;
   return;
end


k;
n;
centers = rand (k,n);
point_center = ones (1,p);


%Initialzing the clusters

init_start_flag = init_start_flag

if (init_start_flag == 1)

for i=1:k

centers(i,:) = data (init_start(i),:);

end

end 

if (init_start_flag == 2)

  if ( (size(init_start, 1) ~= k) | (size(init_start,2) ~= n) )

    size_init_start_1 = size(init_start, 1) 

    size_init_start_2  = size(init_start, 2) 

      return;

 end


  centers = init_start;

end

%centers

%sum(centers,2)

disp (sprintf ('I am here2'));

% End of Code added to init the centers of the clusters from  the data points


threshold = 10^(-2);

flag = 1;
iter = 0;

while ( (flag > 0) & (iter < 200) )

    iter= iter +1
              
    flag   

    %step 2: clustering the samples to their nearest centers
    
    for i =1:p
      
        min_norm = 10000000;
        center = 1;
        
        for j=1:k
          
            temp1 = sum(data (i,:).*centers (j,:));
            temp2 = norm (data (i,:));
            temp3 =  norm (centers (j,:));
          
            if ( (temp2 == 0) | (temp3 == 0) )
              temp_norm = 0;
            else
             temp_norm = temp1/(temp2*temp3);
            end
            
            temp_norm = 1 - temp_norm;

            temp_norm = norm (data(i,:) - centers(j,:));
            
           if (temp_norm < min_norm)
               min_norm = temp_norm;
               center = j;
           end
        end      
        
        point_center (i) = center;
    end
    
    
    %step 3: calculating new centers
    
    
    new_centers = zeros (k,n);  
    for i=1:k

        counter = 0;

        for j=1:p
            if (point_center (j) == i)
                new_centers(i,:) = new_centers(i,:) + data(j,:);
                counter = counter + 1;
            end
        end

%fixed_points_flag, fixed_points_numbers, fixed_points_matrix


       if (fixed_points_flag == 1)
	 
	 if (i == 1)
	   fixed_idx1 = 1; 
	 else
           fixed_idx1 = sum(fixed_points_numbers(1:i-1));
        end

	 fixed_idx2 = fixed_points_numbers(i);


         for m = fixed_idx1+1:fixed_idx2

	   new_centers(i,:) = new_centers(i,:) + fixed_points_matrix(m,:);
           counter = counter + 1;

         end

       end

        if (counter > 0)            
           new_centers(i,:) = new_centers(i,:) / counter;
        end

    end
    
    %step 4: checking if there is any point that closer to its center than
    %in the former iteration

    flag = 0;
    avg_diff = 0;
    
    for i=1:p
       
        temp1 = sum(data(i,:).* centers (point_center(i),:));
        
        temp2 = norm (data (i,:));
        
        temp3 =  norm  (centers(point_center(i),:));
        
        
        
         if ( (temp2 == 0) | (temp3 == 0) )
              temp_norm = 0;
            else
             temp_norm = temp1/(temp2*temp3);
         end
        
      
        temp_norm = 1 - temp_norm;

        temp_norm = norm (data(i,:) - centers(point_center(i),:));
         
        former_dist  = temp_norm;
        
        
        temp1 = sum(data(i,:).* new_centers (point_center(i),:));
        
        temp2 = norm (data (i,:));

        temp3 =  norm  (new_centers (point_center(i),:));

         if ( (temp2 == 0) | (temp3 == 0) )
              temp_norm = 0;
         else
             temp_norm = temp1/(temp2*temp3);
         end
        
         temp_norm = 1 - temp_norm;

temp_norm = norm (data(i,:) - centers(point_center(i),:));

        current_dist = temp_norm;

        if (abs (former_dist - current_dist) > threshold)
            flag = flag + 1;
            abs (former_dist - current_dist);
            avg_diff = avg_diff + abs (former_dist - current_dist);
        end
    end
    
    flag;
    avg_diff/p;
    centers = new_centers;
        
end


%out_file = sprintf ('point_center_%d_%d_%d', k, p, ser_num)
%print_matrix(point_center, out_file);

 % step 5: Calculating the avarge distance (l1 norm) between each point and its center 
 %This step is for check and is not an integral part of the algorithm
 
 avg_my_center = 0;
 
 for i = 1:p
    
     temp1 = sum(data(i,:).* centers (point_center(i),:));
        
     temp2 = norm (data (i,:));
        
     temp3 =  norm  (centers(point_center(i),:));
        
     
     if ( (temp2 == 0) | (temp3 == 0) )
            temp_norm = 0;
         else
           temp_norm = temp1/(temp2*temp3);
       end
        
     
      temp_norm = 1 - temp_norm;

temp_norm = norm (data(i,:) - centers(point_center(i),:));

     avg_my_center = avg_my_center + temp_norm;
 end
 
 avg_my_center = avg_my_center/p;
 
 %Step 6: Calculating the avarage distance between each point and the other
 %centenrs.
 %This step is for check and is not an integral part of the algorithm
 
 avg_other_center = 0;
 for i = 1:p
     for j = 1:k
         if (point_center(i) ~= j)
           
             temp1 = sum(data(i,:).* centers (j,:));
        
             temp2 = norm (data (i,:));
        
             temp3 =  norm  (centers(j,:));
        
             
             if ( (temp2 == 0) | (temp3 == 0) )
                 temp_norm = 0;
             else
               temp_norm = temp1/(temp2*temp3);
             end
             
             temp_norm = 1 - temp_norm;

                        temp_norm = norm (data(i,:) - centers(j,:));
                                  
             avg_other_center = avg_other_center + temp_norm;
         end
     end
 end
 
 mechane = p*(k-1);
 avg_other_center = avg_other_center / mechane;
 
 %w_measure = zeros (p,p);
 
     

         
         
         
     
 
     
     
     
     
