function print_matrix(matrix,filename)

  fid = fopen (filename,'w');

  size1 = size(matrix,1);

  size2 = size(matrix,2);

  
  size1
 
  size2

  for(i=1:1:size1)

    for(j=1:1:size2)

   % i

  %  j

%      matrix(i,j)

      fprintf (fid, '%f    ', matrix(i,j));


      

   end

  fprintf (fid, '\n');

  end


  fclose (fid);
