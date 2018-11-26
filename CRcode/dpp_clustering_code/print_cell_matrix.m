function print_cell_matrix (cell_matrix, full_cl_clusters,file_name)

%we assume this is a line vector

fid = fopen(file_name,'w');

size1 = size(cell_matrix,2);
size2 = size(cell_matrix{1},2);

   for i=1:1:size1

     fprintf (fid, 'avg_results{%d}     ',i);
   
     num = size(unique(full_cl_clusters(i,:)));

     for j=1:1:size2

	     fprintf (fid,'%f      ',cell_matrix{i}(1,j));
 
     end

     fprintf (fid,'%d      ',num);
     fprintf (fid,'%d', size1-i+1);
     fprintf (fid,'\n\n');
end
  
