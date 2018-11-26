function [conf_mat] = create_confusion_matrix_zeros (vec1, vec2)

  size_vec1 = max(size(vec1));

  size_vec2 = max (size(vec2));
 
  size_vec1

  size_vec2 


 if (size_vec1 ~= size_vec2)
   return;
 end

  max_vec1 = max(max(vec1))

  max_vec2 = max(max(vec2))

  conf_mat = zeros(max_vec1, max_vec2);

  for(i=1:size_vec1)

    if (vec1(i) ~= 0)

    conf_mat (vec1(i), vec2(i)) = conf_mat (vec1(i), vec2(i)) + 1;

   end

  end
