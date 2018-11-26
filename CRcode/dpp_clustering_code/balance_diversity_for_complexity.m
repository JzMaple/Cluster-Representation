function [out_matrix] =  balance_dicersity_for_complexity (in_matrix,alpha)


  dd = diag (in_matrix);

  size_dd = max (size (dd));

  q = zeros (1,size_dd);

  for i = 1:1:size_dd

    if (dd(i) < 0)

      dd(i)
      return;

    end

  q(i) = sqrt (dd(i));

  end

  for i = 1:1:size_dd
    for j = 1:1:size_dd

    if (i < j)
      out_matrix(i,j) = in_matrix (i,j)*(q(i)^(alpha-1))*(q(j)^(alpha-1));
    else 
     out_matrix(i,j) = in_matrix (i,j)*(q(j)^(alpha-1))*(q(i)^(alpha-1));
    end


   end
  end
