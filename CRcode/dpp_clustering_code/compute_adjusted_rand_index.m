function [rand,numerator,denum] = compute_adjusted_rand_index (conf_mat)

  n = sum(sum(conf_mat))

  sum_rows = sum (conf_mat,2)

  sum_columns = sum(conf_mat,1)

size_rows = size(conf_mat,1);
 
    size_columns = size(conf_mat,2);

   % return;

  denum = 0;

  numerator = 0;

  for i=1:1:size_rows

    denum = denum + my_nchoosek(sum_rows(i),2)/2;

    for j=1:1:size_columns

      if (i == 1)

            denum = denum + my_nchoosek(sum_rows(j),2)/2;

       end

    conf_mat (i,j)

    sum_rows(i)

      sum_columns (j)

      numerator = numerator + my_nchoosek (conf_mat (i,j), 2) - my_nchoosek(sum_rows(i),2)*my_nchoosek(sum_columns(j),2)/my_nchoosek(n,2);


      denum - denum -   my_nchoosek(sum_rows(i),2)*my_nchoosek(sum_columns(j),2)/my_nchoosek(n,2);

    end

 end


       rand = numerator/denum;
