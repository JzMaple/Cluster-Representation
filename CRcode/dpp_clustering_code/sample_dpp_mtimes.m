function [Y,stat,probs] = sample_dpp_mtimes(L,k,m)

  Y = zeros (m,k);

  for i=1:1:m 

    if (rem(i,1000) == 0)
      i = i
    end


    [temp,ret_flag] = sample_dpp (L,k);

    if (ret_flag == 0)
      Y(i,:) = temp';
      matrix = L.M (Y(i,:),Y(i,:));
      probs (i) = det (matrix);
    end

  end


  a = max(max(Y));

  stat = zeros(a,2);

  for i=1:1:a
    
    stat(i,1) = i;

    stat (i,2) = size(find(Y == i),1);

  end
