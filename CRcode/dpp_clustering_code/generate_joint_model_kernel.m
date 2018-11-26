function [l12]  = generate_joint_model_kernel (matrix1,matrix2)

  size1 = size(matrix1,1);

  size2 = size(matrix1,2);

  size3 = size(matrix2,1);

  size4 = size(matrix2,2);


  if ( (size1 ~= size3) | (size2 ~= size4))
    fprintf (2,'ERROR !\n');
    return;
  end


  q1 = sqrt(diag(matrix1))' %'

  q2 = sqrt(diag(matrix2))' %'

  for i=1:1:size1
    for j=1:1:size2
    
      if (q1(i)*q1(j) > 0)
       s1(i,j) = matrix1(i,j)/(q1(i)*q1(j));
      end

      if (q2(i)*q2(j))
       s2(i,j) = matrix2(i,j)/(q2(i)*q2(j));
      end

   end
end


   s12 = s1*s2;

   s21 = s2*s1;

   S =  s1*s2*s2*s1;


size (find (s1 ~= s1')) %'

size (find (s2 ~= s2')) %'


  size (find (S ~= S')) %'

 size (find (S == S')) %'


S = (S + S')./2 ; %'


%sum(sum(S - S')) %'

  size (find (S ~= S')) %'

 for i=1:1:size1
    for j=1:1:size2

	      l12(i,j) = S(i,j)*(q1(i)^2)*(q1(j)^2)*(q2(i)^2)*(q2(j)^2);

    end
end

%  size (find (l12 ~= l12')) %'

	      l12 = (l12 + l12')./2;%'

  size (find (l12 ~= l12')) %'

%	[v,d] = eig(l12);

%	diag(d)'
