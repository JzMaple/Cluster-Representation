function r = Representation(iCluster,jCluster,score,consistency,N)
    r = [iCluster(1) jCluster(1)];
    scr = score(1,1);
    con = consistency(1,1);
    for i = iCluster
        for j = jCluster
            if i<N && j <N
                if score(i,j) == 1
                    if scr == 1 && con < consistency(i,j)
                        r = [i j];
                        con = consistency(i,j);
                    end
                    if scr ~= 1
                        r = [i j];
                        con = consistency(i,j);
                        scr = 1;
                    end
                else
                    if score(i,j)>0.9 && consistency(i,j)>0.5 && score(i,j)*0.6 + consistency(i,j) * 0.4 > scr*0.6 + con*0.4 
                        r = [i j];
                        con = consistency(i,j);
                        scr = score(i,j);
                    end
                end
            end
        end
    end
end