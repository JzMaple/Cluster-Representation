function [result] = calculate_missClass_results (kmeans_trans_mat, index, gold_missClass_1_100, gold_missClass_101_200, gold_missClass_201_300, gold_missClass_301_400, gold_missClass_401_500, gold_missClass_501_600, gold_missClass_601_692 )


 

     conf_mat = create_confusion_matrix (kmeans_trans_mat{1}, gold_missClass_1_100{index});

     [a,b,result(1)] = calculate_hungarian_munkers (conf_mat(1:4,1:5));



     conf_mat = create_confusion_matrix (kmeans_trans_mat{2}, gold_missClass_101_200{index});

     [a,b,result(2)] = calculate_hungarian_munkers (conf_mat(1:4,1:5));




     conf_mat = create_confusion_matrix (kmeans_trans_mat{3}, gold_missClass_201_300{index});

     [a,b,result(3)] = calculate_hungarian_munkers (conf_mat(1:4,1:5));



     conf_mat = create_confusion_matrix (kmeans_trans_mat{3}, gold_missClass_201_300{index});

     [a,b,result(3)] = calculate_hungarian_munkers (conf_mat(1:4,1:5));



     conf_mat = create_confusion_matrix (kmeans_trans_mat{4}, gold_missClass_301_400{index});

     [a,b,result(4)] = calculate_hungarian_munkers (conf_mat(1:4,1:5));



     conf_mat = create_confusion_matrix (kmeans_trans_mat{5}, gold_missClass_401_500{index});

     [a,b,result(5)] = calculate_hungarian_munkers (conf_mat(1:4,1:5));



     conf_mat = create_confusion_matrix (kmeans_trans_mat{6}, gold_missClass_501_600{index});

     [a,b,result(6)] = calculate_hungarian_munkers (conf_mat(1:4,1:5));



     conf_mat = create_confusion_matrix (kmeans_trans_mat{7}, gold_missClass_601_692{index});

     [a,b,result(7)] = calculate_hungarian_munkers (conf_mat(1:4,1:5));




