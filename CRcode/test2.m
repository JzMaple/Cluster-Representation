for j = 1:length(testPara)
    acc_ave_cutmatch_gm_200_50{j}=0;
    acc_ave_cutmatch_gm_150_30{j}=0;
    acc_ave_cutmatch_gm_100_20{j}=0;
    acc_ave_cutmatch_gc_200_50{j}=0;
    acc_ave_cutmatch_gc_150_30{j}=0;
    acc_ave_cutmatch_gc_100_20{j}=0;
    acc_ave_rrwm_gm{j}=0;
    acc_ave_ipfp_gm{j}=0;
    acc_ave_rawcut_gc{j}=0;
    acc_ave_ibgp_gm{j}=0;
end

for j = 1:length(testPara)
    for i = 1:nTest
        acc_ave_cutmatch_gm_200_50{j} = acc_ave_cutmatch_gm_200_50{j} + acc_matching_200_50{i,j};
        acc_ave_cutmatch_gm_150_30{j} = acc_ave_cutmatch_gm_150_30{j} + acc_matching_150_30{i,j};
        acc_ave_cutmatch_gm_100_20{j} = acc_ave_cutmatch_gm_100_20{j} + acc_matching_100_20{i,j};
        acc_ave_cutmatch_gc_200_50{j} = acc_ave_cutmatch_gc_200_50{j} + acc_cut_200_50{i,j};
        acc_ave_cutmatch_gc_150_30{j} = acc_ave_cutmatch_gc_150_30{j} + acc_cut_150_30{i,j};
        acc_ave_cutmatch_gc_100_20{j} = acc_ave_cutmatch_gc_100_20{j} + acc_cut_100_20{i,j};
        acc_ave_rrwm_gm{j} = acc_ave_rrwm_gm{j} + acc_matching_rrwm{i,j};
        acc_ave_rawcut_gc{j} = acc_ave_rawcut_gc{j} + acc_cut_raw{i,j};
        acc_ave_ibgp_gm{j} = acc_ave_ibgp_gm{j} + acc_matching_IBGP{i,j};
        acc_ave_ipfp_gm{j} = acc_ave_ipfp_gm{j} + acc_matching_IPFP{i,j};
    end
    acc_ave_cutmatch_gm_200_50{j} = acc_ave_cutmatch_gm_200_50{j}/nTest;
    acc_ave_cutmatch_gm_150_30{j} = acc_ave_cutmatch_gm_150_30{j}/nTest;
    acc_ave_cutmatch_gm_100_20{j} = acc_ave_cutmatch_gm_100_20{j}/nTest;
    acc_ave_cutmatch_gc_200_50{j} = acc_ave_cutmatch_gc_200_50{j}/nTest;
    acc_ave_cutmatch_gc_150_30{j} = acc_ave_cutmatch_gc_150_30{j}/nTest;
    acc_ave_cutmatch_gc_100_20{j} = acc_ave_cutmatch_gc_100_20{j}/nTest;
    acc_ave_rrwm_gm{j} = acc_ave_rrwm_gm{j}/nTest;
    acc_ave_rawcut_gc{j} = acc_ave_rawcut_gc{j}/nTest;
    acc_ave_ibgp_gm{j} = acc_ave_ibgp_gm{j}/nTest;
    acc_ave_ipfp_gm{j} = acc_ave_ipfp_gm{j}/nTest;
end


testPara = testPara; % testPara = testPara(end:-1:1); % comment this line if non-negative

acc_ave_cutmatch_gm_200_50_v = cell2mat(acc_ave_cutmatch_gm_200_50);
acc_ave_cutmatch_gm_150_30_v = cell2mat(acc_ave_cutmatch_gm_150_30);
acc_ave_cutmatch_gm_100_20_v = cell2mat(acc_ave_cutmatch_gm_100_20);
acc_ave_ibgp_gm_v = cell2mat(acc_ave_ibgp_gm);
acc_ave_rrwm_gm_v = cell2mat(acc_ave_rrwm_gm);
acc_ave_ipfp_gm_v = cell2mat(acc_ave_ipfp_gm);
figure; hold on; title('Matching performance on synthetic data');
p1 = plot(testPara,acc_ave_cutmatch_gm_200_50_v,'LineWidth',2,'LineStyle','-.','Marker','o','MarkerSize',6,'Color','b');
p2 = plot(testPara,acc_ave_cutmatch_gm_150_30_v,'LineWidth',2,'LineStyle','-','Marker','s','MarkerSize',6,'Color','y');
p3 = plot(testPara,acc_ave_cutmatch_gm_100_20_v,'LineWidth',2,'LineStyle','-.','Marker','d','MarkerSize',6,'Color','k');
p4 = plot(testPara,acc_ave_ibgp_gm_v,'LineWidth',2,'LineStyle','--','Marker','*','MarkerSize',6,'Color','g');
p5 = plot(testPara,acc_ave_rrwm_gm_v,'LineWidth',2,'LineStyle','-','Marker','p','MarkerSize',6,'Color','r');
p6 = plot(testPara,acc_ave_ipfp_gm_v,'LineWidth',2,'LineStyle','--','Marker','+','MarkerSize',6,'Color','m');
legend([p1,p2,p3,p4,p5,p6],'CutMatch-200-50','CutMatch-150-30','CutMatch-100-20','IBGP','RRWM','IPFP','Location','Best');
ylabel('\fontname{times new roman}Matching accuracy','FontSize',15);
xlabel('\fontname{times new roman}Feature permutation \it\rho','FontSize',15);
text(0.02,0.5,'\fontname{times new roman}graph deformation \it\sigma=0.1','FontSize',15);
text(0.02,0.44,'\fontname{times new roman}overlapping \it\gamma=0.1','FontSize',15);
text(0.02,0.38,'\fontname{times new roman}feature noise \it\mu=0.1','FontSize',15);
% text(100,100,'\fontname{times new roman}test \it\sigma = 0.1');
set(gca,'xtick',testPara);
hold off;

acc_ave_cutmatch_gc_200_50_v = cell2mat(acc_ave_cutmatch_gc_200_50);
acc_ave_cutmatch_gc_150_30_v = cell2mat(acc_ave_cutmatch_gc_150_30);
acc_ave_cutmatch_gc_100_20_v = cell2mat(acc_ave_cutmatch_gc_100_20);
acc_ave_rawcut_gc_v = cell2mat(acc_ave_rawcut_gc);
figure; hold on; title('Cuts performance on synthetic data');
p1 = plot(testPara,acc_ave_cutmatch_gc_200_50_v,'LineWidth',2,'LineStyle','-.','Marker','o','MarkerSize',6,'Color','b');
p2 = plot(testPara,acc_ave_cutmatch_gc_150_30_v,'LineWidth',2,'LineStyle','-','Marker','s','MarkerSize',6,'Color','y');
p3 = plot(testPara,acc_ave_cutmatch_gc_100_20_v,'LineWidth',2,'LineStyle','-.','Marker','d','MarkerSize',6,'Color','k');
p4 = plot(testPara,acc_ave_rawcut_gc_v,'LineWidth',2,'LineStyle','--','Marker','*','MarkerSize',6,'Color','g');
legend([p1,p2,p3,p4],'CutMatch-200-50','CutMatch-150-30','CutMatch-100-20','Raw Cuts','Location','Best');
ylabel('\fontname{times new roman}Cuts accuracy','FontSize',15);
xlabel('\fontname{times new roman}Feature permutation \it\rho','FontSize',15);
text(0.02,0.87,'\fontname{times new roman}graph deformation \it\sigma=0.1','FontSize',15);
text(0.02,0.88,'\fontname{times new roman}overlapping \it\gamma=0.1','FontSize',15);
text(0.02,0.89,'\fontname{times new roman}feature noise \it\mu=0.1','FontSize',15);
set(gca,'xtick',testPara);
hold off;

a = 1;