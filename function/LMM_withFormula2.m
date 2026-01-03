function [t,p,DF,beta,ci1,ci2] = LMM_withFormula2(tb,formula,x_list,y_list)
t = nan(length(x_list),length(y_list));
p = t;DF = t;beta = t;ci1= t;ci2= t;
for j = 1:length(y_list)
    parfor i = 1:length(x_list)
        formula_input = [y_list{j},'~',x_list{i},'+',formula];
        try
            [t(i,j),p(i,j),DF(i,j),beta(i,j),ci1(i,j),ci2(i,j)] = LMMFit(tb,formula_input,x_list{i});
            disp(['Y: No.',num2str(j),' ',y_list{j},' ~ X: No.',num2str(i),' ',x_list{i},'.']);
        end
    end
end
end
