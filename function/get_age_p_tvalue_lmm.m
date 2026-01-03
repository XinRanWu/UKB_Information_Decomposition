function [value_cell] = get_age_p_tvalue_lmm(tableName,column_name_cell)
value_cell = cell(length(column_name_cell),3);
for i = 1:length(column_name_cell)
    value_cell{i,1} = column_name_cell{i};
end
for i = 1:length(column_name_cell)
    formula = sprintf('%s~Sex_0_0+AgeAttend_2_0+BMI_2_0+HeadMotion_2_0+Race_1+Race_2+Race_3+Race_4+Vol_WB_TIV_2_0+(1|Centre_0_0)', column_name_cell{i});
    model = fitlme(tableName, formula);
    value_cell{i,2} = model.Coefficients{4,6}; % pvalue
    value_cell{i,3} = model.Coefficients{4,4}; % tvalue
end
end

