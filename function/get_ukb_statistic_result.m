function get_ukb_statistic_result(ukb_RedSyn_table,UKB_UsedData_Brain,formula,XListAll,tablePath)
cd(tablePath);
UKB_BrainXY = innerjoin(UKB_UsedData_Brain,ukb_RedSyn_table(:,[1,8:end]));
UKB_BrainXY(UKB_BrainXY.HeadMotion_2_0>0.2,:) = [];
YList = ukb_RedSyn_table.Properties.VariableNames(8:end)'; 
[RedSynT,RedSynP,RedSynDF,RedSynBeta,~,~] = LMM_withFormula2(UKB_BrainXY,formula,XListAll,YList);
for i = 1:length(YList)
    RedSynLMMRes = [table(XListAll),array2table([RedSynBeta(:,i),RedSynT(:,i),RedSynDF(:,i),RedSynP(:,i)])];
    RedSynLMMRes.Properties.VariableNames(1:5) = {'X','Beta','TVal','DF','PVal'};
    writetable(RedSynLMMRes,[YList{i},'_forManhattan.csv']);
end
end

