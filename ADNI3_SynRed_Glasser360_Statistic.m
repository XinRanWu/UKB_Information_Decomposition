load('ADNI3_tabluar.mat')
PTDEMOG_03Sep2024 = removevars(PTDEMOG_03Sep2024, {'VISCODE','VISCODE2','PHASE','RID','PTSOURCE','PTDOBYY','PTWORK','PTCOGBEG',...
    'PTADDX','PTIDENT','PTENGSPK','PTNLANG','PTENGSPKAGE','PTCLANG','PTLANGSP','PTLANGWR','PTSPTIM','PTSPOTTIM','PTLANGPR1',...
    'PTLANGSP1','PTLANGRD1','PTLANGWR1','PTLANGUN1','PTLANGPR2','PTLANGSP2','PTLANGRD2','PTLANGWR2','PTLANGUN2','PTLANGPR3',...
    'PTLANGSP3','PTLANGRD3','PTLANGWR3','PTLANGUN3','PTLANGPR4','PTLANGSP4','PTLANGRD4','PTLANGWR4','PTLANGUN4','PTLANGPR5',...
    'PTLANGSP5','PTLANGRD5','PTLANGWR5','PTLANGUN5','PTLANGPR6','PTLANGSP6','PTLANGRD6','PTLANGWR6','PTLANGUN6','PTLANGTTL',...
    'PTETHCATH','PTBORN','PTBIRPL','PTIMMAGE','PTIMMWHY','PTBIRPR','PTBIRGR','USERDATE2','DD_CRF_VERSION_LABEL','LANGUAGE_CODE',...
    'HAS_QC_ERROR','USERDATE','update_stamp'});
MMSE_03Sep2024 = removevars(MMSE_03Sep2024, {'PHASE','RID','VISCODE','VISCODE2','DONE','NDREASON','SOURCE','MMDATE','MMYEAR',...
    'MMMONTH','MMDAY','MMSEASON','MMHOSPIT','MMFLOOR','MMCITY','MMAREA','MMSTATE','WORDLIST','WORD1','WORD2','WORD3','MMTRIALS',...
    'MMD','MML','MMR','MMO','MMW','MMLTR1','MMLTR2','MMLTR3','MMLTR4','MMLTR5','MMLTR6','MMLTR7','WORLDSCORE','WORD1DL','WORD2DL',...
    'WORD3DL','MMWATCH','MMPENCIL','MMREPEAT','MMHAND','MMFOLD','MMONFLR','MMREAD','MMWRITE','MMDRAW','ID','SITEID','USERDATE',...
    'USERDATE2','DD_CRF_VERSION_LABEL','LANGUAGE_CODE','HAS_QC_ERROR','update_stamp'});


MAYOADIRL_MRI_QUALITY_ADNI3_03Sep2024 = removevars(MAYOADIRL_MRI_QUALITY_ADNI3_03Sep2024, {'RID',...
    'SERIES_DATE','STUDYINSTANCEUID','STUDY_OVERALLPASS','STUDY_COMMENTS','SERIESINSTANCEUID',...
    'SERIES_COMMENTS','PROTOCOL_COMMENTS','LONI_STUDY','LONI_SERIES','LONI_IMAGE','STUDY_MEDICAL_ABNORMALITIES',...
    'SERIAL','MEDICAL_EXCLUSION','RELEASE_FROM_QUARANTINE','PAY_SITE','QUALIFICATION','FIELD_STRENGTH',...
    'T1_ACCELERATED'});

APOERES_03Sep2024 = removevars(APOERES_03Sep2024, {'PHASE','RID','VISCODE','APTESTDT','APVOLUME','APRECEIVE',...
    'APAMBTEMP','APRESAMP','APUSABLE','ID','SITEID','USERDATE','USERDATE2','update_stamp'});


PTDEMOG_03Sep2024 = removevars(PTDEMOG_03Sep2024, 'VISDATE');
PTDEMOG_03Sep2024 = removevars(PTDEMOG_03Sep2024, {'ID','SITEID'});
% PTDEMOG_03Sep2024(1,:) = [];
PTDEMOG_03Sep2024 = sortrows(PTDEMOG_03Sep2024,'PTID','ascend');
PTDEMOG_03Sep2024(PTDEMOG_03Sep2024.PTGENDER==NaN,:) = [];
PTDEMOG_03Sep2024(PTDEMOG_03Sep2024.PTGENDER==-4,:) = [];
PTDEMOG_03Sep2024(isnan(PTDEMOG_03Sep2024.PTGENDER),:) = [];
PTDEMOG_03Sep2024 = removevars(PTDEMOG_03Sep2024, 'PTWORKHS');
PTDEMOG_03Sep2024.PTADBEG(PTDEMOG_03Sep2024.PTADBEG==-4) = NaN;



ADNI3_with_PET_9_03_2024.Properties.VariableNames{2} = 'PTID';
ADNI3_with_PET_9_03_2024.PTID = string(ADNI3_with_PET_9_03_2024.PTID);


ADNI3_fMRI_info.Properties.VariableNames{1} = 'PTID';
ADNI3_fMRI_info = removevars(ADNI3_fMRI_info, {'Type','Format','Downloaded','visit_yr'});
ADNI3_fMRI_info = removevars(ADNI3_fMRI_info, 'Modality');


NEUROPATH_02_06_23_03Sep2024 = removevars(NEUROPATH_02_06_23_03Sep2024, {'RID','update_stamp'});
GDSCALE_03Sep2024 = removevars(GDSCALE_03Sep2024, {'PHASE','RID','VISCODE','VISCODE2','SOURCE','ID','SITEID','USERDATE','USERDATE2','DD_CRF_VERSION_LABEL','LANGUAGE_CODE','HAS_QC_ERROR','update_stamp'});
GDSCALE_03Sep2024 = removevars(GDSCALE_03Sep2024, {'GDUNABL','GDSATIS','GDDROP','GDEMPTY','GDBORED','GDSPIRIT','GDAFRAID','GDHAPPY','GDHELP','GDHOME','GDMEMORY','GDALIVE','GDWORTH','GDENERGY','GDHOPE','GDBETTER'});
NPIQ_03Sep2024 = removevars(NPIQ_03Sep2024, {'PHASE','RID','VISCODE','VISCODE2','NPDATE','SOURCE'});
NPIQ_03Sep2024 = removevars(NPIQ_03Sep2024, {'NPIA','NPIASEV','NPIB','NPIBSEV','NPIC','NPICSEV','NPID','NPIDSEV','NPIE','NPIESEV','NPIF','NPIFSEV','NPIG','NPIGSEV','NPIH','NPIHSEV','NPII','NPIISEV','NPIJ','NPIJSEV','NPIK','NPIKSEV','NPIL','NPILSEV'});
NPIQ_03Sep2024 = removevars(NPIQ_03Sep2024, {'SPID','ID','SITEID','USERDATE','USERDATE2','DD_CRF_VERSION_LABEL','LANGUAGE_CODE','HAS_QC_ERROR','update_stamp'});
CDR_03Sep2024 = removevars(CDR_03Sep2024, {'PHASE','RID','VISCODE','VISCODE2','CDSOURCE','CDVERSION','SPID','SITEID','USERDATE','USERDATE2','DD_CRF_VERSION_LABEL','LANGUAGE_CODE','HAS_QC_ERROR','update_stamp','ID'});
CDR_03Sep2024 = removevars(CDR_03Sep2024, {'CDMEMORY','CDORIENT','CDJUDGE','CDCOMMUN','CDHOME','CDCARE','CDGLOBAL'});



CDR_03Sep2024.PTID = string(CDR_03Sep2024.PTID);
NPIQ_03Sep2024.PTID = string(NPIQ_03Sep2024.PTID);
GDSCALE_03Sep2024.PTID = string(GDSCALE_03Sep2024.PTID);

ADNI_INFO = outerjoin(PTDEMOG_03Sep2024,APOERES_03Sep2024,'MergeKeys',1);
ADNI_SCORES = outerjoin(outerjoin(outerjoin(MMSE_03Sep2024,...
    GDSCALE_03Sep2024,'MergeKeys',1),...
    NPIQ_03Sep2024,'MergeKeys',1),...
    CDR_03Sep2024,'MergeKeys',1);

ADNI_INFO(ADNI_INFO.PTID=="",:) = [];

ADNI3_fMRI_info.PTID = string(ADNI3_fMRI_info.PTID);


ADNI3_fMRI_data_qc.mean_fd = mean_fd_qc;
ADNI3_fMRI_data_qc.Properties.VariableNames{1} = 'PTID';
ADNI3_fMRI_data_qc.PTID = string(ADNI3_fMRI_data_qc.PTID);

ADNI_INFO_fMRI = innerjoin(ADNI_INFO,ADNI3_fMRI_info);
% ADNI_INFO_fMRI = innerjoin(innerjoin(ADNI_INFO,ADNI3_fMRI_info),ADNI3_fMRI_data_qc);

load('/public/home/zhangjie/ZJLab/ADNI_Project/data/imaging_data/ADNI3_fMRI_redunancy_synergy.mat')
load('/public/home/zhangjie/ZJLab/UKBiobank_Project/project/Detecting_Compensation_Effects_of_Aging_Using_PID/data/atlas/HCPMMP_atlas_info.mat', 'Yeo7MMP')






synergy_mean = mean(icatb_mat2vec3d(synergy),2);
redundancy_mean = mean(icatb_mat2vec3d(redundancy),2);

ADNI3_fMRI_info = ADNI3_with_PET_9_03_2024(strcmp(ADNI3_with_PET_9_03_2024.Modality,'fMRI'),:);
ADNI3_fMRI_info = sortrows(ADNI3_fMRI_info,'Subject','ascend');
Age = table2cell(ADNI3_fMRI_info(:,5));
ADNI3_fMRI_info = removevars(ADNI3_fMRI_info, 'ImageDataID');
for i = 1:size(ADNI3_fMRI_info,1)
    AgeNum(i,1) = str2num(Age{i});
end
ADNI3_fMRI_info.Age = AgeNum;



ADNI3_Subject_info = ADNI3_with_PET_9_03_2024;
ADNI3_Subject_info = removevars(ADNI3_Subject_info, {'Description','Type','AcqDate','Format','Age'});
ADNI3_Subject_info = removevars(ADNI3_Subject_info, {'ImageDataID','Modality','Downloaded','Visit'});
ADNI3_Subject_info = unique(ADNI3_Subject_info);


date_col = datetime(ADNI3_fMRI_info.AcqDate, 'InputFormat', 'yyyy-MM-dd');
ADNI3_fMRI_info.visit_yr = year(date_col);

ADNI3_fMRI_info = unique(ADNI3_fMRI_info);

removeValues = {
    'Axial 2D PASL'
    'Axial 3D PASL (Eyes Open)'
    'Axial MB DTI'
    'Axial MB DTI_ADC'
    'Axial MB DTI_FA'
    'Axial MB DTI_TENSOR_B0'
    'Axial MB DTI_TRACEW'
    'MoCoSeries'
    'Perfusion_Weighted'
    'SOURCE - Axial 2D PASL'
    'WIP SOURCE - Axial 3D pCASL (Eyes Open)'
    'relCBF'
};

rowsToKeep = ~ismember(ADNI3_fMRI_info.Description, removeValues);
ADNI3_fMRI_info = ADNI3_fMRI_info(rowsToKeep, :);

ADNI3_fMRI_info.visit_yr(452) = 2020;
ADNI3_fMRI_info(513,:) = [];
ADNI3_fMRI_info.visit_yr(650) = 2020;
ADNI3_fMRI_info(719,:) = [];
ADNI3_fMRI_info.visit_yr(859) = 2020;
ADNI3_fMRI_info.visit_yr(860) = 2021;
ADNI3_fMRI_info.visit_yr(1041) = 2020;
ADNI3_fMRI_info(1048,:) = [];
ADNI3_fMRI_info.visit_yr(1179) = 2020;
ADNI3_fMRI_info.visit_yr(1182) = 2020;
ADNI3_fMRI_info.visit_yr(1219) = 2021;



for i = 1:size(ADNI3_fMRI_data_qc,1)
    if sum(strcmp(ADNI3_Subject_info.Subject,ADNI3_fMRI_data_qc.Subject(i)))>0
        
        matchingRows(i,1) = find(strcmp(ADNI3_Subject_info.Subject,ADNI3_fMRI_data_qc.Subject(i)));
    else
        matchingRows(i,1) = NaN;
    end
    
end


ADNI3_fMRI_data_qc_group = outerjoin(ADNI3_fMRI_data_qc,ADNI3_Subject_info,'MergeKeys',1);
ADNI3_fMRI_data_qc_group = ADNI3_fMRI_data_qc_group(~isnan(ADNI3_fMRI_data_qc_group.visit_yr),:);


red_group(1,1) = mean(redundancy_mean(strcmp(ADNI3_fMRI_data_qc_group.Group,'CN')));
red_group(2,1) = mean(redundancy_mean(strcmp(ADNI3_fMRI_data_qc_group.Group,'SMC')));
red_group(3,1) = mean(redundancy_mean(strcmp(ADNI3_fMRI_data_qc_group.Group,'EMCI')));
red_group(4,1) = mean(redundancy_mean(strcmp(ADNI3_fMRI_data_qc_group.Group,'MCI')));
red_group(5,1) = mean(redundancy_mean(strcmp(ADNI3_fMRI_data_qc_group.Group,'LMCI')));
red_group(6,1) = mean(redundancy_mean(strcmp(ADNI3_fMRI_data_qc_group.Group,'AD')));


syn_group(1,1) = mean(synergy_mean(strcmp(ADNI3_fMRI_data_qc_group.Group,'CN')));
syn_group(2,1) = mean(synergy_mean(strcmp(ADNI3_fMRI_data_qc_group.Group,'SMC')));
syn_group(3,1) = mean(synergy_mean(strcmp(ADNI3_fMRI_data_qc_group.Group,'EMCI')));
syn_group(4,1) = mean(synergy_mean(strcmp(ADNI3_fMRI_data_qc_group.Group,'MCI')));
syn_group(5,1) = mean(synergy_mean(strcmp(ADNI3_fMRI_data_qc_group.Group,'LMCI')));
syn_group(6,1) = mean(synergy_mean(strcmp(ADNI3_fMRI_data_qc_group.Group,'AD')));

Group = {'CN','SMC','EMCI','MCI','LMCI','AD'};
Yeo7Names = {'Visual';'Somatomotor';'DorsalAttention';'VentralAttention';'Limbic';'Frontoparietal';'Default'};
for i = 1:6
    X = mean(synergy(:,:,strcmp(ADNI3_fMRI_data_qc_group.Group,Group{i})),3);
    subplot(2,3,i);plot_connmatrix(X([181:360,1:180],[181:360,1:180]),...
        Yeo7MMP,[Yeo7Names],[0,0.2],'jet');colormap(jet(200));
    title(Group{i});
end

for i = 1:6
    X = mean(redundancy(:,:,strcmp(ADNI3_fMRI_data_qc_group.Group,Group{i})),3);
    subplot(2,3,i);plot_connmatrix(X([181:360,1:180],[181:360,1:180]),...
        Yeo7MMP,[Yeo7Names],[0,0.12],'jet');colormap(jet(200));
    title(Group{i});
end










roi_tmp = parcel_to_surface(1:360,'glasser_360_fsa5');
for i = 1:360
    roi_size(i,1) = sum(roi_tmp==i);
end

parfor i = 1:length(ts_qc)
        [L,R] = AsymCalc_HemiGS(ts_qc{i},roi_size',181:360,1:180);
        adni3_glasser_LI(i,:) = L - R;
        adni3_glasser_dLI(i,:) = std(DynLaterIndex(ts_qc{i},roi_size,181:360,1:180,30,1));
        disp(['Temporal Lateralization (Cort) of ADNI3: sub-',...
            ADNI3_fMRI_data_qc.Subject{i},' ses-',num2str(ADNI3_fMRI_data_qc.visit_yr(i)),'.'])
end
