%[text] Scientific Programming - TCGA Ovarian Cancer dataset
%[text] 01 - Import and clean the dataset
clear;
clc;
close all;

if ~exist("figures", "dir")
    mkdir("figures");
end
%%
%[text] Import clinical files and gene expression data
clinical_patient = readtable( ...
    "data/data_clinical_patient.txt", ...
    "FileType","text", ...
    "Delimiter","\t", ...
    "CommentStyle","#");

clinical_sample = readtable( ...
    "data/data_clinical_sample.txt", ...
    "FileType","text", ...
    "Delimiter","\t", ...
    "CommentStyle","#");

expression_data = readtable( ... %[output:group:18723482] %[output:4c3b5eca]
    "data/data_mRNA_unified.txt", ... %[output:4c3b5eca]
    "FileType","text", ... %[output:4c3b5eca]
    "Delimiter","\t"); %[output:group:18723482] %[output:4c3b5eca]
%[text] **Note: Gene-expression preprocessing**
%[text] The unified TCGA mRNA expression profile was used. The individual platforms were normalized separately: the Affymetrix platforms were processed using RMA, while the Agilent platform used LOWESS normalization followed by log2 transformation. The platform-specific expression values were then standardized across patients and integrated across the three platforms using factor analysis, after which the unified estimates were rescaled using gene-specific variability estimates.
%[text] Batch effects were detected in the individual platforms but were substantially reduced in the unified dataset. The authors checked this using permutation-based batch testing, PCA and hierarchical clustering, and chose not to apply further batch correction because empirical-Bayes adjustment produced only modest changes and could remove genuine biological differences between batches. 
%[text] Source: [https://doi.org/10.1038/nature10166](https://doi.org/10.1038/nature10166) (Supplementary Information)
%%
%[text] Inspect the datasets: variable names, number of rows and columns.
size(clinical_patient) %491 x 7 %[output:284bb327]
size(clinical_sample) %488 x 15 %[output:0fce9570]
size(expression_data) %11838 x 497 %[output:95544765]

clinical_patient.Properties.VariableNames %[output:64520e29]
clinical_sample.Properties.VariableNames %[output:8de36f0c]
expression_data.Properties.VariableNames(1:5) %[output:6d41e9e2]
%%
%[text] Inspect the first rows: get patient ID format, sample ID format, how to match the samples from different files
head(clinical_patient) %[output:5061b037]
head(clinical_sample) %[output:4d1ee97a]
expression_data(1:5,1:6) %[output:052a0b64]

% Patient ID column
clinical_patient.PATIENT_ID(1:5) %[output:6991c56d]

clinical_sample.PATIENT_ID(1:5) %[output:6b515ea6]

clinical_sample.SAMPLE_ID(1:5) %[output:21cf3a32]
%[text] Inspect if there are unique IDs and how much overlap there is.
% Uniqueness of patient IDs
numel(unique(clinical_patient.PATIENT_ID)) %[output:902b3b5b]
height(clinical_patient) %[output:42cb7785]

% Uniqueness of sample IDs
numel(unique(clinical_sample.SAMPLE_ID)) %[output:4a0337b5]
height(clinical_sample) %[output:13c0eb10]

% Uniqueness of patient IDs in sample table
numel(unique(clinical_sample.PATIENT_ID)) %[output:408bc634]
height(clinical_sample) %[output:57306e8a]

%no repeated patient IDs or sample IDs in any of the files

%check which patients dont have matching sample
patientOnly = setdiff( ...
    clinical_patient.PATIENT_ID, ...
    clinical_sample.PATIENT_ID);

patientOnly %[output:407125ac]

%3 patients do not have corresponding sample IDs
%[text] There are in total 3 patients that don't have a corresponding sample ID, so they are excluded from further analysis
clinical_patient = clinical_patient( ...
    ismember(clinical_patient.PATIENT_ID, clinical_sample.PATIENT_ID), ...
    :);

% Keeps only patients with a corresponding sample-level record

%check size now
size(clinical_patient) %488 x 7 %[output:5c48fc7c]
%[text] From 491 samples it decreased to 488 samples, this confirms that the 3 patients that were unmatched were removed from the dataset.
%%
%[text] Match and merge the patient datasets by patient\_ID
clinical_merged = innerjoin( ...
    clinical_sample, ...
    clinical_patient, ...
    "Keys","PATIENT_ID");

% Check dimensions
size(clinical_merged) %488 x 21 %[output:79cd82ac]
%%
%[text] Match clinical samples to gene expression samples (sample\_ID)
expressionSampleIDs = string( ...
    expression_data.Properties.VariableNames(2:end));

% Convert underscores to TCGA barcode format
expressionSampleIDs = strrep(expressionSampleIDs, "_", "-");

% Shorten barcode: TCGA-01-0630-01-11 -> TCGA-01-0630-01
expressionSampleIDs_short = extractBefore(expressionSampleIDs, 16);

% Check that shortening did not create duplicate IDs
numel(unique(expressionSampleIDs_short)) %[output:1a539e95]
numel(expressionSampleIDs_short) %[output:60a29bfb]

clinicalSampleIDs = string(clinical_merged.SAMPLE_ID);

% Match clinical samples to expression samples
[isMatched, expressionOrder] = ...
    ismember(clinicalSampleIDs, expressionSampleIDs_short);

sum(isMatched) %[output:33f8b243]

% Check unmatched samples
clinicalSampleIDs(~isMatched) %[output:360f3e3d]

% Extract expression values
expressionValues = expression_data{:,2:end};

% Reorder expression columns so they match clinical_merged
expressionValues_matched = ...
    expressionValues(:, expressionOrder);

% Store gene IDs
geneIDs = expression_data.Entrez_Gene_Id;

% Check result
size(expressionValues_matched) %[output:63df1b45]
%%
%[text] Check for missingness and outcome availability
% Check missing values in all clinical variables

for i = 1:width(clinical_merged) %[output:group:666b529f]
    variableName = clinical_merged.Properties.VariableNames{i};
    nMissing = sum(ismissing(clinical_merged.(variableName)));

    fprintf('%s: %d missing\n', variableName, nMissing); %[output:029f1d89]
end %[output:group:666b529f]

% Inspect categorical values in important clinical variables

unique(clinical_merged.TUMOR_STAGE_2009) %[output:9af91628]
unique(clinical_merged.GRADE) %[output:03e136de]
unique(clinical_merged.TUMOR_RESIDUAL_DISEASE) %[output:0c4f6c2d]
unique(clinical_merged.PRIMARY_THERAPY_OUTCOME_SUCCESS) %[output:08333ba5]
unique(clinical_merged.PLATINUM_STATUS) %[output:51fadab2]
unique(clinical_merged.DFS_STATUS) %[output:61c686a1]
unique(clinical_merged.OS_STATUS) %[output:0fccea96]

tabulate(clinical_merged.TUMOR_STAGE_2009) %[output:5534b29f]
tabulate(clinical_merged.GRADE) %[output:506858a2]
tabulate(clinical_merged.TUMOR_RESIDUAL_DISEASE) %[output:058ac4e1]
tabulate(clinical_merged.PRIMARY_THERAPY_OUTCOME_SUCCESS) %[output:65aadeb2]
tabulate(clinical_merged.PLATINUM_STATUS) %[output:6c2fd3bc]
tabulate(clinical_merged.DFS_STATUS) %[output:90d9db0e]
%[text] Check missingness for DFS feature
dfsStatusAvailable = ~ismissing(clinical_merged.DFS_STATUS);
dfsTimeAvailable   = ~ismissing(clinical_merged.DFS_MONTHS);

% Both status and time available
sum(dfsStatusAvailable & dfsTimeAvailable) %[output:53886646]

% Status available but time missing
sum(dfsStatusAvailable & ~dfsTimeAvailable) %[output:38a61561]

% Time available but status missing
sum(~dfsStatusAvailable & dfsTimeAvailable) %[output:496debd0]

% Both missing
sum(~dfsStatusAvailable & ~dfsTimeAvailable) %[output:46b844fb]
%[text] Since DFS is one of the features I want to look at the most, I will only keep the samples that have both time and status available (395 samples). 
%DFS analysis cohort

dfsAvailable = ...
    ~ismissing(clinical_merged.DFS_STATUS) & ...
    ~ismissing(clinical_merged.DFS_MONTHS);

clinical_dfs = clinical_merged(dfsAvailable, :);
expression_dfs = expressionValues_matched(:, dfsAvailable);

size(clinical_dfs) %[output:98520319]
size(expression_dfs) %[output:06fc23f2]
%[text] Create event variable for DFS status
dfsEvent = strcmp(clinical_dfs.DFS_STATUS, ...
                  '1:Recurred/Progressed');

sum(dfsEvent)      % 283 events %[output:8dea4a8c]
sum(~dfsEvent)     % 112 censored %[output:1a4860f9]
%%
%[text] Visualization of missing data
%%\ Clinical missingness overview

% Variables used later in the survival analysis
stageStr = string(clinical_dfs.TUMOR_STAGE_2009);
gradeStr = string(clinical_dfs.GRADE);
residualStr = string(clinical_dfs.TUMOR_RESIDUAL_DISEASE);

% Count both true missing values and text encoded as "NA"
missingStage = ...
    sum(ismissing(stageStr) | stageStr == "NA");

missingGrade = ...
    sum(ismissing(gradeStr) | gradeStr == "NA");

missingResidual = ...
    sum(ismissing(residualStr) | residualStr == "NA");

missingDFS_Time = ...
    sum(isnan(clinical_dfs.DFS_MONTHS));

missingDFS_Status = ...
    sum(ismissing(string(clinical_dfs.DFS_STATUS)));

features = { ...
    'DFS time', ...
    'DFS status', ...
    'Stage', ...
    'Grade', ...
    'Residual disease'};

missingCounts = [ ...
    missingDFS_Time, ...
    missingDFS_Status, ...
    missingStage, ...
    missingGrade, ...
    missingResidual];

missingPercent = ...
    100 * missingCounts / height(clinical_dfs);

%%\ Visualize missing clinical data

figure; %[output:1bcf3488]

bar(missingPercent) %[output:1bcf3488]

xticklabels(features) %[output:1bcf3488]
xtickangle(30) %[output:1bcf3488]

ylabel('Missing values (%)') %[output:1bcf3488]
title('Missingness in clinical variables') %[output:1bcf3488]

grid on %[output:1bcf3488]

saveas(gcf, "figures/missingness_clinical.png"); %[output:1bcf3488]

missingTable = table( ...
    features', ...
    missingCounts', ...
    missingPercent', ...
    'VariableNames', ...
    {'Variable','MissingCount','MissingPercent'});

disp(missingTable) %[output:97ab26cf]
%%
%[text] Save data 
save("results/cleaned_data.mat", ...
    "clinical_merged", ...
    "clinical_dfs", ...
    "expressionValues_matched", ...
    "expression_dfs", ...
    "expression_data", ...
    "geneIDs", ...
    "dfsEvent");

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:4c3b5eca]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Column headers from the file were modified to make them valid MATLAB identifiers before creating variable names for the table. The original column headers are saved in the VariableDescriptions property.\nSet 'VariableNamingRule' to 'preserve' to use the original column headers as table variable names."}}
%---
%[output:284bb327]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["491","7"]]}}
%---
%[output:0fce9570]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["488","15"]]}}
%---
%[output:95544765]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["11838","497"]]}}
%---
%[output:64520e29]
%   data: {"dataType":"matrix","outputData":{"columns":7,"header":"1×7 cell array","name":"ans","rows":1,"type":"cell","value":[["'PATIENT_ID'","'PRIMARY_THERAPY_OUTCOME_SUCCESS'","'PLATINUM_STATUS'","'OS_STATUS'","'OS_MONTHS'","'DFS_STATUS'","'DFS_MONTHS'"]]}}
%---
%[output:8de36f0c]
%   data: {"dataType":"matrix","outputData":{"columns":15,"header":"1×15 cell array","name":"ans","rows":1,"type":"cell","value":[["'SAMPLE_ID'","'PATIENT_ID'","'TUMOR_STAGE_2009'","'GRADE'","'TUMOR_RESIDUAL_DISEASE'","'ACGH_DATA'","'MRNA_DATA'","'SEQUENCED'","'COMPLETE_DATA'","'CANCER_TYPE'","'CANCER_TYPE_DETAILED'","'ONCOTREE_CODE'","'SAMPLE_TYPE'","'SOMATIC_STATUS'","'TMB_NONSYNONYMOUS'"]]}}
%---
%[output:6d41e9e2]
%   data: {"dataType":"matrix","outputData":{"columns":5,"header":"1×5 cell array","name":"ans","rows":1,"type":"cell","value":[["'Entrez_Gene_Id'","'TCGA_01_0630_01_11'","'TCGA_01_0631_01_11'","'TCGA_01_0633_01_11'","'TCGA_01_0636_01_11'"]]}}
%---
%[output:5061b037]
%   data: {"dataType":"text","outputData":{"text":"       <strong>PATIENT_ID<\/strong>       <strong>PRIMARY_THERAPY_OUTCOME_SUCCESS<\/strong>    <strong>PLATINUM_STATUS<\/strong>     <strong>OS_STATUS<\/strong>      <strong>OS_MONTHS<\/strong>       <strong>DFS_STATUS<\/strong>        <strong>DFS_MONTHS<\/strong>\n    <strong>________________<\/strong>    <strong>_______________________________<\/strong>    <strong>_______________<\/strong>    <strong>____________<\/strong>    <strong>_________<\/strong>    <strong>_________________<\/strong>    <strong>__________<\/strong>\n\n    {'TCGA-61-1895'}         {'Stable Disease'   }          {'Tooearly'}      {'0:LIVING'}      1.41       {'0:DiseaseFree'}       1.41   \n    {'TCGA-13-1510'}         {0×0 char           }          {'Tooearly'}      {'0:LIVING'}      2.56       {'0:DiseaseFree'}       2.56   \n    {'TCGA-23-1809'}         {0×0 char           }          {'Tooearly'}      {'0:LIVING'}      0.53       {'0:DiseaseFree'}       0.53   \n    {'TCGA-20-1687'}         {'Complete Response'}          {'Tooearly'}      {'0:LIVING'}      2.66       {'0:DiseaseFree'}       2.66   \n    {'TCGA-13-1404'}         {'Complete Response'}          {'Tooearly'}      {'0:LIVING'}       6.6       {'0:DiseaseFree'}        6.6   \n    {'TCGA-13-1506'}         {'Complete Response'}          {'Tooearly'}      {'0:LIVING'}       4.4       {'0:DiseaseFree'}        4.4   \n    {'TCGA-23-1116'}         {0×0 char           }          {'Tooearly'}      {'0:LIVING'}      4.57       {'0:DiseaseFree'}       4.57   \n    {'TCGA-20-1686'}         {'Complete Response'}          {'Tooearly'}      {'0:LIVING'}      2.92       {'0:DiseaseFree'}       2.92   \n\n","truncated":false}}
%---
%[output:4d1ee97a]
%   data: {"dataType":"text","outputData":{"text":"         <strong>SAMPLE_ID<\/strong>            <strong>PATIENT_ID<\/strong>       <strong>TUMOR_STAGE_2009<\/strong>    <strong>GRADE<\/strong>       <strong>TUMOR_RESIDUAL_DISEASE<\/strong>      <strong>ACGH_DATA<\/strong>    <strong>MRNA_DATA<\/strong>    <strong>SEQUENCED<\/strong>    <strong>COMPLETE_DATA<\/strong>       <strong>CANCER_TYPE<\/strong>                <strong>CANCER_TYPE_DETAILED<\/strong>            <strong>ONCOTREE_CODE<\/strong>    <strong>SAMPLE_TYPE<\/strong>    <strong>SOMATIC_STATUS<\/strong>    <strong>TMB_NONSYNONYMOUS<\/strong>\n    <strong>___________________<\/strong>    <strong>________________<\/strong>    <strong>________________<\/strong>    <strong>______<\/strong>    <strong>__________________________<\/strong>    <strong>_________<\/strong>    <strong>_________<\/strong>    <strong>_________<\/strong>    <strong>_____________<\/strong>    <strong>__________________<\/strong>    <strong>____________________________________<\/strong>    <strong>_____________<\/strong>    <strong>___________<\/strong>    <strong>______________<\/strong>    <strong>_________________<\/strong>\n\n    {'TCGA-13-1510-01'}    {'TCGA-13-1510'}        {'IIIC'}        {'G3'}    {'No Macroscopic disease'}     {'Yes'}      {'Yes'}      {'Yes'}        {'Yes'}       {'Ovarian Cancer'}    {'High-Grade Serous Ovarian Cancer'}      {'HGSOC'}      {'Primary'}     {'Matched'}            2.9667     \n    {'TCGA-23-1809-01'}    {'TCGA-23-1809'}        {'IIC' }        {'G3'}    {'1-10 mm'               }     {'Yes'}      {'Yes'}      {'No' }        {'No' }       {'Ovarian Cancer'}    {'High-Grade Serous Ovarian Cancer'}      {'HGSOC'}      {'Primary'}     {'Matched'}               NaN     \n    {'TCGA-20-1687-01'}    {'TCGA-20-1687'}        {'IV'  }        {'G3'}    {'No Macroscopic disease'}     {'Yes'}      {'Yes'}      {'No' }        {'No' }       {'Ovarian Cancer'}    {'High-Grade Serous Ovarian Cancer'}      {'HGSOC'}      {'Primary'}     {'Matched'}               NaN     \n    {'TCGA-13-1404-01'}    {'TCGA-13-1404'}        {'IIIC'}        {'G3'}    {'1-10 mm'               }     {'Yes'}      {'Yes'}      {'Yes'}        {'Yes'}       {'Ovarian Cancer'}    {'High-Grade Serous Ovarian Cancer'}      {'HGSOC'}      {'Primary'}     {'Matched'}            1.4667     \n    {'TCGA-13-1506-01'}    {'TCGA-13-1506'}        {'IIIC'}        {'G3'}    {'NA'                    }     {'Yes'}      {'Yes'}      {'Yes'}        {'Yes'}       {'Ovarian Cancer'}    {'High-Grade Serous Ovarian Cancer'}      {'HGSOC'}      {'Primary'}     {'Matched'}           0.73333     \n    {'TCGA-23-1116-01'}    {'TCGA-23-1116'}        {'IIIC'}        {'G3'}    {'>20mm'                 }     {'Yes'}      {'Yes'}      {'Yes'}        {'Yes'}       {'Ovarian Cancer'}    {'High-Grade Serous Ovarian Cancer'}      {'HGSOC'}      {'Primary'}     {'Matched'}            1.5667     \n    {'TCGA-20-1686-01'}    {'TCGA-20-1686'}        {'IIIC'}        {'G3'}    {'No Macroscopic disease'}     {'Yes'}      {'Yes'}      {'No' }        {'No' }       {'Ovarian Cancer'}    {'High-Grade Serous Ovarian Cancer'}      {'HGSOC'}      {'Primary'}     {'Matched'}               NaN     \n    {'TCGA-13-1505-01'}    {'TCGA-13-1505'}        {'IIIC'}        {'G3'}    {'1-10 mm'               }     {'Yes'}      {'Yes'}      {'Yes'}        {'Yes'}       {'Ovarian Cancer'}    {'High-Grade Serous Ovarian Cancer'}      {'HGSOC'}      {'Primary'}     {'Matched'}               1.6     \n\n","truncated":false}}
%---
%[output:052a0b64]
%   data: {"dataType":"tabular","outputData":{"columnNames":["Entrez_Gene_Id","TCGA_01_0630_01_11","TCGA_01_0631_01_11","TCGA_01_0633_01_11","TCGA_01_0636_01_11","TCGA_01_0637_01_11"],"columns":6,"dataTypes":["double","double","double","double","double","double"],"header":"5×6 table","name":"ans","rows":5,"type":"table","value":[["2","0.9968","1.1713","1.9008","0.3202","0.9608"],["9","0.7437","0.9164","0.5575","0.6509","0.7654"],["10","0.1070","0.3039","0.0823","0.4172","-0.0750"],["12","0.7266","4.3839","3.0339","-0.6705","5.9939"],["13","-1.6075","-1.6908","-0.5149","-1.7002","-1.7817"]]}}
%---
%[output:6991c56d]
%   data: {"dataType":"matrix","outputData":{"columns":1,"header":"5×1 cell array","name":"ans","rows":5,"type":"cell","value":[["'TCGA-61-1895'"],["'TCGA-13-1510'"],["'TCGA-23-1809'"],["'TCGA-20-1687'"],["'TCGA-13-1404'"]]}}
%---
%[output:6b515ea6]
%   data: {"dataType":"matrix","outputData":{"columns":1,"header":"5×1 cell array","name":"ans","rows":5,"type":"cell","value":[["'TCGA-13-1510'"],["'TCGA-23-1809'"],["'TCGA-20-1687'"],["'TCGA-13-1404'"],["'TCGA-13-1506'"]]}}
%---
%[output:21cf3a32]
%   data: {"dataType":"matrix","outputData":{"columns":1,"header":"5×1 cell array","name":"ans","rows":5,"type":"cell","value":[["'TCGA-13-1510-01'"],["'TCGA-23-1809-01'"],["'TCGA-20-1687-01'"],["'TCGA-13-1404-01'"],["'TCGA-13-1506-01'"]]}}
%---
%[output:902b3b5b]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"491"}}
%---
%[output:42cb7785]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"491"}}
%---
%[output:4a0337b5]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"488"}}
%---
%[output:13c0eb10]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"488"}}
%---
%[output:408bc634]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"488"}}
%---
%[output:57306e8a]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"488"}}
%---
%[output:407125ac]
%   data: {"dataType":"matrix","outputData":{"columns":1,"header":"3×1 cell array","name":"patientOnly","rows":3,"type":"cell","value":[["'TCGA-25-2390'"],["'TCGA-29-2435'"],["'TCGA-61-1895'"]]}}
%---
%[output:5c48fc7c]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["488","7"]]}}
%---
%[output:79cd82ac]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["488","21"]]}}
%---
%[output:1a539e95]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"496"}}
%---
%[output:60a29bfb]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"496"}}
%---
%[output:33f8b243]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"488"}}
%---
%[output:360f3e3d]
%   data: {"dataType":"text","outputData":{"text":"\nans = \n\n  0×1 empty <a href=\"matlab:helpPopup('string')\" style=\"font-weight:bold\">string<\/a> array\n\n","truncated":false}}
%---
%[output:63df1b45]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["11838","488"]]}}
%---
%[output:029f1d89]
%   data: {"dataType":"text","outputData":{"text":"SAMPLE_ID: 0 missing\nPATIENT_ID: 0 missing\nTUMOR_STAGE_2009: 0 missing\nGRADE: 0 missing\nTUMOR_RESIDUAL_DISEASE: 0 missing\nACGH_DATA: 0 missing\nMRNA_DATA: 0 missing\nSEQUENCED: 0 missing\nCOMPLETE_DATA: 0 missing\nCANCER_TYPE: 0 missing\nCANCER_TYPE_DETAILED: 0 missing\nONCOTREE_CODE: 0 missing\nSAMPLE_TYPE: 0 missing\nSOMATIC_STATUS: 0 missing\nTMB_NONSYNONYMOUS: 172 missing\nPRIMARY_THERAPY_OUTCOME_SUCCESS: 94 missing\nPLATINUM_STATUS: 1 missing\nOS_STATUS: 0 missing\nOS_MONTHS: 1 missing\nDFS_STATUS: 25 missing\nDFS_MONTHS: 93 missing\n","truncated":false}}
%---
%[output:9af91628]
%   data: {"dataType":"matrix","outputData":{"columns":1,"header":"8×1 cell array","name":"ans","rows":8,"type":"cell","value":[["'IIA'"],["'IIB'"],["'IIC'"],["'IIIA'"],["'IIIB'"],["'IIIC'"],["'IV'"],["'NA'"]]}}
%---
%[output:03e136de]
%   data: {"dataType":"matrix","outputData":{"columns":1,"header":"3×1 cell array","name":"ans","rows":3,"type":"cell","value":[["'G2'"],["'G3'"],["'NA'"]]}}
%---
%[output:0c4f6c2d]
%   data: {"dataType":"matrix","outputData":{"columns":1,"header":"5×1 cell array","name":"ans","rows":5,"type":"cell","value":[["'1-10 mm'"],["'11-20 mm'"],["'>20mm'"],["'NA'"],["'No Macroscopic disease'"]]}}
%---
%[output:08333ba5]
%   data: {"dataType":"matrix","outputData":{"columns":1,"header":"5×1 cell array","name":"ans","rows":5,"type":"cell","value":[["''"],["'Complete Response'"],["'Partial Response'"],["'Progressive Disease'"],["'Stable Disease'"]]}}
%---
%[output:51fadab2]
%   data: {"dataType":"matrix","outputData":{"columns":1,"header":"5×1 cell array","name":"ans","rows":5,"type":"cell","value":[["''"],["'Missing'"],["'Resistant'"],["'Sensitive'"],["'Tooearly'"]]}}
%---
%[output:61c686a1]
%   data: {"dataType":"matrix","outputData":{"columns":1,"header":"3×1 cell array","name":"ans","rows":3,"type":"cell","value":[["''"],["'0:DiseaseFree'"],["'1:Recurred\/Progressed'"]]}}
%---
%[output:0fccea96]
%   data: {"dataType":"matrix","outputData":{"columns":1,"header":"2×1 cell array","name":"ans","rows":2,"type":"cell","value":[["'0:LIVING'"],["'1:DECEASED'"]]}}
%---
%[output:5534b29f]
%   data: {"dataType":"text","outputData":{"text":"  Value    Count   Percent\n   IIIC      353     72.34%\n   IIIB       20      4.10%\n     NA        5      1.02%\n     IV       79     16.19%\n    IIC       17      3.48%\n   IIIA        7      1.43%\n    IIB        4      0.82%\n    IIA        3      0.61%\n","truncated":false}}
%---
%[output:506858a2]
%   data: {"dataType":"text","outputData":{"text":"  Value    Count   Percent\n     G3      419     85.86%\n     G2       57     11.68%\n     NA       12      2.46%\n","truncated":false}}
%---
%[output:058ac4e1]
%   data: {"dataType":"text","outputData":{"text":"                   Value    Count   Percent\n                 1-10 mm      223     45.70%\n  No Macroscopic disease       89     18.24%\n                      NA       57     11.68%\n                   >20mm       89     18.24%\n                11-20 mm       30      6.15%\n","truncated":false}}
%---
%[output:65aadeb2]
%   data: {"dataType":"text","outputData":{"text":"                Value    Count   Percent\n    Complete Response      276     70.05%\n  Progressive Disease       37      9.39%\n       Stable Disease       24      6.09%\n     Partial Response       57     14.47%\n","truncated":false}}
%---
%[output:6c2fd3bc]
%   data: {"dataType":"text","outputData":{"text":"      Value    Count   Percent\n  Sensitive      197     40.45%\n    Missing      143     29.36%\n  Resistant       90     18.48%\n   Tooearly       57     11.70%\n","truncated":false}}
%---
%[output:90d9db0e]
%   data: {"dataType":"text","outputData":{"text":"                  Value    Count   Percent\n  1:Recurred\/Progressed      349     75.38%\n          0:DiseaseFree      114     24.62%\n","truncated":false}}
%---
%[output:53886646]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"395"}}
%---
%[output:38a61561]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"68"}}
%---
%[output:496debd0]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0"}}
%---
%[output:46b844fb]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"25"}}
%---
%[output:98520319]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["395","21"]]}}
%---
%[output:06fc23f2]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"ans","rows":1,"type":"double","value":[["11838","395"]]}}
%---
%[output:8dea4a8c]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"283"}}
%---
%[output:1a4860f9]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"112"}}
%---
%[output:1bcf3488]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAp0AAAGSCAYAAACha3p2AAAQAElEQVR4AeydC9xVU\/rHn91FpTtdldFFFFFuFULJkHIJ0VUlpCFShDIuMRpk+FeMoVANSZoQKkYoQxfRRVTSxagpXXSnlPr3Xdlv53077\/2c9+xz9q9P6z377L323mt91z57\/faznvXsQvv0TwREQAREQAREQAREQATiTKCQ6Z8IiIAIiECCCej0IiACIpD6BCQ6U7+NVUMREAEREAEREAERSDiBwIvOhBNSAURABERABERABERABPJNQKIz3wh1ABEQARFIeQKqoAiIgAjkm4BEZ74R6gAiIAIiIAIiIAIiIALZEZDozI5Qdtu1XQREQAREQAREQAREIFsCEp3ZIlKGZCfw1VdfWYMGDaxWrVpp6YknnjikWv\/973\/t7LPPTstD\/o4dO9r27dtd3n\/84x9p21h2K+P8h\/NQDhLLcT5dKA5\/5513prXjBx98kPA6RyvPunXr7LzzznPl5JPveS1o5PUfeT3n9Xg52S8R58xJuaLliUVZuUfAlt8p9xqOGe1cketoU9qWffjke+R2LYtAKhKQ6EzFVlWdsiXw5ZdfpolJP\/O3335ra9eu9b\/qUwREQAREQAREIIYEJDpjCDOYh1KpohH47rvvbM2aNek2ffrpp7Zv37506yK\/9OzZ05YvX+4Sy5Hb4rXMeQr6nPGqS1COi5XbZ3rBBRcEpVjpylGpUiWbNm2au9b45Hu6DLn4ctJJJ9n8+fPdscaMGWOlSpXKxd6pn1V8Ur+NVcPgEJDoDE5bqCQFQIAO96ijjrINGzbY119\/nXbGLVu22Lx586x48eJ25JFHpq3XggiIgAjEhIAOIgIiYBKdughCRcDzPGvSpImrc6RlE39OrF8Iztq1a7vtGf\/gU4n\/FYllf\/svv\/xio0aNsubNmxv7sr1u3brWrl07mz179iHW09WrV9s999xjp5xyivPZIz\/LrGObf1w+OQ\/bSSyzjhTpB\/j+++\/bxIkTrWnTpu54WG4eeugh27p1K1nTpc8\/\/9yuvvpqq1OnjtXdX0bOuXDhwqj+g\/il4Z\/GufFXW7p0qd16661uP+rZqlUr++yzzw6pH+d97LHH7IwzznDl4Tw33HCD4b6QrjD7v+SWHce48cYbjTpSLtJZZ51lnI\/z7j9ktv8j2fk+nfjT4VfH8fjkWoChfx7Ywvi3337L9vh+hoxl5Vi9e\/e2jG3s54\/8zFgevrOdT8qXm3JmbEf8DzkWifq8++67RltyTXBceA4fPtxoG\/L4ie+5uc79\/bL65JjdunVz10mjRo2Maywy\/5IlS9KuI\/KRn+185rQske39yiuvWNeuXd3vlOtz1qxZlhWfjG3o88nuetuxY4c9\/PDD7jrlt3LhhRfa5MmTD\/mtUJdoyW8X9mN\/Esu0Fdsi98kNi8j9tCwCiSAg0ZkI6jpnJIECX65Xr56zaGLp3LRpkzv\/ggULbNu2bVazZk37wx\/+4Nbl5M\/u3btd5zJw4ED7\/vvv0zqVX3\/91RB4nTp1stdeey3tUAgOOs9x48YZ1lV\/A8usu+aaaw7peP08mX1y7j59+tj\/\/vc\/l4UOb+TIkXbffffZnj173Dr+0Olde+219sUXXxgdF2XknAihzZs3kyXTBKu2bdsanR774YawePFiY\/gf\/1h\/R8pAp\/7cc8\/Zxo0b3Wryf\/jhh3bZZZe5jtet3P8nt+ywRMNz6tSpRh33H8L9xw+X83Xp0sXWr1\/v1uXnz48\/\/mjt27e3kfsZ+uehXjCeNGlSjg4N66uuusoiy8qx3n77bbvyyivtm2++ydFxssqU33LC\/4EHHrDbbrvNaEuuCc4Hz7\/+9a\/Wq1evNL9n8iKiuNaiXeewp67sn5tUokQJO\/\/8890ujD7MmDHDLft\/EIX+ddSmTRsjf37KMmTIEPvkk0+M65dRDUY9\/HNl\/MzueuMBLFLA+\/vv3LnT\/vSnP9lLL73krlPOhTsPPBHzfPfzRvukfoha2oX9yE9imXPyuyYP+\/KZVbvwe4m8\/7CPkggkkkChRJ5c5xaBRBA4+uij3RD6ypUrbdmyZU6Y+Z3dySef7KwuOS0XnTVCjPxYaugkOeYbb7xh+OHRkY\/cL158MYToYDv5ETFYckh0JJ7nOT9TrDF0MuTJSULUvfjii06s0qkedthhbrc5c+akCVHOP2zYMCMv2xEVWJUQCmXKlHGC2+2UyR8EeYsWLZwLAqKVWf5kpdP96KOPWHQd+YgRI5z\/IOf4v\/\/7P1emadOmGWw496BBg+yHH35w+XPDDh5wQYAULVrU8Muk\/DwsdNtvKeOALCP2WM5PopxYtBAniG0s1hyPMrz33nvueuF7ZokHi7\/97W9OcGA5p9OnrH7b0BbPPvusIRgyO0ZO1ue3nFh5X3\/9dddutCfXbmR9P\/74Y8OKTlnY9q9\/\/cs8z7O+ffsa1yzXMdZ32oOyIMhhRP7cJEYIqlWr5nbh4QTLHV+4tqZMmcKisf20005zy\/kpS6FChdyoAOUfPXq0Va9e3R0z459du3bZ008\/7R6csDLCCuv3vHnz3IgC+efOnevuHyxHJtr1iCOOcA9otPvf\/\/53K126tOP8wgsvuN9EZP6MyzNnzjTKBkseQrmuSbfccovjz\/Xkc8nNbyjjeQ79rjUiEH8CEp3xZ6wzBIgAAothxOOPP96wSCCgEAF0thTzxBNP5CPHCesMnSM7YA3DgoV1kWFpOg86KjqIihUrksVZQ93C\/j8M661atcqKFCli1113nevAyP\/ggw+6zmV\/lhz9v\/76693weOHChd0Q\/6mnnur2o1xYUPlCJ4tIYJnJM1jhyI9l9\/bbbzeEA9syS8ccc4z169fP4Fe+fHk3ROnnxTLGMnX597\/\/zaJxDoZsOQciH+snGxBk06dPZ9H51VJGvmTH7ueffzYse+SFL53\/Tz\/95CbF3H\/\/\/W6SDOywuJEnPwmBcO+99zqhg2UN1wLWcUzOyXXDcmYJCzdlYTttwzAuHBgeZZie9QinlfsfeljOa6JMeS0nDHlYQiBxnLvuusu4Rqlv9+7drUKFCk4kISQRYJTbv4YQ+fxuHn\/8cYM9x6AOXAdYc1nOTUL48UDDPgg5LHoscz5cP1hmO\/lYzk9ZLr74Yqtfv7673rn2Pc\/jkIekYsWKGQ9QtCO\/X8Tq+PHjnfXXf0Dl2uXekXFnz\/OM3xQjKrT7RRddZFdccYXLRn6uD\/clyh\/ahfMg4qtWrWo333yzu8bxRe\/QoYNhmUWMvvXWW0a75Pb+E+WUWiUCBUpAorNAcetkiSZA51GyZEnnl0hZsAbSySFo6GiPPfZYVqdLWX057rjj3JA8eRBd1+0Xj1hkEAP4g9FBsM1P5557ruvw+I71hM4Ugcawmz\/Uz7bcpMgyM2SIlSXj\/nTgvjig00Xo+nnoHKtUqeJ\/jfpJB4g48TeWLVvWuSj43\/lEdGCJZBmxAhsshiSsNKwn+UKC7XT8rMuOHW2GryF5Yfryyy8739zLL7\/cJkyYkDYMzPb8pnLlylnlypXTDkNd6fTTVmSzwExxP0tk20QKGXx9efjx8+XlMz\/lxP8VQcV5Mx6HclE+tiO8KDf5EEuIJx6oEOK4NOBjyrb8JM\/z3MMSDz5Y1LG+czwszQg71jdv3jzdg1hey0LZOXZOEudmmBt\/a1wAEOb\/+c9\/nGtKVvtzrdSoUSMti+d51rhx47Tv\/vWftiJiAdGOfzmriK7RrFkzN\/LCbwixzUMb21asWOF8tnPzG2I\/JRFINAGJzkS3gM6fEAK+RXPRokVuGAzrFUN4DIlHiqvsCsc+DKX64on8dFavvvqqtWzZ0vBDo4NkPQlr14ABAwwRxXcSHQjD3QxBY03EGsL6vCTEJMIz476+4My4Pq\/fOQdiIHJ\/rKpwjFwXbdkvS27ZYS1lQhLWI46L+MRazESRM888000mgj3bYpmoK1xzekyEU07zxjJfbsqJRQ2Bk9PzIwDxB0ZQ8+CGAMLFAb\/F3PxeMjsf1nkehtjOuXDBwC2D76xnO8sktue1LAhCjpFdYogfNwKENdcU94WbbrrJsDDym85u\/4zbI38r\/vWfMQ\/fsV5iTWc5q0T7cf3n9jeU1TGTYJuKmAIEJDpToBFVhdwTwAKFZRPrHJYyjoBwZPg40srF+uwS1hOsllhoEEBYDvFpZD9EEUOsDEfyHcGEeGLyDdY6f8jM8zxnQcFvjsk5sRYukZ0e5YhHQvT452GyDJayaAmx4p8\/N+wY+kWwMwT7\/PPP2yWXXOKGgTkWAgqBgIWZDpl1iUo+g0SdPyfn9TzPuXXkJC\/ihpn7DPlStyeffNL5HNLG0azqOTlmxjwIV6zWrGfkAV9SrPN8\/+Mf\/+h8IlkuiLJwHsqACwTLPAziNnL33Xe72eiIbtbnNVHXzPbl4cZ\/IMUtxXfRyfg7QpAjhDlObn5D5FcSgUQSkOhMJH2du8AJYBk4\/PDD3SQflunEfJHCsLjneXkqk+d5bpgdHyx85RCbWOU4GJaSSGsn6+i8GS5+5JFHjCE7EiF12EaHlyt\/P3bKJuHDyjnJNm\/evHSTWLD2Ir7Zlp\/ErH+EPMfAfxTLJ8skhCZDhP6kDNb5yfNyxw5rFS4JQ4cOdSGp4O13wNQtJ5Yi\/9zx+IwcSqU8kecYOHCgGy6NxiEyX7yX8cvl4YjzYFn\/\/vvvWXQJfvjj0l4Ma2OJ94d1ebCIHDqmnWP1gAQ3Jl5xPEYPsJrznTK4gu3\/w8NFQZQFdxvKsf+UbrKRLxS5pqkz6zNL7MeDkb+de4z\/0Mk6RCKf0RIPvXBnW8Z24R5ywgknuOuHERGOSz6S5+XuN8Q+SiKQCAISnYmgrnMmnAA3dyybfkHoTPGP8r\/n9JNZpnQSJIbjfL9MhCY+Wf5x8JlkyI5JGOQlduXYsWOd+KPzwBfSD1uEqMKP0N83Fp+IToQOx2JWMtZdZtYjOIlHmdWQH\/vkJDHJ4fTTT3dZmW3LjHrqjA8bkyPYQDn8Tjc37OiAcVeAHRNzmDFP+Ulw5jwcH8ub74PI90Qkykc5OTczjbFKUU5EwzvvvMNqi+TgVhTwHyxqrVu3dv7FiDtEHmKO64C3FvlWRgQfvxMe0CgigorIDPCmXYmIwPpYJK5PHsQ4FmXiEyHKepZJWAELoiyMdvhCkwdC6spvGh9Xrm3KklWCi9\/uTJzzr3\/KzsNtZvvSLkw88jzPTXT02wUfXCz5cGEUhfBjnuc5izPXGim7+09m59T62BHQkbInINGZPSPlSEECnudZ5M0fiwodAlVFgPpWQb5nlZhYQbBu8rz55pvumHQAHBsLHOvxf2PCEMPDhEliljBDlQwVIz7oVOlE8GPzPM\/NDGdojX1jlTgnMf7osDh3\/\/79XYB4hAfr\/A42P+ejw2SyBfVBSNPx4jtL3Zhwwjmw9FEWzpMbduxDPFFEBwIdlwUmvJCwKCOGqEePHj0s1oKdsuYmcR3dcccdzm+XsjK5jHISeonv1IG2oE65OW6s82It5kUBnuc5i\/E5RPUX8gAAEABJREFU55zjxDDD57QfD0Zw9jzPxViFL2UYOXKk+e3KBDD\/t4IQRVyTJy+J64che\/94nucZAoz1\/vE8r2DKgvsNgpfzcu1yDRNO7ZlnnrHIhxqEKHkiE+Vn9MRvdz5pd1xrWM7ut83DFdcKx2RCF+3SsGFDI5yU53mGbyn3FLbn5jdEfiURSDQBic5Et4DOnzACWDYRmBQA8YfwZBnR4q\/ne1YJIUUcPoaPGa70O2Y6GCxEBN\/GQkE+joOVj1idCCUsg553YDgfIYLvGBOQeOOO5x1Yzz6xSoSLQTAgJjzPc6FY6MAof7ly5dxp8mspRHDhl8pxse5y0MMOO8zN3iUmJFZA1pFgwrlzyo7yM5GDcE\/+MD7Hob0Q9fgdkod1iU6UAw6Ui7alPHzynfVsZ10iE+KIhwDcFLgmuGYpD+1G+2GJpz1Zh\/CJvHa4zpnRzbXMTHbyMOxM3EiW85oYPuZ3w\/4I9Sa\/vz2M734qiLLwgIj4RiTSbpyb3yu\/Z16o4F9\/+HFjHWa7n7h3IE67devmHjxYT53wQ+Z4fM8q0S4EfOehzW8Xz\/OMsGVPPfWUC+bveQfuD7n9DWV1Xm0TgYIgEGrRydM8N4LOnTtHDbnCkzuhbBjywXrFDYCJH4TCKYjG0TliQwBfSWbd4ozP0CHD1xwZAURcTdYzbOZbMDLLzwQf8pJY5hgkOigsNFg26XTZTlBoOiQmDbGdfH7C\/xArJ8N2dNTkxweUTh7h6XkHOhTycx62k1hmHQmhxjoSFivW+cnfRp2pi7+eTzpxwhlxXoYJmRxBvffu3ctm98YXX3ywL8fgHJHcyBi5jfOxzk+4LnBcXpHJvjBBTCPy\/Tz+J2xyw47f4eDBg51ljmOT8J\/jdxrt+P55Ij8pL\/uRfHa0CcOhrOOT7\/4+LLOObRk5+HmifVIeykXbsi+ffGd9ZP7clCeyLJSJ7\/6xWGYd54osZ2RbRa5nP9oaazfXBNcs+9JutB\/tSB4\/RV47tCm\/GeqCEGM\/riffOpjVOf3jRfvE+ktcTI7HJ9+j5ctNWaLxjTxmZmXlPsFLG2g3ysPvld8zs+mxQLIOYYhIJC9sWcdvhj6DGKb+vtwLcFXwvIO\/7cj2ot347pcrY7vwe8WlBIsr2\/x8fOb2N8Q+SiElEIBqh1Z08nTK0ygdmN\/hRrYHwyZ08gRAJlwHT728xowbM0NS+GdF5teyCASZAOKaIX9EG8Pf+BoyFMrvAP9OgrNTfiy+dKAsK4mACIiACIhALAmETnRi3cQywuvFEJJ0vNGAMtuQcBU4Z\/M0S2w2RCdDhAxpIFixhEbbV+tEIGgEGN7DOkS5+A34Pp2ITH4HrEOMMuxPnpAlVVcEREAERKAACIROdDILEV8bhoII1M3QSjTOzN7Ezw2Hbc87OCRSo0YNYwiUMCJMXoi2r9aJQNAIMATHMCNDoVg78cnzy4gPH8G28dnzffj8bfoUAREQAREQgVgRyFp0xuosATqO53lGWBdCWDB07juJZywis2AJDpxRlDLzl\/AVGfPruwgEnQDCE580Atnjk4f\/GQkfPiaUSHAGvQVVPhEQARFIbgKhE504gb\/wwgvG+3Tz8mYJfOPw5ySEBTN9k7v5VXoREIFkIKAyioAIiEAqEAid6MxPoxE8mTiEDE0S+iIyflx+jqt9RUAEREAEREAERCDVCSS56Cy45uGNKAR8JoA3b3DB0pnd2QmczGQkpZkmBmKga0DXgK4BXQNhvgbQBNnphlTfLtGZgxZm0hC+cPPmzTNCLPHGiOx24+Li\/bgETlbqaGIgBroGUvga6Ki66frWNZDdNYAmQBtkpx9SebtEZzatS+ik9u3bG5ZOgiFfeuml5nkHZ7NntjsXFhOREKkEDVYaYwXNgNcm0j5qg7yzF8O8s+N6F7\/88RPD\/PMTw9gwhGN+EvcCNAHagH4prEmiM4uWx7JJgHje2PLSSy8Z4ZOyyB5tkzEjmPiISk2soBn4b0dRG+SdvRjmnR3Xu\/jlj58Y5p+fGMaGIRzzk\/x7QVShEKKVEp2ZNDaThjCF8yo4Zrsz6z2TrFotAiIgAiIQaAIqnAiIQBAISHRGaQXezoLQXLZsmRHbECvngAEDLDI99thjpnidUeBplQiIgAiIgAiIgAhEISDRGQXKjh07jHicbOJz7NixljFNmjTJdu7cSZakTqlceN60wwQwPlO5nvGsG+zEMO+ExS\/v7Pw9xdAnkfdPMcw7O+0ZWwKhFp2lSpVKm9zCso+WZUQmb2vJLE2bNs0qVark76LPABLgRtulSxfjM4DFS4oiwU4M895U4pd3dv6eIWHoVzcun2IYF6w6aB4IhFp05oGXdhEBERABERABERABEcgDAYnOPEDTLgVIQKcSAREQAREQARFICQISnSnRjKqECIiACIiACMSPgI4sArEgINEZC4o6hgiIgAiIgAiIgAiIQJYEJDqzxKONIpAdAW0XAREQAREQARHICQGJzpxQUh4REAEREAEREIHgElDJkoKARGdSNJMKKQIiIAIiIAIiIALJTUCiM7nbT6UXgewIaLsIiIAIiIAIBIKARGcgmkGFEAEREAEREAERSF0CqhkEJDqhoCQCIiACIiACIiACIhBXAhKdccWrg4uACGRHQNtFQAREQATCQUCiMxztrFqKgAiIgAiIgAiIQGYECmS9RGeBYNZJREAEREAEREAERCDcBCQ6w93+qr0IiEB2BLRdBERABEQgJgQkOmOCUQcRAREQAREQAREQARHIikB+RGdWx9U2ERABERABERABERABEUgjINGZhkILIiACIpCMBFRmERABEUgOAhKdydFOKqUIiIAIiEASE1i1apXNnDkzIWnOnDk2f\/78hJw7ss4wSOImVNFjQCClRWcM+OgQIiACIiACIpAvAoitfv36WceOHROSunTpYn379rXu3bsn5Px+vWEAi3zB1M5JTUCiM6mbT4UXAREQgcATCH0BEVqzZs2yn0\/tbtub3hXKtLPuZQYDWIT+gggxAInOEDe+qi4CIiACIlBwBPYeXsH2VDg+pKluwYHWmQJLQKIzkU2jc4uACIiACIiACIhASAhIdIakoVVNERABERCB6AS0VgREoGAISHQWDGedRQREQAREQAREQARCTUCiM9TNn13ltV0EREAEREAEREAEYkNAojM2HHUUERABERABEYgPAR1VBFKEgERnijSkqiECIiACIiACIiACQSYg0Rnk1lHZsiOg7SIgAiIgAiIgAklCQKIzSRpKxRQBERABERCBYBJQqUQgZwQkOnPGSblEQAREQAREQAREQATyQUCiMx\/wtKsIZEdA20VABERABERABA4QCLXo3Ldvnz3\/\/PPWuXNn2759+wEiEX\/ZPnnyZGvRooXVqlXLTjrpJHvooYds69atEbm0KAIiIAIiIAIiEGACKlpACIRWdO7evdueeeYZGzx4sO3duzdqc0yZMsX69OljZcqUsSeffNLatm1rr7\/+uvXq1cu2bdsWdR+tFAEREAEREAEREAEROJRA6EQn1su5c+faNddc44Tkb7\/9diiV\/WvWr19vQ4YMsYYNG9ro0aOtTZs29sADD9igQYNs9uzZNnHixP259F8EkpyAii8CIiACIiACBUQgdKJz4cKF1q1bN1uwYIFdeOGFbsg8Guv58+fb0qVLrX379la6dOm0LOeee67Vq1fPPvroI9u1a1faei2IgAiIgAiIgAiIQF4IhGWf0IlOz\/Ps9NNPt\/Hjx9vjjz9uJUuWjNrWixcvtlKlSlnt2rXTbWeovWbNmk6QbtmyJd02fREBERABERABERABEYhOIHSis379+vbCCy\/YKaecYoUKZV79lStXOl\/OI488Mh05z\/Pcfkw82rBhg+mfCIhAPAno2CIgAiIgAqlCIHPVlSo1zEc9EKVFihQ55AjHHnus7dmzxzLzB43cYdasWTZz5sy0hJjduXOnKRUMA1wgfvnlF\/HOxzUnhvm7VsUvf\/y4VyY7Q\/qLyH4hzMu0ZxgSfX1k37969erkbvYYlV6iM0YgMzsMk5E6duxofho1apStWbNGqQAYrF271pgQ9uOPP4p3HnmLYf5+q+KXP37cK1OBoUbFDvaQmzdvDsX9mL7e7\/f57Nev30EIIV6S6MxD43\/33XeGBbRw4cLZ7k1IpjFjxpifunbtalWrVlUqAAaVK1e2ihUrWpUqVcQ7j7zFMNvfapbXlvjljx\/3ylRgWK5cuWz7irBkqFChQpa\/Gdo8FRJ9vd\/v89m7d++wNHGW9ZTozARPjRo1jCcyrGSRWQi5RFxPJhnx44ncFm25WrVq1qRJk7TEcYsXL25K8WdQokQJK1asmFjn43oTw\/xdp+KXP37cJ1OFYbT+IYzrMNjQrqme6Osj+\/7GjRuHsbkPqXP8ROchp0quFXXr1nVvKVq2bFm6gvM2ohUrVlidOnWsbNmy6bbpiwiIgAiIgAiIgAiIQHQCEp3RuViDBg2csJwwYUK6tw9Nnz7dFi1aZM2bN3dWtEx212oREAERCAQBFUIEREAEgkIgkKKTcESTJk2ybt26udBGtWrVMj81atTIrWc7+eIFEl9Azj9jxgzr3LmzvfbaazZw4EAbMGCAUYbLLrssXqfWcUVABERABERABEQg5QgESnSuW7fOvWqS4O283xyrIgHY8Z2sXr268cksQNaznXwIwU2bNuWhYbLfpV27djZ06FA3zN6\/f38XUB6xyYz0yLcUZX8k5RABERABERABERCBcBMIhOgk3iXhBVq0aOEsik2bNrXnn3\/eiHGJT+Xs2bMNockn31n\/9NNPO4vjyy+\/bLyacvjw4bZ79+5ctSaTgZhVRmI5486e59nFF19sU6dOteXLl9tXX31ljzzyiJUvXz5jVn0XAREQARHICwHtIwIiEBoCCRedDJH36NHDnnjiCbvhhhsMYTlixAi74IILXLgbz\/PSNYbneW59q1atbPTo0S4\/+2GR7Nmzp7NKpttBX0RABERABERABERABBJOIOGik7fFnHHGGc6SSRwr3m2eGypYHdnv008\/NY7D8XKzf4DzqmgiIAIiIAIiIAIikDIEEi46mbCDhRLxmB+qiFWOw\/HycxztKwIiIAIiIAIHCWhJBEQgVgQSLjpjVREdRwREQAREQAREQAREILgEkkJ04vf5wQcfuJBFzFZnmXXBxRqOkqmWIiACIiACIiACIpBTAoEXnfPnzzdmtTPZiBnuJJaZSLRw4cKc1lP5REAEREAERCAVCahOIpA0BAItOnft2mXDhg2zmjVrurBFS5cuNRIhjI466igjbBJ5koa2CioCIiACIiACIiACISWQcNFJbE3eZx6NP4HhEZmEREJ4Fi5c2Egss47XUZIn2r5aJwImBCIgAiIgAiIgAoEhUCjRJeFtQu3btzfe8pNRfBYrVsyOOOIIFyQeceqXlWUCxLONPP56fYqACIiACIiACASLgEojAj6BhIvOI4880gh19MTKtFIAABAASURBVM9\/\/tPOP\/9890YiP9Zm2bJljVdRjhw50r19iDcPkRo1amTkx7eTPH5l9CkCIiACIiACIiACIhBMAgkXnQyX8z7zjz\/+2K6\/\/nr7y1\/+Ypdccom9\/\/77xusxEZ3jxo2z008\/PY1ggwYN3NuIWrZsmbZOCyKQfARUYhEQAREQAREID4GEi04fdalSpexPf\/qTezNR06ZN7ZZbbrG2bdvavHnzrGHDhjZ8+HC3bfr06eZbPj0v\/Ssy\/WPpUwREQAREQAREQARyRECZCoxAYESnX2PeTEQszg8\/\/NCqVq3qhtdvuukmW758uZ9FnyIgAiIgAiIgAiIgAklGIHCi0+d39NFH2zPPPGNvvfWWrV+\/3i666CL785\/\/bOvWrfOz6FMERCC+BHR0ERABERABEYgZgUCIzn379tncuXPtxhtvNCYKkVhesGCB1atXz8aPH+8E6IwZM1yg+Ggz3WNGRAcSAREQAREQAREQgcAQSJ2CBEJ0TpkyxTp06GBz5sxJI8syk4gmT57sYnNeeOGF9s477zhrJzPXmemOGCV8UtpOWhABERABERABERABEQgkgYSLzm3bttmIESPs2muvtdmzZ6dNFmKZdWwjD\/RKlCjhfDzx92Tbiy++aMT5ZJuSCIhA+AioxiIgAiIgAslDIOGik5icGzZssMaNG1vRokXTyLHMOraRJ23D\/oUyZcpY7969bezYscbEo\/2r9F8EREAEREAEREAERKDgCeT4jAkXnYcffrhVq1bNWTtXrFjhYnPi48kyVs6aNWta6dKlo1YI8Yk4jbpRK0VABERABERABERABAJDIOGik\/icd999tyEyW7RoYXXq1LHatWu7CUOs69OnjzGsHhhiKogIiIAI5JSA8omACIiACKQRKJS2lMAF3jA0depUe\/rpp61z587WtWtXe\/755411bEtg0XRqERABERABERABERCBGBBIlOg8pOhYPFu1amUPPfSQPfDAA3bBBRcY6w7JqBUiIAIiIAIiIAIiIAJJRyDhopOQR1u3bo0JOI7D8WJyMB1EBERABFKegCooAiIgAgVHIOGik5BH7du3t0GDBhmiMS9VZz\/25zgcLy\/H0D4iIAIiIAIiIAIiIALxI5Bw0VmpUiW7\/\/77jQDxTZs2deJz9erV2dbYn+HOe9rPPvts97pMXpPJ8bLdOQcZlEUEREAEREAEREAERCB2BBIuOqlKkyZN3NuGrrrqKhs9erSdc845dsopp9hNN93kJhS9++675icmGLH+1FNPdTPcX331VWvbtq39+9\/\/trPOOovDKYmACIiACKQGAdVCBEQghQgUCkpdiLnJBKIZM2Y4sVmkSBEnJB999FG79dZb0xLfEZjFihWze+65x8jPfuwflLqoHCIgAiIgAiIgAiIgAukJBEZ0+sXiDUPE7eQ1mJ999pmzdCIu8dd88MEH3XfWIzZ79OgR3jcS+cD0KQIiIAIiIAIiIAJJQCBwotNn5nmeValSxYVOQlwyUahLly7uO+s9z\/Oz6lMEREAEREAEEkJAJxUBEcg5gcCKzpxXQTlFQAREQAREQAREQASCTkCiM4sW+u2332zixInGrPpatWrZSSed5PxI169fn8Ve2nSAgP6KgAiIgAiIgAiIwEECEp0HWRyy9OKLL1qfPn2sYsWK9uSTT7pZ8syi5zWdEp6H4NIKERABERCBoBFQeUQgQAQkOjNpjDVr1tiYMWPsjDPOcGGc2rRp417PSVzQJUuW2Pvvv5\/JnlotAiIgAiIgAiIgAiKQkYBEZ0Yiv3\/fsGGD\/fTTT25ovXTp0r+vNSOm6FFHHWVz585NW6eFpCSgQouACIiACIiACBQgAYnOTGAXL17ciAX6\/fffG28\/8rNt3LjRva4zUoj62\/QpAiIgAiIgAiKQGwLKGyYCgRWdTOLZsmVLmuDbvXu3m9QzYMAAGzZsmK1bty6u7VSzZk278MIL7e2337Zx48YZ5Vm9erU98sgjdthhhxlvT4prAXRwERABERABERABEUghAoWCWBcEZceOHe1Pf\/qT7dixwxVx5MiRblLP2LFj7amnnrJOnToZItBtjMMf3oh07733WocOHax\/\/\/5Wp04d93rOlStXugD19evXj8NZdUgROEhASyIgAiIgAiKQSgQCKTqZwDNnzhzj1ZZ79+51vpVvvPGGHXHEEU7wPfzww\/bDDz\/YiBEj4tYW27ZtM96MRFlatWrlZq\/37dvXna9nz542f\/58y8m\/WbNm2cyZM9MSonXnzp2mVDAMdu3aZb\/88ot45+OaE8P8Xavilz9+3CuTneGePXty0l2EIg\/tmWQpT\/0HfX1k3x9PI1kyXTiBE53bt293Aq1hw4Y2ePBgJzy\/+eYbW7ZsmTVu3NiaNWtm7dq1c28mmjdvnjEEHw\/gxOecMmWK3XHHHW44n9nrvXr1sn\/961+uTFhBN23alO2phwwZYlht\/TRq1ChjZrzSmrhzWLt2rRHa6scff4z7uVK1PcUwf9ep+OWPH7+rVGDIxNRsO4uQZNi8eXMo7sf09X6\/z2e\/fv1C0sJZVzNwovPnn392FyQ+lf5kHcQlPp1nnnmmMexNYqIPs8t5As66innbyuz0ypUr2yWXXGKed\/CVm9WqVbPLL7\/cWVpXrVqV7cERzlhL\/USMz6pVq5pS\/BnQfsRY5bWp4p033lky1HWc7e9Y\/PJ23UX+XlOBYbly5bLtK8KSoUKFCtn+biLbP1mX6ev9fp\/P3r17h6WJs6xn4ERnxtIiKhGdiMzjjjvObWY2OcPuiE\/POygI3cYY\/uEc+R0WQaQSZslPNWrUMOqiVDzuHEqUKOEiEIh13lmLYd7Zcd2JX\/74pRLDGHZNSX0o+m3aNdUTfb3f7\/PJSG28Gi6Zjhs40VmyZEk75phj3HA6w9crV650\/pPHH3+8+aJz+fLlNnv2bMMait9nPICffPLJxrAOQ+yIXP8c+GW89dZb7kkNC5q\/Xp8iIAIiIAIiIAIiIAKZEwik6LzoootswYIFxmf79u2N2JjnnXeeFS5c2IUsuuKKK5yvHmGLiKWZefXyvoVznH322c6vtHv37vbmm2\/a3\/\/+d+efySQmTOUM3eb9DNpTBEQgvgR0dBEQAREQgSARCJzoBM4111xjt99+uwvCzsQiJg7dcMMNbLKvv\/7axe5kIk\/Lli3dunj8wZ8UkXnbbbcZr71k5vrQoUONtxG99tprdvHFF8fjtDqmCIiACIiACIiACKQOgYiaFIpYDsxi0aJF7dZbb7XFixfb0qVL7a9\/\/auVKlXKGHonODvhlLp06ZJugk88Cs85Kcdnn31mDOlTnldffdUaNGgQj9PpmCIgAiIgAiIgAiKQsgQCKTozo+15nvPjjNeQembn1XoREAERiAMBHVIEREAEQkUgsKKT105++OGHLiZn3bp1DZ\/OFStWOAso8a8I+B2qllJlRUAEREAEREAERCCJCQRSdO7evdseeOABu\/HGG+3zzz+3X3\/91fjHWwwIEj9w4EC7+eabjbcGsV5JBERABERABERABEQg2AQCKTrHjRtn+E4ShP2LL76wK6+80lE88sgj7bnnnnNvI5o2bZq98sorbr3+iIAIiIAIxJ6AjigCIiACsSQQONG5Y8cOe++994w4mffdd5+VL18+XX2PPvpoe+yxx9xknunTpxuz29Nl0BcREAEREAEREAEREIHAEQik6Pz++++tdu3ahwhOnx5ClO0Eaue1mf76gvvUmURABERABERABERABHJDIHCiMzeFV14REAEREIEQE1DVRUAEkopA4ERn2bJlrU6dOjZv3jzDkhmN5qpVq2zu3LkuH\/mj5dE6ERABERABERABERCB4BAInOgkBmerVq2M8Ej9+\/c3Xjnp4yKM0qJFi+yWW24x3slOPvL72\/WZRkALIiACIiACIiACIhAoAoETndC59NJLrUOHDvbpp5+6+JwTJkyw\/\/3vf9asWTNr3bq1ffXVVy5+J\/nIryQCIiACIiACwSOgEomACEQSCKTo5DWYAwcOtOHDh9tpp51mhQsXNqyce\/fuNQLFDxs2zB5++GEjX2RltCwCIiACIiACIiACIhBMAoEUnaBCaJ5\/\/vn2+uuvu\/ev8+5z3sM+adIkZ+1kO\/mUkpOASi0CIiACIiACIhAuAoEVneFqBtVWBERABERABAqcgE4oAgVKIHCic926dc6Ps1atWpZdOu+884z8BUpMJxMBERABERABERABEcg1gcCJTs\/zrGrVqla9evVDkh8eqUiRInbuuecas9eLFy+e60prBxHIloAyiIAIiIAIiIAIxJRA4ERnxYoVbezYscYrLjMmYnMuWLDAOnfubGvXrnWfZcqUiSkQHUwEREAEREAERCAYBFSK1CIQONGZHd5SpUpZnz59rGTJkvb444\/bnj17sttF20VABERABERABERABBJMIOlEJ7xKly5tNWvWtIULF9qmTZtYpSQCISOg6oqACIiACIhAchFIStGJdXPnzp3Oyrlv377kIq7SioAIiIAIiIAIpAYB1SJXBJJOdO7evdvGjx9vH3zwgdWpU8f8yUW5qrUyi4AIiIAIiIAIiIAIFCiBwIlOQiARCimzcEnHH3+8DRgwwEFq3769FStWzC3rjwiIQKAIqDAiIAIiIAIikI5A4ESn52UeMskPo0S4pNGjR1uLFi3SVUZfREAEREAEREAEREAEfALB+gyc6MwqZJIfQmnkyJHWqFEj8zwvWDRVGhEQAREQAREQAREQgagEAic6o5ZSK0VABEQgxgR0OBEQAREQgYIlkHDR+dtvv9n69evd6yzx58xNYj\/2L1hkOpsIiIAIiIAIiIAIiEBuCUQRnbk9RP7yb9y40dq2bWtNmjTJdWI\/9s9fCbS3CIiACIiACIiACIhAvAkkXHQWL17cWrVqZcxEz21iP\/aPNyQdXwREQAQKnIBOKAIiIAIpRiDhorNMmTJ2991326BBg3Kd2I\/9U6xNVB0REAEREAEREAERSDkCCRedeSCatgv+nHojURoOLYiACIiACIiACIhAYAkEVnSuWLHCnnrqKRcInmDwGVOXLl3skksucZOQAktXBRMBERCBlCWgiomACIhA7ggEUnT+5z\/\/scsuu8yGDRtmY8eOjZrI43meYnXmrr2VWwREQAREQAREQAQSQiBwonPPnj1OZPKO9dtuu83effddO+mkk6xly5Y2bdo0e\/jhh61q1apupvu4ceOMYPLxJPftt9+6SU68553EZCfWZXVObRMBERABERABERABEUhPIHCic+vWrbZ8+XInKm+66SarV6+enXDCCbZy5UorWbKkderUyZ5++mlbsmSJvf\/+++lrE+NvH3zwgV111VW2Zs0a+\/Of\/2z33HOPfffdd64MCxcujPHZdDgREAEREIEYEtChREAEAkYgcKITS+eOHTusQoUKVqJECYfr2GOPtQ0bNpgfkxPL55lnnmmTJk2yXbt2uTyx\/rNp0yaAplGrAAAQAElEQVR75plnrEqVKjZmzBjr2rWrXX\/99W4d5xw\/frxpElOsqet4IiACIiACIiACqUogcKLT8zwrUqRIOt7HHHOMbd++3TZv3uzWs7148eK2dOlS27Jli1sX6z9ff\/21LVq0yBjir1atWtrhTzzxRKtfv76zeG7bti1tfVItqLAiIAIiIAIiIAIiUMAEAic6ibtZs2ZNY\/ia11zCo3Llyla0aFHzfSmxMO7du5dNcUvz5s2zI4880ho0aJDuHKVKlXKWz5dfftkoa7qN+iICIiACIiACOSSgbCIQNgKBE53FihWzyy+\/3FkxmbTDRCKG2rE2jhw50r755hubOnWqffLJJ8a6ww8\/PC5thg9p2bJljaF+fDkZ0q9du7ZdccUV9vnnn8flnDqoCIiACIiACIiACKQqgcCJTkAzU71Hjx723\/\/+10aMGGEISyYQMcGI2Jxs++mnn6x169aG5ZF9YpkQmuvWrTMmNeHHyTA7E4n69OljP\/zwg3Xv3t0I2WT6FycCOqwIiIAIiIAIiECqEQik6GQonVdcfvzxx9a3b1\/DfxOr59ChQw2fShLLrItHgzB8z4Sm\/\/3vf3b22Wcbk4batWtnvXr1sjfffNMqVarkAtfnxKdz1qxZNnPmzLSEBXXnzp2mVDAMmPT1yy+\/iHc+rjkxzN+1Kn7548e9MtkZ0p\/Eo6+K+zHjcALaMwyJvj6y71+9enUcaCbfIQMnOnm1JTE6Qcnw+TnnnGMMuRcuXNhZNt9++20jYeVkHflinRC5RxxxhBO711xzjfMn9c9RvXp1O\/fcc10IJ6ye\/vrMPocMGWIdO3ZMS6NGjXIhmAjDpLQmrizWrl3r3lj1448\/xvU8qdyOYpi\/a1T88seP31YqMCT6SmZ9RNjWb968ORT3Y\/r6yL6\/X79+YWvqqPUNnOgkLBLD69ddd53znUSERi15HFcyO7506dJObCJAM56KiU08ueakbIMHD3YTjwi7RCL0EsHtlaq6IP\/x5EA78fIAwl7F+DxxL3tQyiuG+btOxS9\/\/PgdpALDcuXKZexGQvudORq0a6on+nr6fD\/17t07tG0eWfHAiU6G1rkop02bZgxpn3LKKTZo0CAraNM0YZEI04SVLBIYQ+\/Lli1zvqRMNIrcFm0Za22TJk1csHs+a9So4SyoiFml4nFlQZxXrOTinHfOYph3dlx34pc\/fqnEMFr\/EMZ1GHVo19ik\/F9f8SoHfT19vp8aN24cxuY+pM6BE53ly5e31157zc1Ov+WWW9wkIiYTMczeokUL51+JGDykJjFewRD6UUcdZaNHj7ZI300mM82YMcOYzc4TeIxPq8OJgAiIgAiIgAiIQEoSCJzo9CljIbzjjjvs008\/NcIm8TpKZpTfddddhvUz3sPvRx99tPEaztmzZ1vnzp2dEOb1m1hfcWrv2bOn8zX1y6tPEQgSAZVFBERABERABIJGILCi0wfFZKF69eoZvpFffvml4ZzLKzAZfr\/zzjvTXo3p54\/lJ2GaXnzxRXdIQiYhOhs2bGivvPKK8ek26I8IiIAIiIAIiIAIHEpAazIQKJThe2C\/4tM5fPhw69+\/v4uR6Xme1fjdPzJehfY8z8466yx76623XLD6xYsXu7ihxx13XLxOqeOKgAiIgAiIgAiIQEoSCLTo3LRpk7MqXnjhhS5M0RNPPGHMGsfXc\/r06c7qqVdRpuR1qUqlOgHVTwREQAREIHQEAic6idE5efJku\/rqq61Ro0Z23333uTcTtWrVyvlV4uOJryc+n6FrLVVYBERABERABERABGJEoKAPEzjRiXXz0UcfNfw3TzjhBMO6OWfOHBs2bJidccYZho9nQUPS+URABERABERABERABPJHIHCikzidzBZn+BxfyiuvvNLFxMxfNbW3CIiACOSGgPKKgAiIgAjEmkDgRCdxOm+88UbT8Hmsm1rHEwEREAEREAEREIHEEci16ExcUXVmERABERABERABERCBZCUg0ZmsLadyi4AIhJmA6i4CIiACSUdAojPpmkwFFgEREAEREAEREIHkI5B6ojP52kAlFgEREAEREAEREIGUJyDRmfJNrAqKgAiIQMET0BlFQAREICMBic6MRPRdBERABERABERABEQg5gQCJzq3bt1qjz32mA0YMCBH6amnnrJvvvnG9u3bF3M48TmgjioCIiACIiACIiAC4SMQONG5c+dOmzRpko0dOzZHiTcVXXrppe51mbxCM3xNqBqLgAiIgAjkmoB2EAERKHACgROdlSpVsvHjx1ujRo2sevXq9swzz9i8efNs+fLl7pPvNWvWtAYNGti\/\/\/1ve\/XVV+3kk092n2+\/\/XaBA9QJRUAEREAEREAEREAEsicQONHJMPmzzz5rq1atslGjRtnFF19sZcqUcTXhk++sZxj+lVdeceL0+eeftzp16jgL6a5du1xe\/ckzAe0oAiIgAiIgAiIgAjEnEDjRuWnTJps5c6YTkzVq1IhaYSygp5xyistH\/ooVK1r9+vVt6dKltmXLlqj7aKUIiIAIiIAIJA8BlVQEUo9A4ETnnj17bMeOHVaoUCHzPC9L4uQjf5aZtFEEREAEREAEREAERCDhBAolvAQZCnD44YdbtWrVbO7cuW6IPcNm95Whd7aTj\/wITyYgHXHEEVasWDGXR39Sl4BqJgIiIAIiIAIikHwEAic6S5UqZa1bt7aVK1faLbfcYp9\/\/rn5s9L55Dvr2U6+kiVLuglFn3zyiR177LFp\/p\/J1xQqsQiIgAiIgAgkDQEVVARyTSBwopMaXHPNNXb99dfbwoULrV27dnb88cdbrVq13CffWc928jHEPmrUKDvssMOsU6dO2Q7Jc3wlERABERABERABERCBgiUQSNFZtGhR69+\/v73zzjt2ySWXGGGUPM+zChUquO+EVGI7+fbu3WsXXnihjRs3zho2bFiw9HQ2EYhGQOtEQAREQAREQAQOIRBI0UkpPc+zevXq2dChQ90s9WXLltns2bPdd2aue96BSUaEUerevbvVrFmT3ZREQAREQAREQAREwIQgeAQCKzpBtX37dlu3bl2maf369fbbb7+RVUkEREAEREAEREAERCDABAIpOrdt22a33Xabe+tQkyZNLLPUtm1b27hxY4DxqmgiEEQCKpMIiIAIiIAIFDyBQIrOF154wflzEg7pj3\/8o7Vv3z5qatWqlRUvXrzgqemMIiACIiACIiACIpAfAiHcN3CikzcKTZs2zfloTp482Z577jkbNGhQ1HT33XcrRFIIL1pVWQREQAREQAREIPkIBE507tq1y3766SdjshCvu0w+pCqxCIhAPglodxEQAREQgRQkEDjRyRuFeLMQoZD27duXgshVJREQAREQAREQAREIOoHYly9worNs2bJ23nnnufBIvHUo9lXWEUVABERABERABERABAqaQOBEJwB42xCB3m+44QYjEPzq1aujhk1SyCRoKYmACBQ0AZ1PBERABEQg9wQCJzqJy8lbiCZNmmQrVqywu+66y84555yoYZMUMin3Da49REAEREAEREAERCARBGIsOvNfBUIgEQopszBJkevJR\/78nzVnR5gwYYLVrl3bPvjgg5ztoFwiIAIiIAIiIAIiIAKOQOBEJ6+1JBRSZmGSIteTj\/yuJnH+8+2339pjjz1mmtwUZ9A6vAiIQP4J6AgiIAIiEEACgROdAWRkv\/zyiz3xxBO2efPmIBZPZRIBERABERABERCBwBNIuOhE0H300UdGYpnE8rvvvmvZJfKRPxeU85T15Zdftnnz5lmPHj3ytL92EgEREAEREAEREIGwE0i46Ny2bZs9+OCDLrFMenD\/91tvvdWyS+QjfzwbEbE5fPhw6927tzVo0CCep9KxRUAERCAkBFRNERCBMBJIuOgsXbq0E5wISJZJLA8bNsyyS+Qjf7waDkH75JNPGuGbrrzyynidRscVAREQAREQAREQgZQnkHDRWaJECWvevLlLLJP43rp1a8sukY\/88WglJgwxrL5kyRK78847La\/nmTVrls2cOTMtrVy50nbu3Jlp0rbYstm1a5fzyRXXvHMVw7yz47oTv\/zxSwWGe\/bsiUc3lZTHpD3DkOjrI\/t+4o0nZYPFuNAJF505rQ8ikLid+HnOnj3bdu\/endNd85SPczz77LPGDPnjjjsuT8dgpyFDhljHjh3T0qhRo2zNmjVKBcBg7dq1xgsEfvzxR\/HOI28xzN9vVfzyx497ZYIYxvSesWHDBroDpf0EmJBLu6Z6oq+P7Pv79eu3v\/b6H0jRuX37dhs4cKB16dLFWKaZpkyZYhdffLHz8yRWJ28tYvibbbFOCBXOz+s4L7300nwdfvDgwTZmzJi01LVrV6tatapSATCoXLmyVaxY0apUqSLeeeQthvn7rYpf\/vhxr0wFhuXKlctXP5JKO1eoUCEU92P6+si+n3khqdSOea1LIEXnW2+9ZaNHj7ZNmzbZr7\/+aojLESNGuDr27NnTeBPRp59+ai+88IJbF+s\/8+fPt8WLF7vZ88cff7zVqlXLJX\/2Op9MKvrqq6+yPXW1atXSvU2pRo0aVrx48eRNSVR2XCKKFSsm1vloMzHM329V\/PLHj3tlqjDMtrMISYYiRYqE4p5MX9+kSZO0\/r9x48amf2aBE507duyw9957z+rUqWMvvviiHXHEEUZg9qVLl1rTpk3dLHKskLwakyFw3xIay8Y88cQTo05iuummm9xp+MSCWb16dfddf0RABERABESgIAnoXCKQjAQCKTq\/\/\/57q1+\/vhsaBSphixCXp59+umG54skXEz2OuT\/\/\/DNZYpoY0ok2iem0005z5+HzwgsvtPLly7vv+iMCIiACIiACIiACIpA1gcCJzozFZdbf3LlzrWjRok6IZtyu7yKQnoC+iYAIiIAIiIAIBJFA4ETn4Ycf7pyMme3H24aY4YaPZc2aNQ3\/SiAy0WfhwoV2zDHHWMmSJVmlJAIiIAIiIAIiEBQCKocIRCEQONFZqlQp53j7ySef2M0332x9+\/Y1htFxyC1TpoxNnjzZmBWGj+dFF11UoKLzggsusOXLlxufpn8iIAIiIAIiIAIiIAI5JhA40UnJCYdEeKRp06bZF1984V4\/yYxxYnMS+4qA7YRN0luCoKWUZARUXBEQAREQAREIJYFAik5ebckrMOfMmWO80Wf8+PF21FFHuTALnTp1sg8++MD+8pe\/5PktQaFsaVVaBERABERABETgdwL6SASBQIpOQHie58IlEdy7cOHCrDLiexGsHf9Oz\/PcOv0RAREQAREQAREQAREIPoHAis6tW7e6oXVmr4ORAPH333+\/nXvuucbQ+ueff85qJREQgRgT0OFEQAREQAREIB4EAik6mSSETycB2Hfu3Gn79u2zp556yl5++WVbtWqVERSeSUbz5s2LBxMdUwREQAREQAREQAQSSSAlzx040YnAfOWVV2zjxo3uDUQMra9du9Y++ugjY1h96tSpXEP3aQAAEABJREFU9vrrr7vGIB\/53Rf9EQEREAEREAEREAERCCyBwIlOhtWJy0mIJGax8\/ahr7\/+2v773\/\/aqaeearzPlE9eg8n6TZs2BRauCiYCIhAHAjqkCIiACIhAUhIInOjctWuX\/fTTT8ZrLhGcUGUYHYvm2WefbZ7nuVSoUCHbsWOH+T6f5FMSAREQAREQAREQARGIP4G8nCFwotPzPDdL3a\/M9u3b7csvvzTCKNWuXdutRmji68nbiJjR7lbqjwiIgAiIgAiIgAiIQGAJBE508tYhfDfnzp3rJg0xhM4rL0888UTn0wnJGTNmGG8sqlOnjpGfdUoiIAIiEAwCKoUIiIAIiEA0AoETncWKFXMhkXj1JeGROnTo4IbRW7dubbyRqHPnztatWzdXF5Zl6XQo9EcEREAEREAEREAEAk2gQEVnTkm0aNHCnn32WcO6eeyxx9ojjzzihGjRokVt7969hoXzueeeszPOOCOnh1Q+ERABERABERABERCBBBIIpOj0PM+aN29ub7\/9tr3\/\/vtOcBI6qVSpUvbSSy\/ZlClTjNntCeSmU4uACIhAshJQuUVABEQgIQQCKTqzIsHwe1bbtU0EREAEREAEREAERCB4BBIuOtetW2fnnXeeSyyT+F6rVi3LLpGP\/DHDqgOJgAiIgAiIgAiIgAjEhUDCRafneVa1alWXPO9ADE6+V69e3bJL5PM8Ly5gdFAREAEREIHEENBZRUAEUpNAwkVnxYoVbezYsS6xTOL79OnTLbtEPvKnZtOoViIgAiIgAiIgAiKQOgQSLjpTB2VB1ETnEAEREAEREAEREIHkJJBw0fnbb7\/Z+vXrDd\/M3Cb2Y\/\/kRK9Si4AIiIAIJCUBFVoERCBPBBIuOjdu3Ght27Z1IZAIg5SbxH7sn6eaaycREAEREAEREAEREIECI5Bw0RlZU96lnt3kocjtmkgUSS8QyyqECIiACIiACIiACEQlkHDRWbx4cRcaieDvO3bsMN46dO2119pbb72liURRm0wrRUAEREAERCArAtomAsEkkHDRWaZMGfeWoS+\/\/NKeeOIJ461Djz76qJ1++ul2+eWX24QJE2zbtm3BpKdSiYAIiIAIiIAIiIAI5IhAwkWnX8rSpUvblVde6Sycc+bMsYceesh++eUX69evn5166ql29dVX26RJk2z79u3+LvoUgVwT0A4iIAIiIAIiIAKJIRAY0RlZ\/fLly1unTp3s\/fffd0Psffr0sTVr1livXr2cBfSGG26wTz75xHbv3h25m5ZFQAREQAREQASCT0AlDCmBQkGvd7Vq1ezmm292IpNg8TfeeKMtXrzYunbtai1btnThloJeB5VPBERABERABERABMJOIPCi02+gvXv3Ot\/On376yXbt2uVW79mzx\/bt2+eW9UcEUoKAKiECIiACIiACKUog0KKTwO+LFi2ye++910455RRr3bq1vfrqq8bko\/79+9v48eOtUqVKKdo0qpYIiIAIiIAIiEAiCOic8SEQONEZKTSbNGmSTmjeeeedbph96tSpxjB7xYoV40NFRxUBERABERABERABEYgpgUCIzsyEJrE7b7nlFjeZiIlD+Hbi4xlTAjqYCIhALggoqwiIgAiIgAjkjUDCRef69evt0ksvTbNoIjRvuukmw5o5Y8YMu+OOO4y3EHmel7ca5mMvQjYNHz7czjrrLBfAvm7dutahQwf79ttv83FU7SoCIiACyUVg1apVNnPmzIQkQujNnz8\/IeeOrDMMkqvVVNqUJpCklUu46GQiEG8igl\/Dhg2tWbNmtmXLFkPs4cs5YMAAyyw99thjtnXrVnaNeSIeaM+ePe3xxx93cUKffPJJF7Jp6dKlLmbo559\/HvNz6oAiIAIiEDQCiC3iJXfs2NESkbp06WJ9+\/a17t27J+T8fp1hAIugtY\/KIwLJRCDhojMS1rx582zcuHE2duzYHCWCxe\/cuTPyEDFb\/uKLL9yTNTe7YcOGWZs2bZzofP31141A9s8884wLXh+zE+pAIiACsSCgY8SYAEJr1qxZ9vOp3W1707tCmXbWvcxgAIsY49XhRCBUBBIuOo888kg3Cz1yGCOny8xeZ\/94tNiSJUusXLly1rRpU\/O8g0P7NWrUsEaNGtmKFStcCKd4nFvHFAEREIGgEdh7eAXbU+H4kKa6QWsOlUcEAk4gevESLjrx4WQWOqGPcpvYj\/2jVy1\/a3v06OGebE866aR0B\/r111\/jNqSf7kT6IgIiIAIiIAIiIAIpRKBQCtWlQKryzTffGP6cDRs2tCOOOKJAzqmTiIAIpA4B1UQEREAEwkpAojMXLb969Wq766677LDDDrPrrrvOihQpkou9lVUEREAEREAEREAEwksgQKIz2I1AaCfCN\/3www\/20EMPGZbOnJQY5\/NIH9WVK1cak5+UdhYIh127drkJX+Kdd95imHd2XHfJzm\/Pnj05udWFIg\/tmZckhgcvj7zwS8Z96Osj+36MVgcphHdJojMHbc+koa5duxqz6wcPHmwtW7bMwV4HsgwZMiRdmI9Ro0bZmjVrlAqAwdq1a42HhR9\/\/FG888hbDPP3W01KfhmulQ0bNhy4memvbd68OU\/3EjE8ePHklWGy9Zv09X64LT4JuXWQQniXJDqzaXtCJ7Vv396JlxEjRrhA9p53cDZ7NrsbInXMmDHmJ8Rr1apVTSn+DCpXrmxMNqtSpYp45\/GaE8P8XaepwK9cuXLZ3eZCs71ChQp5upeUK1cuNIyyq2heGSZbn0lf7\/f7fPbu3Ts7NKHYLtGZRTNj2SRAfLFixeyll15y4ZOyyB51U7Vq1axJkyZpqUaNGla8eHGlAmBQokQJo+3EO+\/XmxjmnR3XXarwi3pzC+FK\/Php17ykEOKKWuX8MMwL90TtQ18f2fc3btw4Ko+wrZTozKTF8b\/AHF6mTBl74YUXrH79+pnk1GoREAEREIGCI6AziYAIJCsBic4oLcerORGay5YtMywVL7300iGv4oznKzijFEmrREAEREAEREAERCCpCUh0Rmk+3gVPPE428RnttZzxfAUn581L0j4iIAIiIAIiIAIiEFQCEp1RWqZUqVKG0Fy+fLlllqZNm2a8QSnK7lolAiIgAiIQXgKquQiIQCYEJDozAaPVIiACIiACIiACIiACsSMg0Rk7ljpSdgS0XQREQAREQAREILQEJDpD2\/SquAiIgAiIQBgJqM4ikCgCEp2JIq\/zioAIiIAIiIAIiECICEh0hqixVdXsCGi7CIiACIiACIhAvAhIdMaLrI4rAiIgAiIgAiKQewLaI2UJSHSmbNOqYiIgAiIgAiIgAiIQHAISncFpC5VEBLIjoO0iIAIiIAIikLQEJDqTtulUcBEQAREQAREQgYInoDPmlYBEZ17JaT8REAEREAEREAEREIEcE5DozDEqZRQBEciOgLaLgAiIgAiIQGYEJDozI6P1IiACIiACIiACIpB8BAJbYonOwDaNCiYCIiACIiACIiACqUNAojN12lI1EQERyI6AtouACIiACCSMgERnwtDrxCIgAiIgAiIgAiIQHgK+6AxPjVVTERABERABERABERCBAicg0VngyHVCERABEciMgNaLgAiIQOoSkOhM3bZVzURABERABERABEQgMASSRnQGhpgKIgIiIAIiIAIiIAIikGsCEp25RqYdREAERCC0BFRxERABEcgzAYnOPKPTjiIgAiIgAiIgAiIgAjklINGZU1LZ5dN2ERABERABERABERCBTAlIdGaKRhtEQAREQASSjYDKKwIiEFwCEp3BbRuVTAREQAREQAREQARShoBEZ8o0ZXYV0XYREAEREAEREAERSBwBic7EsdeZRUAEREAEwkZA9RWBEBOQ6Axx46vqIiACIiACIiACIlBQBCQ6C4q0zpMdAW0XAREQAREQARFIYQISnSncuKqaCIiACIiACOSOgHKLQPwISHTGj62OLAIiIAIiIAIiIAIi8DsBic7fQehDBLIjoO0iIAIiIAIiIAJ5JyDRmQW7ffv22eTJk61FixZWq1YtO+mkk+yhhx6yrVu3ZrGXNomACIiACIiACMSJgA6bxAQkOrNovClTplifPn2sTJky9uSTT1rbtm3t9ddft169etm2bduy2FObREAEgkJg1apVNnPmzISkOXPm2Pz58xNy7sg6wyAo7aFyiIAIhJeARGcmbb9+\/XobMmSINWzY0EaPHm1t2rSxBx54wAYNGmSzZ8+2iRMnZrKnVgeFwNq1a13b8RmUMsW1HHE4OOy4\/vmMw+HjfkjEVr9+\/axjx44JSV26dLG+ffta9+7dE3J+v94wgEXcgcfhBIV+3mDFF79lfMbh8KE4JOzEMBRNHfhKSnRm0kRYJ5YuXWrt27e30qVLp+U699xzrV69evbRRx\/Zrl270tZrIXgEEEqjRo0yPoNXuuQoEeySmSFCa9asWfbzqd1te9O7CjztrHuZa+hEnZ86UwYYwMIVJsn+FPp5437ROXG\/6NyYZCUPTnHFsGDbQmfLnIBEZyZsFi9ebKVKlbLatWuny8FQe82aNQ1BumXLlnTb9EUERCCYBPYeXsH2VDg+AamuA5K481PnA2VwBdEfERABEUggAYnOTOCvXLnS+XIeeeSR6XJ4nmeFChWy7du324YNG0z\/REAEckpA+URABERABMJMQKIzi9ZHXBYpUuSQHMcee6zt2bPHfvvtt0O2ZVzBsFakQ7+WC25Cx+rVq11zqA3yzjzZGfrlL7JhsRXZsKTAU6GfDzyYFknQ+Yvsr7Nfhrz+DsRwyf6h\/QPtKIZ5\/w3l9zpM9r7T\/x25TinRfxJ4fonOOMGvXr26NW7c2E1G8p359VmwkzmYPEHzMiFM7PPGPtkZ+uUvvniilfrP4wWeDv\/yRS5B55OYiPNzTr8Mef0diOHjJob5\/+3kl2Gy38P5HaEJ0AbuphDSPxKdeWj47777zrCAFi5cONO9ubAGDx5sY8aMURIDXQPJcQ2ondROugZ0DcTtGkAToA0yFQ4h2CDRmUkj16hRwzZv3mw\/\/vhjuhwEjN+7d6+bZFShQoV02zJ+4eJq0qSJKYmBrgFdA7oGdA3oGgj3NYAmyKgTwvY9Z6IzbFT217du3bpustCyZcv2fzv4n7cRrVixwurUqWNly5Y9uEFLIiACIiACIiACIiACmRKQ6MwETYMGDZywnDBhQrq3D02fPt0WLVpkzZs3t2LFimWyt1aLgAiIQOwJ6IgiIAIikMwEJDozab2KFStat27dbMaMGda5c2d77bXXbODAgTZgwABr1KiRXXbZgaDPmeyu1SIgAiIgAiIgAiIgAhEEUkR0RtQohovt2rWzoUOHumH2\/v372\/jx453YZBZo5FuKYnhKHUoEREAEREAEREAEUpKARGcWzep5nl188cU2depUW758uX311Vf2yCOPWPny5bPYS5tEQAREIKQEVG0REAERyIKARGcWcLRJBERABERABERABEQgNoxqvHwAABAASURBVAQkOmPDMbujaLsIiIAIiIAIiIAIhJqARGeom1+VFwEREIEwEVBdRUAEEklAojOR9HVuERABERABERABEQgJAYnOkDR0dtXUdhEQAREQAREQARGIJwGJznjS1bFFQAREQAREIOcElFMEUpqARGdKN68qJwIiIAIiIAIiIALBICDRGYx2UCmyI5CH7fv27bOffvrJ1q9fb7t3787DEbQL3D755BP3Nq4333zTfvnlF0HJJQFdh7kEFiW7GEaBkstVYphLYMoeFwISnXHBqoMmkgA3148++siaNWtmp59+ujVu3NhOO+00Gz16tP3222+JLFpSnRuGZ555pnXt2tVGjRplffv2tX79+tmvv\/6aVPVIVGF1HeaM\/Pbt223r1q1RM4thVCy5Whlrhrk6uTKLQAYCEp0ZgOhr8hN45ZVX7JZbbrHjjz\/e3nnnHXv11Vetdu3azlr3wgsvGDfh5K9lfGvwwQcfOIYI9y+++MKWLl1qI0aMsAceeMAOO+yw+J48RY6u6zD7hlyxYoV769tf\/\/pX27NnzyE7iOEhSHK9QgxzjUw7xJGARGcc4erQBU\/ghx9+sOeee85at25tzzzzjJ1wwgnO0vniiy9agwYN7Mknn7T\/\/Oc\/cShY6hxy165dNnbsWDv33HMNMVC+fHkrXLiwnX\/++VaxYsXUqWgca6LrMGdwjz76aGvatKm98cYbNmPGjHQ7iWE6HHn6IoZ5wqad4khAojOOcHXogiewZMkS27hxo11zzTVWtGhRV4Bvv\/3Wunfvbl999ZW1a9fOiU+3QX\/SCGBx+t\/\/\/ue+48f5888\/2969e531ieHPadOm2Z133ukEAhZkhtq3bdvm8uuPGdbzSNcNXYc5uyqKFCliPXr0sAoVKtizzz5rXGv+nmLok8j7ZygZ5h2X9iwAAhKdBQBZp4gPAfzAmCQU2dlHnmn16tV28803u+E71r\/77rt2\/\/33GzdihotZp2T2+eef22WXXWbDhg1zIrNUqVLWvHlzmzp1qp144ol28skn23XXXWcTJkxw\/pxly5Y1JhXhqiB+5iaqdejQwf7v\/\/7Pic+MTHQdZiSS\/nvNmjWtc+fONmvWLHv\/\/ffTb\/z9mxj+DiKLD0YodD\/MApA2BYKARGcgmkGFyA0BbqyIyVNOOcUNnSOIWMcxEER8DhgwwFq2bOkE5siRI238+PF23HHHuSE8OjgEFflClDKtasOGDR0rRDn+m2Ts1q2b84HFJaFNmzb2xBNPOFEwe\/ZsJ0Zhz3KkZYr9wpJ40EGsDxw40IYMGWKbN2+2cePGOd9XGOg6hELWKZLhqlWrDDeOv\/\/974bAZE8xhEL2iZEc7oH169d398MLL7zQFi5c6HYUQ4dBfwJEQKIzQI2homRPgJA99957ryF4br31Vrv00kvt008\/tccff9yFRapbt67Vq1fPvv\/+ezdsx0Sic845x\/kkcnSGhJmwcOyxx\/I1ZRPDvQjtTZs2pasjHb0v0P0NuCHceOONVqJECSegEJKsQ5zja4cf7JVXXnmIP+cRRxxhxYsX9w8Tmk\/cD5hQ1alTJ\/vyyy9ty5YthmsCXF966SVnLdZ1mPXlkJEh7LDULV++3PkTc\/2KYdYM2frhhx\/a1Vdf7ULD\/eUvf7H+\/ftbuXLlrFixYmw2MXQYAvYn3MWR6Ax3+ydF7RFKK1euNPwOGRrHKvL666\/b7bff7ixwCKK3337bZs6caaVLl7Y\/\/elPVqhQISdMd+7cmVZHBNiYMWOsTp06boJR2oYUXFi0aJE99thjzsLrVw+xfc8999gll1xiWEf89XzCBD\/YjEOciAOGPGFHPgTpU089ZRwf\/1h88lgfpoRFE1eDv\/3tb\/bmm286twR8Xi+++OK0CTG6DrO+IjIyZPIfAqpRo0ZOdH799df6LWeN0Ll18LDdqlUrx4zfIw+P\/\/rXv9w9jt11HUJBKUgEJDqD1BohLAuCEjGJZSNa9RGYHTt2dDOnW7Ro4SYC8fReo0YNlx2LHEPtkRMRzjvvPOvSpYt99tlnzp8TqxQJUYDYQowdddRRbv9U\/VO1alWrXr26i6\/JUNvw4cNt3bp1zncTy5xvkfPr73meY8Ykoeeff95Z7tj21ltvORHPTHYScTvp1IYOHeomFZEnVikZjoM1jvil+Lk2b97cPM9zxWZoGAs816U\/IUbXoUNzyJ\/MGBIZgWgJWOkQoTzwiOEh+NJWrF271tasWWPcFxmlYANs582bZ\/yG8dHmNy+GkFEKCgGJzqC0RAjLgeWNjrp9+\/ZpvnA+BsQobxNiCJ0hcTqjtm3buqFLbrY7duzws1rNiIkIiCSE6N13321Y5LDE\/fOf\/3STYPDxxCKKn2Lazim6gAjCEsywJZOECIFEB9SkSRPnkhDpv+kjoNNnyJhJVsQ25UHgqquucnFOmVCEf9htt91mvKHoggsuSBNcFqJ\/CCFm9mNJz1htHmR4QPKtxboOMxI68D0rhvyWmc3+73\/\/241ciOEBZvyOiZ3Lb4\/RBtYyjF6mTBl7+OGHXWQJLJ74dTLy8+ijj7r7Hw\/f7Bv2+yG8lHJMIK4ZJTrjilcHz4oAgvCKK65wE3zoXPy8WCMRSgxf8uSOlY6ho0ceecS4oeJHN3fuXD+7+0QcYaXDgsdEBOJKcgxicuInRrikBx980CpVquTyp\/ofBCazqeng8bukY2LCEJwZgsMywgQYvwPzeeD\/WqVKlbRJMZ7nmT\/kiWBnXzo6P3+qfnIN8cCDVR3XA64dnyW+rMQ\/9F0OIhnwBizEgD8hJszXYWYMS5YsadWqVTN+lwiiSH4sY03nGsOXGMZhZkiEDkZpzjrrLOejztvBGNnhQZwYp0yY5OEcdw94nnTSSS62LuL0z3\/+s3uYxzIfZoZcU0rBISDRGZy2CEVJsJ7RGfmV5RWVWCKxcGD5ZD0WI6xqWNsI7k4HxXpfMLGNwO\/ceFlPwkqHSEAMYNXjPKwPa4IHkwrwYaWTHzlyZNp70xFR0fw3YYXIp4PCF5bhOYQW68OSuG7owLGKMzsd\/1d4cY2O3M+QB6WLLrrIuR8wSY38kWy4RuGHAAjrdQiTrBh6nudePIDgJB\/5Ixny++ahaMGCBTZp0qTITaFa5j7JKBCjEkwSQkjypjUmTjJ8DjdchngIJ+oEftbw5AGdeybuH9wH+C2HCpwqG2gCEp2Bbp7UKxxDt7169XKTgvzaMVROsHFullg2iBPJZCBet7hhw4Y0sUR+wh5dd911LrbkxIkTWZWWmjVrZvgvffzxx\/bjjz+mrQ\/DApZNLMJYPuCCgMdCiain42KiC8NzsPC8g\/6biHf8ZlmPwMRSDGMso8SeRGSxLZUS19igQYPSBSL36zdlyhTDQsTLBLAQESqK1wjiwzlq1Ch33fIGHb7zlis6fH9fRACMK1eu7OKbMsnIZ+vnCcNnThg2b97chfeBacY3ESGecK9hqJgRCh6EwsAtso5E6SAcF79JXIa4NyIkuS55EGciFvdS9uFBB+s6eflOgh9ileNwH2CdkggEgUChGBRChxCBHBNADNKp8KTOTRJrx+GHH+6G2HGA9y0bvu8hnToWj8gTZBxK97dhHWFIjhsyQ8T++lT+hCFik8kE+L1iXevTp48hHhFBnucZopOhOCa4wBseWECwDGNNadOmjSH6se7hS8dDAa4JDOl53oGJMuyTKgmRg+WXzjyyTjz8YF2nk+7Zs6cLIcV2HnwY5iQsEgHxmRH84IMPuggJ3bp1c\/FMcT3gjU0MqzP8yTIPAtOnT+cQoUk5ZYjLx1133eVC+9x0000uVBcRKl577TW74447jFiT\/I65pplYFBqAv1eU3ykz+C+\/\/HLnisBqBCQPQbgfcU36IxE8YDKL\/eyzzzauTx52+M1zLcIS8c7+SiIQBAISnUFohRCVwReTxH8kxmbr1q1t8eLF5gtJrEdYh3h6x38QIZnR9xDBhCWUoXRusogrHyFWUvbxv6fCJ4KGsEUZ64LgxN+LmeQIHcQ5w2yIRSyYzFpnHywkbCfcFLE7fV4tW7Y0LKNYkhBg+Cpi1cP3k\/1SNZ1xxhlGnFauNUQObOjk6bxx2eAa4vqDE0OaxCutXbu24e+KFRnGTEZjFj8TrEaPHm29e\/d28WKxoMKVWLE8+GT0PY4v04I\/OoxgyO\/6H\/\/4h5vol1OGXGcIS9qC3\/j555\/vrMx\/\/OMf3UMQYtPzUu+hJyethCgn8RsnP\/fEbvsfcLgP4OqBewwPiISJYzSCER78YPF7R2iyDw9QTCTyvHAyhIFS8AhIdAavTVK6RDytM2z+66+\/2saNGw0rHR00QhLLGzdXxACdGb503Fz92cCRYJg9zY2Wp346uchtqbSMhe3ll182JgUhjCLrxjAmPlwIHUQ4YolZ6wh5mPjWZPbBcsRQOxY+huXgzDAznRJDmMuWLTNinzK0Tv5kTlw7XBdERsACzMMNnLj2qBfXGvUmVBdD6XxH4OD7ysxpfGE5BsL9hhtucP6HWNCxACNMsSAxAYuQVFiW58+fb3PmzHFvu+JcnucZD0Tw5c1NnDMVE9cmdcfKhoUTZrlliC83bQM\/BBQiHesm13IqMousE4KSh0DcNWrVquVeN0s4N64trsnJkycboxYw5qGcUGf4aGOJ5zcOb65LrjNEP4KUeyX5+S3zYBR5Pi2LQBAIhEJ0BgG0ynCAAFY1bpJYNZh1idg5sMWs2e8+mViQsNJ53kHfQwQUw5t+XqyZWEfwCeMJ31+fap\/UjaFurJQMmePnhdChnohF\/AcRk3ynE6ITe+ihh4w4nXRC+BiyjU4caxzCC0scMTfvv\/9+Q\/yzPVUSwpzJFkywwIrGJApENsPdDOfSoePeQWeN0MS1g1BccMZixH5YhhHnDF9iNSL8DPyOOeYY91YmIiLQ+cMbbggtfOgQpHxH8NIGTIjDd5F1yZSoC5OoCI\/FdULipQK4Yvj1wBWG3zCC8w9\/+IPh\/8qDT14Yep5nWNmJLAFn\/xyp\/Ml1SBgjHnD4\/cKRBxTuc1je2c7DOfdKXDl4wMHCyYM4XLh2mbCGwMT\/mvx8R6yGhSEclJKPgERn8rVZ4EqMkMRKgeM61iO\/M45WUG6KWIqwZhDiaMSIEe7NGuRFSNJxsczNF0sA+bE2IRwYLoo8Nvm50ZI\/lROzpbEGM3yJjyCdPBwQRzCk7nQ6iCqse1hHsI7CDo50XOTBQkIehjPxo8NCSsfGtlRIiCJCyvBmG6zmDIMz05dXphJqCyHO9UIYJK5VhCP1hhHXGsskhCvCnmFeJm6wD+u\/++47Y9gd5gharknWIxgQDliTSVihWI8PLcKT5WRJWMCJNcrEFYQk5Ua4U18YInJYB4NTTz3VxdLE9YCXM7DeT5kwtMwY+vuF5RN3FkYqcI0hNBwPPjzoIES5ZrF48hBDexAJgWuR3zN8uFZxT+JhlAcDHta5F7JNSQSCTkCiM+gtFODyMbyG5Y1hHDoqOnqGM\/Fjg68WAAAQAElEQVQ98oVOZsXnBoqYxIKHTxwiirz4eTGUFGmloyPH2oK4TTXLHHXOKmF1evrpp128PcQPwpHkeZ4hsBBMdDhYR2DmDwkzcQiWiAQEFHw9zzP8E8lH54Z1L6tzJ9M26od\/Lx00LwW4\/fbbza8fn0zAwH8VVp7nOYslfLjW8I2bNm1aWnU5Fh0+iWU28DCF\/yziAL9OBIPvikDHTyJgPCITsYnFmaFj9k2WxLUCD3ys4cUkPh5yGLLFyo4vLA881AeO\/N6ZVMVoBUO\/rPcT3OBHYpn1WTFke1gSrgjvvfeei3+LYKTeMGKCG\/6xWN+xqrOO3zzWY4Qo91QSls2PP\/7YeCDiOkd8ko\/jKIlA0AlIdAahhZKwDFh5COyORY1hNnyIiFuIryU3T2ZQZlctxCT+TFhROB75Pc+za6+91okCBBU3WYaL8P1ENDCsRL6wJDoTZpZjrWOWORY6LHiR9WfyD8OdiB34e55niFUsVFjmcEOgM4vcJ9WWuU6wtlN\/rEI5qZ\/nHbjWsNIh1jkG+yGksOQhRrkW8fNk2J3A24h2RACJvCQEGJ3\/1KlTDWsVZYA725IlYZlkqHfv3r2GywpWTb8OnucZ\/tM8\/OAz7NeJWdGI1G+++cY+++wzf7X7zC1Dt1NI\/iAmEeO4FDBZCHcMfNe7detmXHc80CAoifSBrycjHYjO0047zUgsc2\/0H3pCgk3VTBECEp0p0pAFWQ06Z4ZpsXRiAWG4nJsh8ePovBne5AaaXZkQk1jo8DPEf45hI\/bBQuQPs2\/evJlVFtnJuxUp+IfOiIkUWIAZAmZyADPX6cAJmO3P5oc3w+k+AoYsEQMMm\/vraAMsSwwhM3xH6BV\/Wyp+IrKxhMPFv478eiKomKyCXyIiCR9NhuLZ7l9rhKfBOsk6RCSzgHmtIAIfUcUQ+uDBg41rljyplrB+M+qA8MQCnFn9YMs1yax+z\/OctQ1\/Q0QQ9wV\/v2Rl6Jc\/p59cc9HyZvZbZj1iHmsm1xYWctwxEKE8GHL\/RIwSsouRDCYEYtnEDYnv3Cd5uEfwRzuv1olA0AlIdAa9hQJYPm6WdNJ0UFgqI4vIDZVOyPNyFqYDscrwEFYlZgdjpWLWJnHmEAqIgsjjp+oynTkWSXzpEElYNbHqIgYYjqPeTCLAIsJEFtazjoQoxfeLiQYrV650MQ8ZPu\/UqZOzjNCpIQLIm6oJ0Y3bAWKIYXCEEZN5uD7xt2QiERZ1RD0CCb5Yh+GBuCT0ERY+tuMvR3vQwWOBhzfXKFZn8qdaQvAwOoFrAA860eqHWMKaySQ0JqQxEQuRyW8dqxzhuhBNkfviQpOqDHnI8R\/m\/AcYv+5cO1n9lvld4y5EhAO4cxweeJhI5HmeYW1GzOIWgs811x2jSVzXYfgt+xz1mZoEJDpTs11jViuskMOHDzdujHQ8HBjLGsNCTCDge8aEBZQOHcsJN8+M2yO\/Y8HEqokPIsKTGy4Cghut5+VMuEYeL1mXGVJDDDG5B9cErHHMPCdUkm9d8zzPCPXDEBx5\/dn8WEZZh5WEzgnxiaWONzclK4\/clpvrCB9D3hREB44wGjlypHszFT7HWIsY\/mYInmUsSEzM4vrkesOFgYDbWEIRqExay20ZkjU\/bhhEksAnFY4Z64GoIiIArgb4biLiEZlwJi+iHe642mBdZ12qJsQm90HuT4zsMBLz\/fffp6tuTn7LuGswmuF5nuEny\/H8gzDcjn8tPrSMcPjr9SkCqUBAojMVWjEOdeDGimXjkksuMWaaYwGio+ZUDONifWOCCqKUjhvrJzOnuVEyPIdPGDdW3pLhz4Jl32iJWdgMGTFhAQteZtaWaPumyjomssANS5LnHRDbiIC1a9e6cDR0\/NQV6xEiEwscE7AIm4TAZzYs7cFwOj5ffrxI9glLwlKHEGKSBhY2On+GxrFcYvXFak5YHqzJ+F0SOgkRBR8soswkZgiTYyBCWR+GxAMkQ7pYPBGgkXXmPsAELSzI119\/veG2wSQjWGIZRmRGinbycpzIY6TKMtcL\/pVcT9QJn3S48DIGvvspJ79lfsdcY\/DmmsMyygMRn7gm8aDEQ5R\/zLRPLYhAkhOQ6EzyBoxH8ZkBTEgd37KBtY0OhuEezocopJNmOI1hScL5XH311YZ1CGGKeGJ4HAsIgvS+++4zOif2zSwx5MSNOPKJP7O8ybAe30Jeg8hnTsrLEDkuCzBETMIetrwFh86eqAC4HXAsxDzfsXbiq4lVlOF1hpdpG1iSL4yJ6wc3BIYhmWiBu0dmHLDqed4Bge95nnv1IkOY\/jBnZvul2no41KhRw3jAwcIWWT\/P8wzxg9hkIh\/3AB4S8TnE0o7LAvm57rAS89tHqLIu1RITzhju5ndH3bnWqCOTHBGQPIDzPae\/ZR582JfjIjaJAILvJp\/4bXMcjqckAqlEQKIzlVozn3VhWJxOl6d5fI1wcqezyejPxnAvQ7m8zQVhRGJCEBYALG\/4YjJMzGzXgQMH2po1aw6Z3ZrPombcPRDf6WyZZc6sfiwVdMQIGDokxHhWhcR6iRUOq\/KDDz7oeDHsxlAvw8F0THRIHAf\/TIbPYY6gZ0Y7lqesjq9tBwlwjWKhYjIGVrqDW8K7RBB7gtzz0JPxQQlGPBD6dPBZxJrHfaBZs2ZuNeKea5KQUiy7lSn2Bwa4AvEADSd+7wh27m\/8Bj\/++GNX49z8lrlH4EbDG5kY6SFMFf7YsHUH0x8RSDECEp0p1qB5rQ7D2vgD4pfFZB6GJvGLw7IR7ZjcFJlAwIQNEuIIoZWxw8HqRF7EUrTjpMo6rByI8A4dOrhXIBJmBwGJ3xc+mghH\/MEyqy\/WI1wMGB7HosmwMG2CsGRoGKsxPoech2PQAcIcv83M2oh8YUxwxjLPA4\/PCw6spy0YvkTgE4EB0cC2sCeuUwQQwodrEEEVjQnriVgxadIkw\/WGAPHR8qXqOqyThI\/CwouLi+cdCLvF7xELKJOrcvtb9rwDb2TiGBnvn8nHUSUWgawJSHRmzSc0WxleY7gIkUhnzY0zY+XpwHmqz7g+q+\/4gmI5wQ80q3zJvI2OmHA6+FQixLGiYa3AdxD\/QN8qDN+s6smwOB0\/k4JY9vPSkfFaQtwYNLHAp5L5Jy8QwOKM8MRHllBJCHT8i7E682Yihi+jXeOZHzW1t2A9Z2QD305mpvMqWq7ryFrz+2fyFdu5RplcFCaRhF81v+3TTz\/dGBXCbxirL6MMWECZXIXLEcz4\/eq3DAklEUhPQKIzPY\/QfkNsEo6HmZMIpcgOh2V\/UhGWN0SQDwrrEYIooyWTfRg2YtY7nT0zXv19Uu2T2fxYKZk8hfWMDtyvI50yPpg5tUYSLonhOQK94zdGJ4elk2U6Njoz\/9j6jE4A\/rxGEIsmAbaxSuHygVUT6zNuIFzv0fcO71rcDXibE76EPDRxPSOiGDqGJ0PwPFzhO4ygD4toR1jiLkSYI7jw28RXnWgbTP7himEkIuMMfv2WIaMkAukJSHSm5xGKb1geSZGV9TzPGBpm8gWhZhg6YjsilM6HN7EwKxVBFCmg6MDbt2\/vQvngK8fkGYbo6PCZjETHREdPR8bxUiFRR7j4dSFkCkPf+MIiMv31kZ90UhmZR273l+ngEPB07kzWgiOWOSxPuC\/4+fSZNQFcOu6\/\/37jmuR6ZIY\/8TaxPiNKs947vFu55hDpDJtjuSMWLxOJnnvuOcPnE+sew+sMBYeFEvdD3I4Qnvhe4nfJUDoPgIxi8HvF7xXLL\/cBQpbx0K3fciCvEBUqwQQkOhPcAAV1eobGsE4w5MMTOYnliRMnulcmUg5CzuCLyI2TzgUfRWYBY8ljiBI\/T75HDhPzHdFJB8VMao6JeKWjJ6YkQ84IT46f7AmGdDbMHGeY1rf4MuRN2BlC9dDRUE9EJpMtGLJkMhFD4zBnX986Qr6MiTbAp44QSJyLvIShIQh3xrz6nj0BHgIQSJEPStnvFe4cDBfzWltmsnNNE98U0c4nv\/fI33+ykmKExv+tZlUHQmrhv4q\/O7PzqbvneYZv51\/+8heDC\/dQjoGvJ\/7wjBTh3qHfMlSURCA9AYnO9DxS7hs3V4bH8G3DV5MKYrlkyJvAxrwzul+\/fuZb4Yg9x1A4+zBDk5soYhMxGs1CxDpuvgylYxHBIsqsdYbcmQmfCp09VgvfvYDwRPin4teGdQOeDEtiIcIaCVeENxMzmOTDsDtiFeFNJ\/Tf\/\/7XmKAVaSnlGJGJjo33fNOxIZgit8VtWQcWgQwEsBbjhoAITYXfsV89hr1xWWEEht+2vz7aJw+PBM7noRL3jMg83FO5P2L5xd+TBxzuk57nGQ+K\/O71W44kpmURMJPoTOGrAGFDCBPEIDdMJgEQOJzgzgwBYbnAP4sZvYgphs\/paBBFfBLo\/fHHHzcEU1aYPM8zOqY+ffoYrx7EvylVhtNhiEhHTBMQn5ilWDJwQ\/CZILzxecPCS0ezceNGY6IF\/q9Yioi\/iW8mw7uIVcS+bx3xj6FPERCBgiFQrlw542GOh0QslVmdlSF0JlfhU+0\/mPv5EZmM7hBCCZGJgOUBlBEeJlBy7\/Dz6lMEckog1fNJdKZoCzOUzlCYPzSOxQ0rJk\/efpURk0z04XWWCCk\/ztxpp51m7Iu1krA9fv4wfTI7FUstHHyGvFkpY8xSnwlWT+IU4kOIPyxDcsyYxlIUyRwrcpUqVYxA3P6++hQBESg4AohFLJKckUlBPCiyHC3hOtOwYUPjjWukjHkqV65sPKAjYBcuXGie5xm+7NOnTzf5YGekpe8iIEtnyl4DZ555pvGEjrAknh7WuGiVZTsTLtj24osvumF2rKIMDROQnNAyDB2xPUwJf0z8URHkxCGlk8qMYW64\/PDDD4YvKMIzN\/sprwiIQOwI+BZJfuf+w3ZmR2cEA+sooxkMzUfmW7lypXFfwMrJgzvbEKGsY1lJBEQgPQFZOtPzSJlvWNj8EEj4ZnJTzKxy+CDim8STvD\/cxPAxwpNXYiK+sto\/s+Mm83r8KZkYQKeCUMxYF6wjTPTBTxa3hMjtiEq2R65jmQ6LWekM2eEPxjolERCBgifgeZ61b9\/ejj76aGPSJL\/NzErBvZRX0S5ZssSI3rFo0SLjN4\/LDL9njoMPvP\/wntlxtF4EUoJAPish0ZlPgEHd3fM8F8aobt26RggUhnwzKyvDTQy9cyP1RSd5ma3J\/ohO\/BJZF5aEpYKO5rDDDjP8MBGS1D1yYhZ+sHRA\/ja2w4+ZrvjF+hZifMGId0i8Trbj94qlhfxKIiACiSHAKA8P1ohJhsczK4Xneca9AFHJQzguNzyUd+3a1YhIwTbEa6QbTWbH0noRCDsBic4UvgJwlufJ3A+BlFWIkBo1ajjfJCx7PhL2JyA5ceiYEOOvT5VPrLdMAhg9erQRFJswJ5GMmGTFRCzCQeGjycQAvvsTs4jdhxUZTj4TZrYTRJvZ\/Oecc47VqlXLdUzEO8Q\/jJh\/4yLYHAAAD9BJREFUiE8\/vz5ThoAqElACjDoMHz7cmBjJ77FNmzbmPxAyEQg\/a37L3AsyqwKCEpGJryYjHFhHmYjJZ6pMmsys7lovArEkINEZS5oBPBYWTML5RL49IzfFZBiYcEF33HFHbnYLfF4EJLPSiZv54IMPuiE2hsmoL7ExqYDneXbttdcawdkZRkcs+pOKmJh11llnuYkD5PUT\/lwMuY0YMcJ4gwnxNfGpRdgyKQtXBj+vPkVABOJLgBGeK664wpgwSbxcfu+M2gwYMMD5r\/N75cGaB3PEJKM9lsU\/JgzieoO1k4gdnudlkVubREAEMhKIv+jMeEZ9L1AC3FQJ3cMwsf\/2jGgF4ObMMDA35sjtDL0zISlyXTIvY+GAhy8gGRrDYjF58mTDB5Y4mjfeeKPxnXrSsWDhYFidobSsYpaSn4RVhCF2rJpYRoYOHWoIf9azXUkERCD+BPzfOpOA+I0jKkk9e\/Y0wpfxEEgpeBC89NJLDb9M\/NpZpyQCIhAfAhKd8eEaqKNGhkCKFh+SGy2BknGYZ2g4UIWPUWEYYqPDwZKJ1Ze4m8Qsvffee12MUeJq4qOJQGSInfBH\/nAbQ+oITvy5CBQdoyLpMCJQoATCdjIeEImJyygF\/pvUn1i7iFFca7gfEH+TB2seNPHjJloHD9\/kVRIBEYg9AYnO2DMN3BEJgUTIH4aJ\/bdnUEhursSp4+0chQoVMkKC8NYctqVamjNnjpsQhF\/lG2+8YXRE0d6ygijlFZdr1qwxAsEz3MaQGvx27tzphuERsKnGR\/URgWQnwEMiItKvBxP8GOFBTLIOsckoB\/eC6667zhjdwScbAcoDNyMdkRZQ9lESARGILQGJTost0KAejWFirHtY63jC5y1EzZo1M95OxFDwhAkTUjqYMUNohIVCTOLPmVk7eZ5nV199tZv8QzD4+fPnu6z4ceEfm5O4fm4H\/REBESgwArzIghc3MLGHB0VOjP8mD45YOfFLx4d73bp17hWVt912m\/uN+5E5PO9AtA9GPIjLywM5x1ASARGILQGJztjyDPTR\/BBIvAITax5WP5bx9eTGHOjC57NwDKHxphD8U5nkk1VcPiyb+HHS8SA8OTX7Y+3kHcyIdqworFcSARGIEYEoh+EhGQtmlE3pVvHWIEYpIv0y8U9nxjphzQhhRvg3LJtYNbn38eKHjRs32uOPP25E6CAKBW9y4\/6AL3y6E+iLCIhATAhIdMYEY3IchJsq4X4qVapkvOIxsxnYyVGb3JfyxBNPdAGhFy9ebLw9hGG1zI5y6qmnulnr5N21a5fLRmxNhDvWT+JuupX6IwIiEBcCiE2iR+BvHunSwqS+jA+NPBT6fplDhgxxM9P9Qr3yyivGMPvDDz9sPFCyHt9sjk+4JH7PhERjPfdITfiDhJIIxIeARGd8uMbyqDE9FkPEvEkDS57v6xTTEwT4YJ534C0kxOpjZjk+XZkVF4smvrA\/\/\/yzMbGIfJ7nuTeS9OvXzxiqY52SCIhAfAgQO5jfGTPP\/WgSDJ3fc889RhiyjG4yWDDxy5w1a5b5M9MpGb9fhCqJ7whYwprxZrDnnnvOGJpv164dm5REQATiTECiM86Ag3Z4nuKxCgStXAVVHtwImExAXD4mVdGJRTs36\/fu3WuHH364RfJiohVx\/bAWR9tP60RABGJDwPMO+lky5I11kwdBfDe3bNliGX+\/nncwPy4w\/qQihtnZt1u3boZbUfPmzY0XQfiTCbGC5qzEyiUCIpBfAhKd+SWo\/ZOOANZehtWYxT5jxoxDys+wO64HdFoEeMcickgmrRABEYg7AYa7easar6r0XWKYFJhZXE3yE42DUQz8N\/ktY\/3EOoofNuHSWrVqZfz2sYzGvQI6gQiIQDoCEp3pcOhLXggk2z5MEsBaSbkJGcWEIZZJDMXxyrwnn3zSzj77bONtJqxXEgERiC0BRhNI0Y7Kb5LfIQ+IWDlLlixpvksMIw+Z+W9yrHPOOceqVKli48aNc2GR\/PwLFiywr776yu6\/\/36LFi6NfZVEQATiS0CiM758dfSAEjjzzDOdoPT9v7CIMORGIHiCxDdq1Mi9Ok\/vVQ5oA6pYSU0Af0xmmyMoM1aEmeS4wPCiBh78eD0tr5NliNwfUsdKiQXT\/\/1GHoOJf7gR+XF1eZCM3J7Cy6qaCASegERn4JtIBYwHAXzDCBBdoUIFFzSeIbkOHToYAaSJ7UfnxlBdPM6tY4pA2AngG401EvHJW4IiebzwwgvGjPLRo0e7KBtdunSxf\/3rX8ZkH4bFcYnxvIP+m4R8W7VqlTsEApPf7nHHHWfMVuc3zW\/dbdQfERCBhBOQ6Ex4E6gAcSeQyQmwlnTr1s3w3WQGK8HzP\/zwQ2OiApaSTHbTahEQgXwSwMXlxRdfdKMJDJ37h0OAfvnll0Z4M36f\/noibTAJiMgTvksMD4WEgONBsU2bNm6SENZT\/DZ79erlfsdnnXWWeZ7nH0afIiACCSYg0ZngBtDpE0fA8zxjOJ2g8YhNfzZr4kqkM4tAeAgwa3zlypXGwx6vn6TmuLng50nkCL5HJkQmvti8xtIPiYTIHDBggDGk\/tZbbxkvfyAuJ8HiI\/fVcjAIqBQiINGpayDUBOjI7r77bjv66KNDzUGVF4FEEGA4nJnmf\/\/7311AdyygvIqStwgRvD1jmRo0aGD8ZsmPhdPzDgyzM0Fo2bJl9vrrrxtD6xn303cREIFgEJDoDEY7qBShJqDKi0A4CTCEnnFC0Pnnn+8sl\/h0IkojyeCfWahQIeP1mGPHjjUso5HbtSwCIhBsAhKdwW4flU4EREAEkp7A1q1bjYDtTZs2NYQmvpjMRvc8z5gohHXTt14Sh\/OPf\/yjTZw40SZMmJBOWM6cOdOKFy9urVu3tnLlypn\/lqGkB6QKBIOAShF3AhKdcUesE4iACIhAeAkwfE6826efftp4O1CzZs2MGejffPONg8JwOQHgGU7Heok1k2DutWvXNvw1u3XrZgSGv+OOO6x\/\/\/6GZXTYsGF2\/fXXG3ndQfRHBEQgKQhIdCZFM6mQIpBQAjq5COSJADE3EZAIy6lTpxpiEYvnZ599Zuedd17aMQkCjxXUD+herVo1Q4BiBSUW51133WWffvqpMdmve\/fuaftpQQREILkISHQmV3uptCIgAiKQNASYnf7dd98ZlkyEJwX3PM9+\/fVXmzdvnvFCBvw2mUBEQPhffvnFCVPW8dagBx54wBYvXux8OBlaJ1C8wplBUSmcBJK\/1hKdyd+GqoEIiIAIBJIA\/pfFihVzVss333zTHnroIcOiefLJJ9uVV15piEiGybGInnbaac5XkzibCMxAVkiFEgERyBcBic584dPOIiACQSCgMiSOADPIiXPbqlUrww+zbt26xpA6E4WOPfZYZ+X86KOPXPD2kSNH2pYtW+ySSy5xbwJDeDJs\/vHHHzv\/zB49eriQSJMmTTLidSauVjqzCIhAPAhIdMaDqo4pAiIgAiEgwDD4kCFDrGfPnk4sPvbYY85ayWsrO3bs6N72deONN9q0adOMEEj4Zy5YsMCGDh1qvEXolltucfsRcxNcNWvWNIK8Dxo0yIlQ1imJgAgkDYFsCyrRmS0iZRABERCB8BH44YcfDJ\/KFi1auGFw\/C+xakaSYBiciUF33nmnYcXkDV9\/+9vfnKhkf2af82rL6tWru2F1\/Do97+BrKefOnWuEU2rYsGHaYXmrkOcdzJO2QQsiIAJJT6BQ0tdAFRABERCBoBNIovIhAhGbBGn\/5z\/\/aYQyQnAyk3zKlCnpasIw+B\/+8AcjJJLnHRCK3377rYvJSQzNWrVqufyIVWJu4sPJ24OwbGLN\/POf\/+x8O\/HndBn1RwREIKUJSHSmdPOqciIgAiKQMwLMHB8+fLizSCIQe\/Xq5WaY83rJl156yUqXLm0IRSyY\/hEZXidIe4kSJQwhSdB3Arez\/d1337X777\/f+XDu3LnTveby66+\/tssvv9zOOeccF\/id7QMHDrSiRYuyi5IIiECKE0i06ExxvKqeCIiACASbAFZI4mYyuefRRx91ohOLZu\/evY2wRZ7n2bnnnut8NdetW2fff\/99WoUQoosWLbIHH3zQWrZsaQR8f+aZZ1ww9+OOO86tv+yyywxrKZbS6dOnOyvoG2+84QLEY\/lUCKQ0nFoQgZQnINGZ8k2sCoqACIhA5gTee+8969q1q+F7iV8mopHg7JF7eJ5nxNIk\/BFC09\/WrFkzw4rJMbByvvPOO0agd19IEo+TsEnsyz7ly5e3Cy64wBo0aCDrJkCURCBkBCQ6Q9bgqq4IiIAIRBLgzUDNmze3zZs3u9Wed8A30335\/c+8efPstddec6+xJCzS76utUaNG7s1CDLMfffTRxjC7v23VqlUu+PtJJ51klStX9lfrUwREIMQEJDqzaXxtFgEREIFUJoBQJED7YYcdZlg5CdTu15fJQAR179q1q+G7iU+nb7UkD\/vefvvtVrZsWevTp4+b7f6f\/\/zHhUdq166dFSpUyO677z7DQkp+JREQgXATkOgMd\/ur9iIgAiJgzB7H9\/Lzzz+3iRMnGn6e+GoSAqlv376GECWo+3XXXed8Mjdt2pRG7YQTTnBvHKpXr54x2x3fzYcfftjOOussGzNmjGUcqk\/bMXcLyi0CIpACBCQ6U6ARVQUREAERyA+BIkWKGNbOqlWr2rPPPmudOnVybw1ipvq9995rBIDHF3PNmjXGZCOG1Qmr5ItPgrozOYhheGJ3Lly40AYPHuwCv+enXNpXBEQgtQhIdCZ7e6r8IiACIhADAgjHzp07u7cIEUuTgO+8nhIxeumll7pwSbyyklntTA7Cqsmsdobcie1JEZjtXqlSJQ2nA0NJBETgEAISnYcg0QoREAERCCcB\/DCZWY4vZtOmTd2M9UgSiEpEJ+LzhhtucD6bS5YssbVr10ZmC+WyKi0CIpA9AYnO7BkphwiIgAiEggAhjXr06OHCIA0bNswIGB+t4ohPXnGJD+jIkSONmJzR8mmdCIiACEQSkOiMpKHlOBDQIUVABJKJALE3CaM0bdo0+\/jjj7MsOjPePe\/QEEtZ7qSNIiACoSUg0RnaplfFRUAEROBQAoRBItA7fpvPP\/+8+ZOFDs2pNUlFQIUVgQAQkOgMQCOoCCIgAiIQJAL169e3q666yhYsWOBeaRmksqksIiACyUtAojN5204ljw0BHUUERCADAc\/zrHv37nb22Wdb9erVXdzODFn0VQREQARyTUCiM9fItIMIiIAIpD4BxObo0aPt4osvNs+T32bqt3iia6jzh4GARGcYWll1FAEREAEREAEREIEEE5DoTHAD6PQikB0BbRcBERABERCBVCAg0ZkKrag6iIAIiIAIiIAIxJOAjh0DAhKdMYCoQ4iACIiACIiACIiACGRN4P8BAAD\/\/1RfhXYAAAAGSURBVAMAXqkGEu76o94AAAAASUVORK5CYII=","height":268,"width":446}}
%---
%[output:97ab26cf]
%   data: {"dataType":"text","outputData":{"text":"          <strong>Variable<\/strong>          <strong>MissingCount<\/strong>    <strong>MissingPercent<\/strong>\n    <strong>____________________<\/strong>    <strong>____________<\/strong>    <strong>______________<\/strong>\n\n    {'DFS time'        }          0                  0    \n    {'DFS status'      }          0                  0    \n    {'Stage'           }          1            0.25316    \n    {'Grade'           }          9             2.2785    \n    {'Residual disease'}         43             10.886    \n\n","truncated":false}}
%---
