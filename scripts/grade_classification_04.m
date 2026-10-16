%[text] Scientific Programming - TCGA Ovarian Cancer dataset
%[text] 04 - Grade classification
clear;
clc;
close all;

load("results/cleaned_data.mat");

if ~exist("figures", "dir")
    mkdir("figures");
end
%%
%[text] 1. Prepare G2/G3 dataset \
gradeMask = ...
    strcmp(clinical_dfs.GRADE, 'G2') | ...
    strcmp(clinical_dfs.GRADE, 'G3');

X_grade = expression_dfs(:, gradeMask)';
gradeLabels = categorical(clinical_dfs.GRADE(gradeMask));

fprintf("G2 samples: %d\n", sum(gradeLabels == 'G2')); %[output:2afa8693]
fprintf("G3 samples: %d\n", sum(gradeLabels == 'G3')); %[output:0a5dfc72]
fprintf("Expression matrix for grade analysis: %d patients x %d genes\n", ... %[output:group:5c1f6fea] %[output:118c24c2]
    size(X_grade,1), size(X_grade,2)); %[output:group:5c1f6fea] %[output:118c24c2]
%%
%[text] 2\. Cross validation scaffold
%%\ 5-fold cross-validation

rng(1);

cv = cvpartition(gradeLabels, 'KFold', 5);
%[text] Five-fold cross-validation was used as a compromise between having enough training samples in each fold and still obtaining repeated held-out evaluations.
%%
%[text] 3\. Compare PCA variance thresholds using linear SVM
varianceTargets = [70 80 90];

predLinear = strings(numel(gradeLabels), numel(varianceTargets));

for i = 1:cv.NumTestSets

    trainIdx = training(cv, i);
    testIdx  = test(cv, i);

    Xtrain = X_grade(trainIdx, :);
    Xtest  = X_grade(testIdx, :);

    Ytrain = gradeLabels(trainIdx);

    %%\ PCA fitted on training data only

    [coeff, scoreTrain, ~, ~, explained, mu] = pca(Xtrain);

    cumulativeExplained = cumsum(explained);

    %%\ Class weights from training data

    nG2 = sum(Ytrain == 'G2');
    nG3 = sum(Ytrain == 'G3');

    costMatrix = [0, nG3/nG2;
        1, 0];

    %%\ Test each PCA variance threshold

    for j = 1:numel(varianceTargets)

        target = varianceTargets(j);

        nPC = find(cumulativeExplained >= target, 1);

        % Training PCA scores
        XtrainPCA = scoreTrain(:,1:nPC);

        % Project test patients into training PCA space
        XtestCentered = Xtest - mu;
        XtestPCA = XtestCentered * coeff(:,1:nPC);

        % Linear SVM
        mdl = fitcsvm( ...
            XtrainPCA, ...
            Ytrain, ...
            'KernelFunction','linear', ...
            'Cost',costMatrix, ...
            'ClassNames',categorical({'G2','G3'}), ...
            'Standardize',true);

        % Predict held-out samples
        predLinear(testIdx,j) = string(predict(mdl,XtestPCA));

    end
end
%[text] 
%[text] Evaluate PCA variance thresholds:
for j = 1:numel(varianceTargets) %[output:group:158e9c21]

    predictions = categorical( ...
        predLinear(:,j), ...
        {'G2','G3'});

    cm = confusionmat( ...
        gradeLabels, ...
        predictions, ...
        'Order',categorical({'G2','G3'}));

    sensitivityG2 = cm(1,1) / sum(cm(1,:));
    sensitivityG3 = cm(2,2) / sum(cm(2,:));

    balancedAccuracy = mean([sensitivityG2, sensitivityG3]);

    fprintf("\n%d%% PCA variance\n", varianceTargets(j)); %[output:48899e97] %[output:9257498d] %[output:4fe855a2]
    fprintf("G2 sensitivity: %.1f%%\n", 100*sensitivityG2); %[output:5d26a749] %[output:60d2921e] %[output:00393f65]
    fprintf("G3 sensitivity: %.1f%%\n", 100*sensitivityG3); %[output:175bf82a] %[output:12623b04] %[output:2d490319]
    fprintf("Balanced accuracy: %.1f%%\n", ... %[output:7b5ad979] %[output:371d078b] %[output:824c2b89]
        100*balancedAccuracy); %[output:7b5ad979] %[output:371d078b] %[output:824c2b89]

    disp(cm) %[output:09223198] %[output:5407ebf5] %[output:16360e11]
end %[output:group:158e9c21]
%[text] PCA thresholds of 70%, 80%, and 90% explained variance were compared within 5-fold cross-validation. The 70% threshold gave the highest balanced accuracy (57.8%) while using substantially fewer PCs, so it was selected for subsequent SVM analysis.
%[text] The comparison of PCA variance thresholds was exploratory. Because the same cross-validation procedure was used both to select the threshold and estimate performance, the reported balanced accuracy should not be interpreted as an independent estimate of predictive performance.
%%
%[text] 4\. Gaussian SVM using 70% PCA variance
varianceTarget = 70;

predRBF = strings(size(gradeLabels));

for i = 1:cv.NumTestSets

    trainIdx = training(cv, i);
    testIdx  = test(cv, i);

    Xtrain = X_grade(trainIdx, :);
    Xtest  = X_grade(testIdx, :);

    Ytrain = gradeLabels(trainIdx);

    %%\ PCA fitted using training samples only

    [coeff, scoreTrain, ~, ~, explained, mu] = pca(Xtrain);

    nPC = find(cumsum(explained) >= varianceTarget, 1);

    XtrainPCA = scoreTrain(:,1:nPC);

    XtestCentered = Xtest - mu;
    XtestPCA = XtestCentered * coeff(:,1:nPC);

    %%\ Class imbalance weighting

    nG2 = sum(Ytrain == 'G2');
    nG3 = sum(Ytrain == 'G3');

    costMatrix = [0, nG3/nG2;
        1, 0];

    %%\ Gaussian SVM

    mdlRBF = fitcsvm( ...
        XtrainPCA, ...
        Ytrain, ...
        'KernelFunction','gaussian', ...
        'KernelScale','auto', ...
        'Cost',costMatrix, ...
        'ClassNames',categorical({'G2','G3'}), ...
        'Standardize',true);

    predRBF(testIdx) = string(predict(mdlRBF, XtestPCA));

end
%[text] 
%[text] Evaluation:
%%\ Gaussian SVM performance

predRBF = categorical(predRBF, {'G2','G3'});

cmRBF = confusionmat( ...
    gradeLabels, ...
    predRBF, ...
    'Order',categorical({'G2','G3'}));

sensitivityG2_RBF = ...
    cmRBF(1,1) / sum(cmRBF(1,:));

sensitivityG3_RBF = ...
    cmRBF(2,2) / sum(cmRBF(2,:));

balancedAccuracy_RBF = ...
    mean([sensitivityG2_RBF, sensitivityG3_RBF]);

disp(cmRBF) %[output:691c7277]

fprintf("Gaussian SVM\n"); %[output:63217dab]
fprintf("G2 sensitivity: %.1f%%\n", 100*sensitivityG2_RBF); %[output:84779bfb]
fprintf("G3 sensitivity: %.1f%%\n", 100*sensitivityG3_RBF); %[output:3f210685]
fprintf("Balanced accuracy: %.1f%%\n", ... %[output:group:823e71ae] %[output:3d5f506d]
    100*balancedAccuracy_RBF); %[output:group:823e71ae] %[output:3d5f506d]
%[text] The Gaussian SVM performed only slightly better than the linear SVM, with balanced accuracy increasing from 57.8% to 58.8%. G2 sensitivity remained low at 27.7%, indicating that allowing a nonlinear decision boundary did not substantially improve discrimination between G2 and G3 tumors.
%[text] At the same 70% PCA representation, the Gaussian kernel improved balanced accuracy only marginally.
%%
%[text] 5\. SIMCA analysis
% SIMCA models the minority G2 class only

rng(1);

cvSIMCA = cvpartition(gradeLabels, 'KFold', 5);
%%
%[text] 6\. SIMCA: compare PCA variance thresholds
varianceTargets = [70 80 90];

acceptedSIMCA = false(numel(gradeLabels), numel(varianceTargets));

for i = 1:cvSIMCA.NumTestSets

    trainIdx = training(cvSIMCA, i);
    testIdx  = test(cvSIMCA, i);

    % Use only G2 samples from the training set
    trainG2 = trainIdx & (gradeLabels == 'G2');

    XtrainG2 = X_grade(trainG2, :);
    Xtest    = X_grade(testIdx, :);

    %%\ PCA using training G2 samples only

    [coeffG2, scoreG2, latentG2, ~, explainedG2, muG2] = ...
        pca(XtrainG2);

    cumulativeExplainedG2 = cumsum(explainedG2);

    for j = 1:numel(varianceTargets)

        target = varianceTargets(j);

        nPC = find(cumulativeExplainedG2 >= target, 1);

        coeffUse = coeffG2(:,1:nPC);
        latentUse = latentG2(1:nPC);
        scoreTrainUse = scoreG2(:,1:nPC);

        %%\ Training G2 distances

        T2_train = sum( ...
            (scoreTrainUse.^2) ./ latentUse', ...
            2);

        XtrainReconstructed = ...
            scoreTrainUse * coeffUse';

        residualTrain = ...
            (XtrainG2 - muG2) - XtrainReconstructed;

        Q_train = sum(residualTrain.^2, 2);

        %%\ Define SIMCA limits from training G2

        T2_limit = prctile(T2_train, 95);
        Q_limit  = prctile(Q_train, 95);

        %%\ Evaluate held-out samples

        XtestCentered = Xtest - muG2;

        scoreTest = XtestCentered * coeffUse;

        T2_test = sum( ...
            (scoreTest.^2) ./ latentUse', ...
            2);

        XtestReconstructed = scoreTest * coeffUse';

        residualTest = ...
            XtestCentered - XtestReconstructed;

        Q_test = sum(residualTest.^2, 2);

        %%\ Accept sample as G2 only if both limits are satisfied

        acceptedG2 = ...
            (T2_test <= T2_limit) & ...
            (Q_test <= Q_limit);

        acceptedSIMCA(testIdx,j) = acceptedG2;

    end
end
%%
%[text] 7\. Evaluate SIMCA thresholds
%%\ Evaluate SIMCA thresholds

for j = 1:numel(varianceTargets) %[output:group:8bdfa94e]

    accepted = acceptedSIMCA(:,j);

    trueG2 = gradeLabels == 'G2';
    trueG3 = gradeLabels == 'G3';

    G2accepted = sum(accepted & trueG2);
    G2rejected = sum(~accepted & trueG2);

    G3accepted = sum(accepted & trueG3);
    G3rejected = sum(~accepted & trueG3);

    sensitivityG2 = ...
        G2accepted / sum(trueG2);

    specificityG3 = ...
        G3rejected / sum(trueG3);

    balancedAccuracy = ...
        mean([sensitivityG2, specificityG3]);

    fprintf("\n%d%% PCA variance\n", varianceTargets(j)); %[output:884a8fd8] %[output:33f7a2c8] %[output:16c851bd]
    fprintf("G2 accepted: %d/%d\n", ... %[output:4a3fa52f] %[output:76ab8af9] %[output:42728a14]
        G2accepted, sum(trueG2)); %[output:4a3fa52f] %[output:76ab8af9] %[output:42728a14]
    fprintf("G3 incorrectly accepted: %d/%d\n", ... %[output:0e2c9483] %[output:33cca909] %[output:9d952df3]
        G3accepted, sum(trueG3)); %[output:0e2c9483] %[output:33cca909] %[output:9d952df3]
    fprintf("G2 sensitivity: %.1f%%\n", ... %[output:7da8fb54] %[output:8f11de79] %[output:38481689]
        100*sensitivityG2); %[output:7da8fb54] %[output:8f11de79] %[output:38481689]
    fprintf("G3 specificity: %.1f%%\n", ... %[output:69da2506] %[output:547488de] %[output:28a9e56d]
        100*specificityG3); %[output:69da2506] %[output:547488de] %[output:28a9e56d]
    fprintf("Balanced accuracy: %.1f%%\n", ... %[output:5832bae0] %[output:1b9d8380] %[output:743b4c6f]
        100*balancedAccuracy); %[output:5832bae0] %[output:1b9d8380] %[output:743b4c6f]

end %[output:group:8bdfa94e]
%[text] Using SIMCA with PCA models built on the training G2 samples and 95th-percentile T² and Q limits resulted in all held-out samples being rejected, including all true G2 tumors. This gave 0% G2 sensitivity and 50% balanced accuracy, showing that the acceptance thresholds were too strict for unseen samples. Therefore, the SIMCA limits need to be calibrated using held-out G2 behavior rather than the same samples used to build the class model.
%%
%[text] 8\. Leave-one-out calibrated SIMCA: PCA sensitivity analysis
%%\ Leave-one-out calibrated SIMCA: compare PCA thresholds

rng(1);

varianceTargets = [70 80 90];

cvSIMCA = cvpartition(gradeLabels, 'KFold', 5);

% Rows = patients, columns = PCA variance thresholds
acceptedSIMCA_calibrated = ...
    false(numel(gradeLabels), numel(varianceTargets));

for fold = 1:cvSIMCA.NumTestSets

    %%\ Outer training and test sets

    trainIdx = training(cvSIMCA, fold);
    testIdx  = test(cvSIMCA, fold);

    % Only G2 samples are used to build the SIMCA model
    trainG2 = trainIdx & (gradeLabels == 'G2');

    XtrainG2 = X_grade(trainG2, :);
    Xtest    = X_grade(testIdx, :);

    nG2train = size(XtrainG2, 1);

    %%\ Test each PCA variance threshold

    for v = 1:numel(varianceTargets)

        varianceTarget = varianceTargets(v);

        %%\ Leave-one-out calibration within training G2

        T2_LOO = zeros(nG2train,1);
        Q_LOO  = zeros(nG2train,1);

        for j = 1:nG2train

            looTrain = true(nG2train,1);
            looTrain(j) = false;

            X_loo_train = XtrainG2(looTrain,:);
            X_loo_test  = XtrainG2(j,:);

            % PCA on remaining G2 samples
            [coeff_loo, ~, latent_loo, ~, explained_loo, mu_loo] = ...
                pca(X_loo_train);

            nPC_loo = find( ...
                cumsum(explained_loo) >= varianceTarget, ...
                1);

            coeff_loo = coeff_loo(:,1:nPC_loo);
            latent_loo = latent_loo(1:nPC_loo);

            % Project left-out G2 patient
            X_loo_centered = X_loo_test - mu_loo;
            score_loo = X_loo_centered * coeff_loo;

            % T2 distance
            T2_LOO(j) = sum( ...
                (score_loo.^2) ./ latent_loo');

            % Q residual distance
            X_loo_reconstructed = ...
                score_loo * coeff_loo';

            residual_loo = ...
                X_loo_centered - X_loo_reconstructed;

            Q_LOO(j) = sum(residual_loo.^2,2);

        end

        %%\ Limits based on unseen G2 behavior

        T2_limit = prctile(T2_LOO,95);
        Q_limit  = prctile(Q_LOO,95);

        %%\ Final G2 PCA model for this outer fold

        [coeffG2, ~, latentG2, ~, explainedG2, muG2] = ...
            pca(XtrainG2);

        nPC = find( ...
            cumsum(explainedG2) >= varianceTarget, ...
            1);

        coeffG2 = coeffG2(:,1:nPC);
        latentG2 = latentG2(1:nPC);

        %%\ Evaluate completely held-out patients

        XtestCentered = Xtest - muG2;
        scoreTest = XtestCentered * coeffG2;

        T2_test = sum( ...
            (scoreTest.^2) ./ latentG2', ...
            2);

        XtestReconstructed = ...
            scoreTest * coeffG2';

        residualTest = ...
            XtestCentered - XtestReconstructed;

        Q_test = sum(residualTest.^2,2);

        %%\ Accept as G2

        acceptedG2 = ...
            (T2_test <= T2_limit) & ...
            (Q_test <= Q_limit);

        acceptedSIMCA_calibrated(testIdx,v) = acceptedG2;

    end
end
%[text] Evaluate the calibrated SIMCA thresholds:
%%\ Evaluate calibrated SIMCA thresholds

trueG2 = gradeLabels == 'G2';
trueG3 = gradeLabels == 'G3';

for v = 1:numel(varianceTargets) %[output:group:78e086ae]

    accepted = acceptedSIMCA_calibrated(:,v);

    G2accepted = sum(accepted & trueG2);
    G3accepted = sum(accepted & trueG3);

    sensitivityG2 = ...
        G2accepted / sum(trueG2);

    specificityG3 = ...
        sum(~accepted & trueG3) / sum(trueG3);

    balancedAccuracy = ...
        mean([sensitivityG2, specificityG3]);

    fprintf("\nCalibrated SIMCA - %d%% PCA variance\n", ... %[output:93f990cc] %[output:27340807] %[output:65cc67d0]
        varianceTargets(v)); %[output:93f990cc] %[output:27340807] %[output:65cc67d0]

    fprintf("G2 accepted: %d/%d\n", ... %[output:636c3bd7] %[output:35b23125] %[output:2d60b208]
        G2accepted, sum(trueG2)); %[output:636c3bd7] %[output:35b23125] %[output:2d60b208]

    fprintf("G3 incorrectly accepted: %d/%d\n", ... %[output:2df79de6] %[output:335a7058] %[output:0262515b]
        G3accepted, sum(trueG3)); %[output:2df79de6] %[output:335a7058] %[output:0262515b]

    fprintf("G2 sensitivity: %.1f%%\n", ... %[output:9c8476e9] %[output:9c4b85fc] %[output:1bc85311]
        100*sensitivityG2); %[output:9c8476e9] %[output:9c4b85fc] %[output:1bc85311]

    fprintf("G3 specificity: %.1f%%\n", ... %[output:35e61fd0] %[output:18c9bc19] %[output:537b03d7]
        100*specificityG3); %[output:35e61fd0] %[output:18c9bc19] %[output:537b03d7]

    fprintf("Balanced accuracy: %.1f%%\n", ... %[output:0321000e] %[output:629c5c16] %[output:8fb42e02]
        100*balancedAccuracy); %[output:0321000e] %[output:629c5c16] %[output:8fb42e02]

end %[output:group:78e086ae]
%[text] Calibrated SIMCA performance was consistent across PCA variance thresholds of 70%, 80%, and 90%. Although most G2 tumors were accepted by the class model, approximately 89% of G3 tumors were also accepted. Balanced accuracy remained close to 50% for all thresholds, indicating that the poor discrimination was robust to the number of retained principal components. Therefore, G2 tumors do not appear to form a sufficiently distinct expression class from G3 tumors for SIMCA classification.
%[text] 
%[text] 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:2afa8693]
%   data: {"dataType":"text","outputData":{"text":"G2 samples: 47\n","truncated":false}}
%---
%[output:0a5dfc72]
%   data: {"dataType":"text","outputData":{"text":"G3 samples: 339\n","truncated":false}}
%---
%[output:118c24c2]
%   data: {"dataType":"text","outputData":{"text":"Expression matrix for grade analysis: 386 patients x 11838 genes\n","truncated":false}}
%---
%[output:48899e97]
%   data: {"dataType":"text","outputData":{"text":"\n70% PCA variance\n","truncated":false}}
%---
%[output:5d26a749]
%   data: {"dataType":"text","outputData":{"text":"G2 sensitivity: 27.7%\n","truncated":false}}
%---
%[output:175bf82a]
%   data: {"dataType":"text","outputData":{"text":"G3 sensitivity: 87.9%\n","truncated":false}}
%---
%[output:7b5ad979]
%   data: {"dataType":"text","outputData":{"text":"Balanced accuracy: 57.8%\n","truncated":false}}
%---
%[output:09223198]
%   data: {"dataType":"text","outputData":{"text":"    13    34\n    41   298\n\n","truncated":false}}
%---
%[output:9257498d]
%   data: {"dataType":"text","outputData":{"text":"\n80% PCA variance\n","truncated":false}}
%---
%[output:60d2921e]
%   data: {"dataType":"text","outputData":{"text":"G2 sensitivity: 10.6%\n","truncated":false}}
%---
%[output:12623b04]
%   data: {"dataType":"text","outputData":{"text":"G3 sensitivity: 93.5%\n","truncated":false}}
%---
%[output:371d078b]
%   data: {"dataType":"text","outputData":{"text":"Balanced accuracy: 52.1%\n","truncated":false}}
%---
%[output:5407ebf5]
%   data: {"dataType":"text","outputData":{"text":"     5    42\n    22   317\n\n","truncated":false}}
%---
%[output:4fe855a2]
%   data: {"dataType":"text","outputData":{"text":"\n90% PCA variance\n","truncated":false}}
%---
%[output:00393f65]
%   data: {"dataType":"text","outputData":{"text":"G2 sensitivity: 10.6%\n","truncated":false}}
%---
%[output:2d490319]
%   data: {"dataType":"text","outputData":{"text":"G3 sensitivity: 95.6%\n","truncated":false}}
%---
%[output:824c2b89]
%   data: {"dataType":"text","outputData":{"text":"Balanced accuracy: 53.1%\n","truncated":false}}
%---
%[output:16360e11]
%   data: {"dataType":"text","outputData":{"text":"     5    42\n    15   324\n\n","truncated":false}}
%---
%[output:691c7277]
%   data: {"dataType":"text","outputData":{"text":"    13    34\n    35   304\n\n","truncated":false}}
%---
%[output:63217dab]
%   data: {"dataType":"text","outputData":{"text":"Gaussian SVM\n","truncated":false}}
%---
%[output:84779bfb]
%   data: {"dataType":"text","outputData":{"text":"G2 sensitivity: 27.7%\n","truncated":false}}
%---
%[output:3f210685]
%   data: {"dataType":"text","outputData":{"text":"G3 sensitivity: 89.7%\n","truncated":false}}
%---
%[output:3d5f506d]
%   data: {"dataType":"text","outputData":{"text":"Balanced accuracy: 58.7%\n","truncated":false}}
%---
%[output:884a8fd8]
%   data: {"dataType":"text","outputData":{"text":"\n70% PCA variance\n","truncated":false}}
%---
%[output:4a3fa52f]
%   data: {"dataType":"text","outputData":{"text":"G2 accepted: 0\/47\n","truncated":false}}
%---
%[output:0e2c9483]
%   data: {"dataType":"text","outputData":{"text":"G3 incorrectly accepted: 0\/339\n","truncated":false}}
%---
%[output:7da8fb54]
%   data: {"dataType":"text","outputData":{"text":"G2 sensitivity: 0.0%\n","truncated":false}}
%---
%[output:69da2506]
%   data: {"dataType":"text","outputData":{"text":"G3 specificity: 100.0%\n","truncated":false}}
%---
%[output:5832bae0]
%   data: {"dataType":"text","outputData":{"text":"Balanced accuracy: 50.0%\n","truncated":false}}
%---
%[output:33f7a2c8]
%   data: {"dataType":"text","outputData":{"text":"\n80% PCA variance\n","truncated":false}}
%---
%[output:76ab8af9]
%   data: {"dataType":"text","outputData":{"text":"G2 accepted: 0\/47\n","truncated":false}}
%---
%[output:33cca909]
%   data: {"dataType":"text","outputData":{"text":"G3 incorrectly accepted: 0\/339\n","truncated":false}}
%---
%[output:8f11de79]
%   data: {"dataType":"text","outputData":{"text":"G2 sensitivity: 0.0%\n","truncated":false}}
%---
%[output:547488de]
%   data: {"dataType":"text","outputData":{"text":"G3 specificity: 100.0%\n","truncated":false}}
%---
%[output:1b9d8380]
%   data: {"dataType":"text","outputData":{"text":"Balanced accuracy: 50.0%\n","truncated":false}}
%---
%[output:16c851bd]
%   data: {"dataType":"text","outputData":{"text":"\n90% PCA variance\n","truncated":false}}
%---
%[output:42728a14]
%   data: {"dataType":"text","outputData":{"text":"G2 accepted: 0\/47\n","truncated":false}}
%---
%[output:9d952df3]
%   data: {"dataType":"text","outputData":{"text":"G3 incorrectly accepted: 0\/339\n","truncated":false}}
%---
%[output:38481689]
%   data: {"dataType":"text","outputData":{"text":"G2 sensitivity: 0.0%\n","truncated":false}}
%---
%[output:28a9e56d]
%   data: {"dataType":"text","outputData":{"text":"G3 specificity: 100.0%\n","truncated":false}}
%---
%[output:743b4c6f]
%   data: {"dataType":"text","outputData":{"text":"Balanced accuracy: 50.0%\n","truncated":false}}
%---
%[output:93f990cc]
%   data: {"dataType":"text","outputData":{"text":"\nCalibrated SIMCA - 70% PCA variance\n","truncated":false}}
%---
%[output:636c3bd7]
%   data: {"dataType":"text","outputData":{"text":"G2 accepted: 41\/47\n","truncated":false}}
%---
%[output:2df79de6]
%   data: {"dataType":"text","outputData":{"text":"G3 incorrectly accepted: 303\/339\n","truncated":false}}
%---
%[output:9c8476e9]
%   data: {"dataType":"text","outputData":{"text":"G2 sensitivity: 87.2%\n","truncated":false}}
%---
%[output:35e61fd0]
%   data: {"dataType":"text","outputData":{"text":"G3 specificity: 10.6%\n","truncated":false}}
%---
%[output:0321000e]
%   data: {"dataType":"text","outputData":{"text":"Balanced accuracy: 48.9%\n","truncated":false}}
%---
%[output:27340807]
%   data: {"dataType":"text","outputData":{"text":"\nCalibrated SIMCA - 80% PCA variance\n","truncated":false}}
%---
%[output:35b23125]
%   data: {"dataType":"text","outputData":{"text":"G2 accepted: 42\/47\n","truncated":false}}
%---
%[output:335a7058]
%   data: {"dataType":"text","outputData":{"text":"G3 incorrectly accepted: 302\/339\n","truncated":false}}
%---
%[output:9c4b85fc]
%   data: {"dataType":"text","outputData":{"text":"G2 sensitivity: 89.4%\n","truncated":false}}
%---
%[output:18c9bc19]
%   data: {"dataType":"text","outputData":{"text":"G3 specificity: 10.9%\n","truncated":false}}
%---
%[output:629c5c16]
%   data: {"dataType":"text","outputData":{"text":"Balanced accuracy: 50.1%\n","truncated":false}}
%---
%[output:65cc67d0]
%   data: {"dataType":"text","outputData":{"text":"\nCalibrated SIMCA - 90% PCA variance\n","truncated":false}}
%---
%[output:2d60b208]
%   data: {"dataType":"text","outputData":{"text":"G2 accepted: 42\/47\n","truncated":false}}
%---
%[output:0262515b]
%   data: {"dataType":"text","outputData":{"text":"G3 incorrectly accepted: 300\/339\n","truncated":false}}
%---
%[output:1bc85311]
%   data: {"dataType":"text","outputData":{"text":"G2 sensitivity: 89.4%\n","truncated":false}}
%---
%[output:537b03d7]
%   data: {"dataType":"text","outputData":{"text":"G3 specificity: 11.5%\n","truncated":false}}
%---
%[output:8fb42e02]
%   data: {"dataType":"text","outputData":{"text":"Balanced accuracy: 50.4%\n","truncated":false}}
%---
