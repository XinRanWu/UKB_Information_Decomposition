% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%  Scripts for Luppi et al.,
%%  Reduced emergent character of neural dynamics in patients with a disrupted connectome
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% This workflow runs through the computation of Emergence Capacity and
% ignition-based hierarchy from functional data, and then uses the
% structural connectomes to fit a dynamic mean-field model and simulate
% brain activity timeseries

% Requires FastDMF (Luppi et al., 2022 Comms Biology; https://www.gitlab.com/concog/fastdmf)
% and PhiID (Mediano, Rosas et al., 2021 arXiv; Luppi et al., 2022 Nature Neuroscience) 
% software packages to be on the MATLAB path: 
% see the corresponding documentation for each, and their dependencies.
% Also requires the Brain Connectivity Toolbox (Rubinov and Sporns, 2010, 2011)
% on the MATLAB path.

%  Warning: the workflow (especially PhiID) can be very time-consuming, in the order
%  of several hours (possibly >12) per subject.
%  Parallelisation should be considered if possible.


%%  Load inputs (empirical data):
clear all

load('Data.mat')
%  User should prepare Data.mat to contain the following:
%  BOLD_data: a cell array with numSubs cells, each cell containing a subject's NxT timeseries
%  deconvolved_BOLD_data: same as BOLD_data but after HRF deconvolution, e.g. using the toolbox from Wu et al (2013), Med Image Analysis.
%  individual_DTI, a cell array with numSubs cells, each cell containing a subject's NxN structural conenctome
%  T:  your number of TRs
%  N: your number of regions
%  TR: repetition time of your data in seconds
%  lowend, highend: limits of bandpass filter, in Hz
%  numSubs: Number of subjects; should be the same in both BOLD and DTI

%  Parameters for functional measures
MMI_CCS = 'CCS'; %  type of information decomposition (definition of redundancy)
discreteYN = true; %  whether discrete or continuous information decompositon
tau = 1; %  time-scale for information decomposition
apply_bandpass_YN = true; %  set to false if providing narrowband-filtered data for this
numIters = 10  %how many simulations to run


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%  RUN WORKFLOW
%  Obtain emergence capacity, hierarchy, and controllability measures for each subject
for sub = 1:numSubs % Number of subjects should be the same in both BOLD and DTI
    
    [hierarchy(sub)] = compute_ignition(BOLD_data{sub}, TR, apply_bandpass_YN);
    
    %  Warning: this is very time-consuming!
    %  It scales quadratically with N, so with 200+ ROIs, this can take ~10h per subject
    [EmergenceCapacity(sub)] = PhiID_EmergenceCapacity_from_timeseries(deconvolved_BOLD_data{sub}, MMI_CCS, discreteYN, tau);
    
    %  Average and modal controllability
    ave_ctrl(sub) = mean(ave_control(individual_DTI{sub}));
    modal_ctrl(sub) = mean(modal_control(individual_DTI{sub}));
    
end

%  Set general parameters for dynamic FC; note that these are hardcoded
window_size = 30;
sliding_by = 3;

%  Obtain FCD from the empirical timeseries, for fitting
%  also combine all FCDs into a group-level, provided in the same format
groupFCD{1} = [];
for sub = 1:numSubs
    
    [empiricalFCD{sub}] = DMF_histogram(BOLD_data{sub}(1:N, 1:T), sliding_by, window_size);
    groupFCD{1} = [groupFCD{1}; empiricalFCD{sub}];
    
end

%  Simulate data using FastDMF with different G parameter values to find the best one
[simFCD, simFiring] = DMF_simulate_FCD(individual_DTI, T, N, TR, lowend, highend, numIters);

%  Now find best-fitting G: either from fitting the FCD, or from looking
%  where the firing rate becomes uncontrolled;
%  In this example we follow the FCD method
[BestG, BestFit] = DMF_find_BestG_individual_DTI(empiricalFCD, simFCD);

%  Now simulate the data with FastDMF using the best-fitting G for each subject; one
%  simulation per subject (more possible)
[simulated_timeseries] = DMF_simulate_timeseries(BestG, individual_DTI, T, N, TR, lowend, highend, 1);

%  Now we can compute the measures on simulated data
for sub = 1:numSubs
    [hierarchy_SIM(sub)] = compute_ignition(simulated_timeseries{sub}, TR, apply_bandpass_YN);

    %  Warning: this is very time-consuming!
    %  It scales quadratically with N, so with 200+ ROIs, this can take ~10h per subject
    [EmergenceCapacity_SIM(sub)] = PhiID_EmergenceCapacity_from_timeseries(simulated_timeseries{sub}, MMI_CCS, discreteYN, tau);
end


%%  Note: 
% This workflow also works for group-level simulation: in that case
% provide individualDTI as a single consensus-connectome (only one cell in the cell array)
%  e.g. using the method of Wang et al (2019) Science Advances,
% and perform the fitting based on this connectome, using as fitting target the groupFCD computed above



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%  FUNCTION DEFINITIONS
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function [EmergenceCapacity] = PhiID_EmergenceCapacity_from_timeseries(ts, MMI_CCS, discreteYN, tau)
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
%  Compute Emergence Capacity from integrated information decomposition (PhiID)
% 
%  Input arguments:
%  ts (matrix of size NxT): timeseries data with N regions and T timepoints
%  MMI_CCS (string): whether to use minimum mutual information or common change in surprisal definitions of redundancy; options are 'ccs'  or 'mmi'
%  discreteYN (boolean): whether to use the discrete version (true) or the continuous one (false)
%  tau (integer): time-step, in terms of number of time-points; default is 1
% 
%  Output:
%  EmergenceCapacity for the whole brain (average of all pairwise)
% % % % % % % % % % % % % % % % % % % % % % % % % % % 

if not(exist('tau', 'var'))
    tau=1; % 1 time-step is default
end

disp(['Running PhiID using ', MMI_CCS])

[num_nodes, num_timepoints] = size(ts);


for row = 1:num_nodes
    for col = 1:num_nodes
        
        % Automatically set to zero any case where all elements of a timeseries
        % are zero or if NaNs found, or if the two timeseries are perfectly identical
        % (indicating bad data, which should not happen)
        if row == col ...
                | sum(abs(ts(row, :))) == 0 ...
                | sum(abs(ts(col, :))) == 0 ...
                | isnan(sum(abs(ts(row, :))) + sum(abs(ts(col, :)))) ...
                | sum(sum(abs( ts(row, :) - ts(col, :)))) == 0
            
            PairwisePsi(row,col) = 0;

        else
            
            if discreteYN == true
                atoms = PhiIDFullDiscrete([ts(row, :); ts(col, :)], tau, MMI_CCS);
            else
                atoms = PhiIDFull([ts(row, :); ts(col, :)], tau, MMI_CCS);
            end
                
        PairwisePsi(row,col) = atoms.sts + atoms.stx + atoms.sty + atoms.str; % following PhiID formula
        end

        clear atoms
    end
end

% Finally take the mean across all pairwise values of emergence capacity
EmergenceCapacity = mean(mean(PairwisePsi));

end




function[hierarchy] = compute_ignition(ts, TR, apply_bandpass_YN)
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
%  Compute hierarchy from intrinsic ignition as per Deco et al 2017 Neuron
%  ts: (NxT): the BOLD signal data
%  TR: seconds
%  apply_bandpass_YN (boolean): true if data need to be filtered (since narrowband
%  filter required); set to false if the data are already filtered in the
%  required band of 0.04-0.07 Hz
% % % % % % % % % % % % % % % % % % % % % % % % % % % 

[N, T] = size(ts);

% Step 1 - narrowband filter
if apply_bandpass_YN == true
    [ts] = bandpass_timeseries(ts, TR, 0.04, 0.07);
end

% Step 2 - Z-score
for roi = 1:N
    Zscored_timeseries(roi,:) = zscore(ts(roi,:));
end

% Step 3 - Binarise (threshold  Z>1)
binary_timeseries = zeros(size(Zscored_timeseries));
for roi = 1:N
    for t = 2:T
        
        % Set to 1 if above threshold (Zscore > 1) and crossing the
        % threshold from below (i.e. the previous one is below threshold)
        if Zscored_timeseries(roi,t) > 1 & Zscored_timeseries(roi,t-1) < 1
            binary_timeseries(roi,t) = 1;
        end
    end
end


% Step 4 - Intrinsic Ignition
window = 3; % Num of TRs within which to compute ignition (3 plus the current one)
for roi = 1:N
    for t = 1:T
        real_window = min(window, T - t); % window can't go beyond the number of TRs available
        
        if binary_timeseries(roi,t) == 1;
            
            % Find which regions are active in the slice of binary timeseries for the window considered
            synch_bin = double(mean(binary_timeseries(:, t : t+real_window) ,2) > 0);
            synch_mat = synch_bin * synch_bin';
            
            % Take size of largest component as measure of intrinsic ignition
            [components, sizes] = get_components(synch_mat);  % BCT function         
            ignition(roi,t) = max(sizes) ./ N;
            
        end
    end
end

% Step 5 - Hierarchy is the standard deviation (variability) across regions
% of the average size of ignition induced by each region
for roi = 1:N
    avg_ignition(roi,1) = mean(ignition(roi,find(ignition(roi,:) ~=0)));
end
hierarchy = std(avg_ignition); % std across nodes: hierarchy (Deco...Tagliazucchi 2017 eNeuro)

end % EOF





function[simFCD, simFiring] = DMF_simulate_FCD(individual_DTI, T, N, TR, lowend, highend, numIters);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%  Takes the DTI from the specified condition and runs simulations for
%  each subject's DTI, across a variety of G values, to obtain simulated FCD
% 
%  Inputs:
%  individual_DTI: cell array with one cell per subject, each cell being an NxN matrix of SC
%  T: number of timepoints (BOLD volumes)
%  N: number of regions (must be the same for all subjects, so avoid missing
%  regions)
%  TR: temporal resolution (repetition time) in seconds
%  lowend, highend: low and high ends for the bandpass filter, in Hz
%  numIters (integer): how many times to simulate each subject, to obtain better estimate of their FCD
% 
%  Outputs:
%  simFCD: cell array of simulated FCD for each subject (cell), with each
%  cell containing another cell array of size numIters x g_values
%  simFiring: same organisation as simFCD, but contains the simulated mean
%  firing rate for each simulation
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

g_vec = 0.1 : 0.1 : 3; % values of the global coupling parameter G to consider
numG = numel(g_vec);

% First get default DMF params (remove those we need to set, just to be sure)
params = DefaultParams();
rmfield(params, 'C');
rmfield(params, 'J');
rmfield(params, 'G');

% Then replace those that need specific values for this simulation
params.N = N;
params.T = T;

% Get filter params
params.TR = TR;
params.lowend = lowend;
params.highend = highend;


% initial and final extra volumes to discard, to minimise effects of
% filtering and equilibration
params.initial_TRs =  floor(80 ./ params.TR); % heuristic
params.final_TRs = 20;

% Standard parameters from Deco 2018
params.sliding_by = 3;
params.window_size = 30;

% Time points (milliseconds); we add a number of TRs equal to the one we remove, so the total we keep is the same as T
params.tmax = (params.T * params.TR * 1000) + (params.initial_TRs * params.TR * 1000) + (params.final_TRs * params.TR * 1000);

%  Parallel computation parameters
params.batch_size = 5000;
nb_steps = params.tmax;


% Assign each individual's params (all identical except for the SC and those derived from it, alphas and strength)
% This way it is possible to use PARFOR to speed things up
numSub = numel(individual_DTI);
for sub = 1:numSub
    
    SC =  individual_DTI{sub}(1:params.N, 1:params.N);
    for g_num = 1:numel(g_vec)
        
        individual_params{sub, g_num} = params;       
        individual_params{sub, g_num}.G = g_vec(g_num);
        individual_params{sub, g_num}.C = SC/max(SC (:))*0.2;
        individual_params{sub, g_num}.alphas = ones(size(individual_params{sub, g_num}.C,1),1).*1.5; %  parameter of the feedback inhibitory control
        individual_params{sub, g_num}.stren  = sum(individual_params{sub, g_num}.C)'./2; %  node strength
        individual_params{sub, g_num}.J      = 0.75*individual_params{sub, g_num}.G*sum(individual_params{sub, g_num}.C, 1)' + 1;
        
    end
    clear SC
end


% Note: the code is organised so that this could be replaced by a PARFOR
% loop if required to speed up processing;
% uncomment the following three lines and comment out the one immediately
% after, to do so

%  delete(gcp('nocreate'))
%  parpool(30)
%  parfor g_num = 1:numG
for g_num = 1:numG
    
    for sub = 1:numSub
        disp(['Simulating sub #', num2str(sub), '/', num2str(numSub) ])
        
        % Repeat multiple times if better estimation desired
        for iter = 1:numIters
                        
            simulated_BOLD = DMF(individual_params{sub, g_num}, individual_params{sub, g_num}.tmax);
            
            % also get the FIRING rate for alternative decision of where to
            % fit the model
            
            [simulated_timeseries] = bandpass_timeseries(simulated_BOLD, individual_params{sub, g_num}.TR, individual_params{sub, g_num}.lowend, individual_params{sub, g_num}.highend);
            simulated_timeseries = simulated_timeseries(:, 1+individual_params{sub, g_num}.initial_TRs : end-individual_params{sub, g_num}.final_TRs);
            
            [cotsamplingsim] = DMF_histogram(simulated_timeseries, individual_params{sub}.sliding_by, individual_params{sub, g_num}.window_size);
            myFCD{g_num}{iter, sub} = cotsamplingsim; % organised like this so it can be run with PARFOR
            
            
            % Now simulate firing rate: we could do both in the same round,
            % but that incurs overhead costs and there is little point in
            % getting such long timeseries for the firing rate
            simulated_firing = DMF(individual_params{sub, g_num}, 21000, 'rate'); % 21000 seems good approximation for what one gets with longer timeseries
            myAvgFiring{g_num}{iter, sub} = mean(median(simulated_firing(:, 1001: end),2));
            
        end
        
        disp(['Finished sub #', num2str(sub), ' for G = ', num2str(g_vec(g_num))])
        
    end
end

% Reshape for saving
for sub = 1:numSub
    for g_num = 1:numG
        for iter = 1:numIters
            
            simFCD{sub}{iter, g_num} = myFCD{g_num}{iter, sub};
            simFiring{sub}{iter, g_num} = myAvgFiring{g_num}{iter, sub};
            
        end
    end
end

end % EOF



function[cotsamplingsim] = DMF_histogram(ts, sliding_by, window_size);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%  Obtain values of functional connectivity between different time-points;
%  useful for FCD fitting of mean-field model
% 
%  Inputs:
%  ts: NxT timeseries data
%  sliding_by: (integer): how many volumes you slide by (greater number means less overlap)
%  window_size: (integer): size of the window for computing dynamic functional connectivity
% 
%  Outputs:
%  cotsamplingsim is a vector of all the values of the FCD, to be used for KS distance fitting of the DMF
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

[N,T] = size(ts);
Isubdiag = find(tril(ones(N),-1));

ts=ts';
overlap = window_size-sliding_by;

win_start = 0:window_size-overlap:T-window_size-1;
nwins = length(win_start);
fcd = zeros(length(Isubdiag),nwins);

for i=1:nwins
    tmp = ts(win_start(i)+1:win_start(i)+window_size+1,:);
    cormat = corrcoef(tmp);
    fcd(:,i) = cormat(Isubdiag);
end
FCD_mat = tril(corrcoef(fcd))-eye(nwins);

cotsamplingsim=FCD_mat(not(triu(true(size(FCD_mat)))));

end % EOF





function[BestG, BestFit] = DMF_find_BestG_individual_DTI(empiricalFCD, simFCD)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%  Find the G value that gives the best FCD fit between empirical and simulated data
%  Suitable for both individual-level simulations and fitting, or group-level (by treating as one "subject")
% 
%  Inputs:
%  empiricalFCD: the empirical FCD of each subject; a cell array with one cell per subject;
%  simFCD: the FCD obtained from simulation; a cell array with one cell per subject; each cell
%  contains another cell array, of size number of iterations x number of G values
%  Example: to access the FCD obtained from a given subject, G value, and iteration, you would use FCD{sub}{iter, g_num}
% 
%  NB: ensure that the empirical and simulated FCD were obtained from the
%  same sliding_by and window_size parameters and TR/filtering
% 
%  Outputs:
%  BestG: the G value that provides the best FCD fit (minimum KS distance)
%  between real and simulated data, averating across iterations, for each subject
%  BestFit: the actual fit value
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Default values of the G global coupling parameter to consider, if not provided
if nargin < 3
    g_vec = 0.1 : 0.1 : 3;
end

% Now compare each subject against the simulated FCD obtained from their own DTI,
% for all G values (averaging over all iterations)

numSub = numel(simFCD);
assert(numel(empiricalFCD) == numel(simFCD)) % there should be as many simulated subjects as empirical ones

for sub = 1:numSub
    [numIter, numG] = size(simFCD{sub});
    assert(numel(g_vec) == numG) % check if the vector of G values is correct

    for g = 1:numG
        for iter = 1:numIter
            
            [~, ~, KS(iter)] = kstest2( empiricalFCD{sub}(:), simFCD{sub}{iter,g}(:) );
            
        end
        
        FCD_fit(sub, g) = mean(KS);
        clear KS
        
    end
    
    % KS fit and corresponding G
    [BestFit(sub,1), BestIdx(sub,1)] = min(FCD_fit(sub,:));
    BestG(sub,1) = g_vec(BestIdx(sub,1));
    
end

end % EOF



function[simulated_timeseries] = DMF_simulate_timeseries(BestG, individual_DTI, T, N, TR, lowend, highend, numIters);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%  Takes the DTI from the specified condition and runs simulations for
%  each subject's DTI, using the best-fitting G
% 
%  Inputs:
%  BestG: the G values optimised or chosen for each subject
%  individual_DTI: cell array with one cell per subject, each cell being an NxN matrix of SC
%  T: number of timepoints (BOLD volumes)
%  N: number of regions (must be the same for all subjects, so avoid missing
%  regions)
%  TR: temporal resolution (repetition time) in seconds
%  lowend, highend: low and high ends for the bandpass filter, in Hz
%  numIters (integer): how many times to simulate each subject, to obtain better estimate of their FCD
% 
%  Outputs:
%  simulated_timeseries: cell array with one cell per subject, containing
%  one simulation per iteration
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% First get default DMF params
params = DefaultParams();
rmfield(params, 'C');
rmfield(params, 'J');
rmfield(params, 'G');

% Then replace those that need specific values for this simulation
params.N = N;
params.T = T;

% Get filter params
params.TR = TR;
params.lowend = lowend;
params.highend = highend;


% initial and final extra volumes to discard, to minimise effects of
% filtering and equilibration
params.initial_TRs =  floor(80 ./ params.TR);
params.final_TRs = 20;

% Standard parameters from Deco 2018
params.sliding_by = 3;
params.window_size = 30;

% Time points (milliseconds); we add a number of TRs equal to the one we remove, so the total we keep is the same as T
params.tmax = (params.T * params.TR * 1000) + (params.initial_TRs * params.TR * 1000) + (params.final_TRs * params.TR * 1000);  

%  Parallel computation parameters
params.batch_size = 5000;
nb_steps = params.tmax;


% Assign each individual's params (all identical except for the SC and those derived from it, alphas and strength)
% This way it is possible to use PARFOR to speed things up
numSub = numel(individual_DTI);
for sub = 1:numSub
    
    SC =  individual_DTI{sub}(1:params.N, 1:params.N);
    
    individual_params{sub, 1} = params;
    individual_params{sub, 1}.G = BestG(sub);
    individual_params{sub, 1}.C = SC/max(SC (:))*0.2;
    individual_params{sub, 1}.alphas    = ones(size(individual_params{sub, 1}.C,1),1).*1.5; %  parameter of the feedback inhibitory control
    individual_params{sub, 1}.stren     = sum(individual_params{sub, 1}.C)'./2; %  node strength
    individual_params{sub, 1}.J         = 0.75*individual_params{sub, 1}.G*sum(individual_params{sub, 1}.C, 1)' + 1;
    
    clear SC
end


% Note: the code is organised so that this could be replaced by a PARFOR
% loop if required to speed up processing;
% uncomment the following two lines and comment out the one immediately
% after, to do so
% parpool(20)
% parfor g_num = 1:numG


for sub = 1:numSub
    disp(['Simulating sub #', num2str(sub), '/', num2str(numSub) ])
    
    % Repeat multiple times if better estimation desired
    for iter = 1:numIters
        
        simulated_BOLD = DMF(individual_params{sub, 1}, individual_params{sub, 1}.tmax);
        
        % also get the FIRING rate for alternative decision of where to
        % fit the model
        
        [simulated_BOLD] = bandpass_timeseries(simulated_BOLD, individual_params{sub, 1}.TR, individual_params{sub, 1}.lowend, individual_params{sub, 1}.highend);
        simulated_timeseries{sub, iter} = simulated_BOLD(:, 1+individual_params{sub, 1}.initial_TRs : end-individual_params{sub, 1}.final_TRs);
        clear simulated_BOLD
        
        
    end
    
    disp(['Finished sub #', num2str(sub)])
end

end % EOF


function[filtered_timeseries] = bandpass_timeseries(ts, TR, lowend, highend)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%  Bandpass BOLD signal timeseries
%  Inputs:
%  ts: NxT BOLD signal data
%  TR = repetition time in seconds
%  lowend = lowest frequency in Hz; default 0.008
%  highend = highest frequency in Hz; default 0.09
% 
%  Outputs:
%  filtered_timeseries: NxT bandpass-filtered timeseries
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

if nargin < 4
    lowend = 0.008;
    highend = 0.09;
end

%  FILTER SETTINGS
fnq=1/(2*TR);                 %  Nyquist frequency
Wn=[lowend/fnq highend/fnq];  %  butterworth bandpass non-dimensional frequency
k=2;                          %  2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   %  construct the filter

for ROI=1:size(ts,1)
    filtered_timeseries(ROI,:) =filtfilt(bfilt,afilt,ts(ROI, :)); %  Band pass filter
end
end% EOF


function [values] = ave_control(A)
%  FUNCTION:
%          Returns values of AVERAGE CONTROLLABILITY for each node in a
%          network, given the adjacency matrix for that network. Average
%          controllability measures the ease by which input at that node can
%          steer the system into many easily-reachable states.
% 
%  INPUT:
%          A is the structural (NOT FUNCTIONAL) network adjacency matrix,
%  	      such that the simple linear model of dynamics outlined in the
% 	      reference is an accurate estimate of brain state fluctuations.
% 	      Assumes all values in the matrix are positive, and that the
% 	      matrix is symmetric.
% 
%  OUTPUT:
%          Vector of average controllability values for each node
% 
%  Bassett Lab, University of Pennsylvania, 2016.
%  Reference: Gu, Pasqualetti, Cieslak, Telesford, Yu, Kahn, Medaglia,
%             Vettel, Miller, Grafton & Bassett, Nature Communications
%             6:8414, 2015.

A = A./(1+svds(A,1));     %  Matrix normalization
[U, T] = schur(A,'real'); %  Schur stability
midMat = (U.^2)';
v = diag(T);
P = repmat(diag(1 - v*v'),1,size(A,1));
values = sum(midMat./P)';
end % EOF


function [values] = modal_control(A)
%  FUNCTION:
%          Returns values of MODAL CONTROLLABILITY for each node in a
%          network, given the adjacency matrix for that network. Modal
%          controllability indicates the ability of that node to steer the
%          system into difficult-to-reach states, given input at that node.
% 
%  INPUT:
%          A is the structural (NOT FUNCTIONAL) network adjacency matrix,
%  	  such that the simple linear model of dynamics outlined in the
% 	  reference is an accurate estimate of brain state fluctuations.
% 	  Assumes all values in the matrix are positive, and that the
% 	  matrix is symmetric.
% 
%  OUTPUT:
%          Vector of modal controllability values for each node
% 
%  Bassett Lab, University of Pennsylvania, 2016.
%  Reference: Gu, Pasqualetti, Cieslak, Telesford, Yu, Kahn, Medaglia,
%             Vettel, Miller, Grafton & Bassett, Nature Communications
%             6:8414, 2015.

A = A./(1+svds(A,1));       %  Matrix normalization
[U, T] = schur(A,'real');   %  Schur stability
eigVals = diag(T);
N = size(A,1);
phi = zeros(N,1);
for i = 1 : N
    phi(i) = (U(i,:).^2) * (1 - eigVals.^2);
end
values = phi;
end % EOF
