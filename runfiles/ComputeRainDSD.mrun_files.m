% Script file that calls for other functions of the model 
% % function instructions:
% RME_SetupModel(path):STEP 1 of RME model. Get the permanent info to run the model, including domain_net and layers
% RME_InitModel(path,events):STEP 2 of RME model. Get the event basic info to run the model
% RME_RunModel(path,events):STEP 3 of RME model. Simulate the raindrop microphysical process for each event
% RME_PostModel(path,events):STEP 4 of RME model. Process the result of the raindrop microphysical simulation
% Tool_GetHourlyRGPairs(path, events):Get hourly Radar Rainfall for each event before evaluation
% RME_EvaModel(path,events):STEP 5 of RME model. Evaluation of the model through radar-gauge comparison

% Modified in April 1,2021, by yaru zhang

% whole process of model using individual event
RME_SetupModel('F:\Zhang\')
RME_InitModel('test_2016-11-11.txt')
RME_RunModel('F:\Zhang\','test_2015-10-06.txt')
RME_PostModel('F:\Zhang\','test_2015-10-06.txt')
Tool_GetHourlyRGPairs('F:\Zhang\','test_2017-06-05.txt');
RME_EvaModel('F:\Zhang\','test_2015-10-06.txt');








%%
% the process using two events, display figures of experiment
RME_EvaModel('F:\Zhang\','test_comparison.txt');




