function cfdPrintCPUTime

global Domain;

deltaTime = double(tic) - Domain.currentTime;

Domain.currentTime = double(tic);
elapsedTime = Domain.currentTime - Domain.startTime;

fprintf('Iteration Time (s): %f\n',deltaTime/1e6);
fprintf('Total Elapsed Time (s): %f\n',elapsedTime/1e6);
fprintf('\n');
fprintf('\n');