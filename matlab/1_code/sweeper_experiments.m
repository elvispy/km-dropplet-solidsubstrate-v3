% This sript will try to sweep simulations according to two
% rules:
% 1) Parameters set in this sweep
% 2) All empty folders will be swept. 

% STEP 1: Define which simulations are to be run. 

ddate = datestr(datetime(), 30); % 30 = ISO 8601
diary(sprintf("../2_output/Logger/%s_%s.txt", mfilename, ddate));
disp("------------");
fprintf("%s \n %s\n", string(datetime("now")), mfilename('fullpath'));


%% Setting simulation parameters
%#ok<*NOPTS>
vars = struct(...    %D = 50  %Quant = 100
    "RhoS", 1, ... % must multiply by x1000
    "SigmaS", 72.20, ... % must multiply by x100
    "R", 0.035, ... % linspace(0.02, 0.05, 5)'; % must multiply by x10 %Ang = 180
    "U", linspace(57, 47, 6)', ... %inspace(59, 39, 6)';
    "modes", 21);%tol = 5e-5

% We check how many outputs we want
numOutputs = length(fieldnames(vars));


idxs = cell(1, numOutputs);
lengths = arrayfun(@(x) length(x{1}), struct2cell(vars), 'UniformOutput', false);
idxinputs = arrayfun(@(x) 1:x{1}, lengths, 'UniformOutput', false);
[idxs{:}] = ndgrid(idxinputs{:});

fnames = fieldnames(vars);
cartesian_product = cell2mat( ...
    arrayfun(@(varidx) vars.(fnames{varidx})(idxs{varidx}, :), 1:numOutputs, ...
    'UniformOutput',false));


% Turn simulations into table
if isempty(cartesian_product) == true; cartesian_product = double.empty(0, length(fieldnames(vars))); end
simulations_cgs = array2table(cartesian_product, "VariableNames", fnames);
% Now you can manually add any simulations that you would like to run, such
% as:
%  simulations_cgs = [simulations_cgs; ...
%      {50, 100, 1, 72.20, 0, 9.78E-3, 1, 72.20, 0.001, 180, 1}];
nrow = size(simulations_cgs, 1);
simulations_cgs.folder = repmat("../2_output/", nrow, 1);



% List of files needed to run the simulation and that do not need to change
% from simul to simul.

force_sweep = false;

% STEP 3: Actually run the simulations. 

% Initial folder to go back to
root = pwd;
% A folder which MUST have all the dependencies needed (and all its
% parents, too. 
safe_folder = fullfile(root, "simulation_code");
addpath(safe_folder, '-begin');
%safe_folder = fullfile(root, "D50Quant100", "rho1000sigma7220nu98muair0", ...
%    "RhoS1000SigmaS7220", "R0350mm", "ImpDefCornerAng180U38");

final_folders = simulations_cgs.folder;

% Some common features of all the simulations to be run


%% Starting simulation
for ii = 1:height(simulations_cgs)
    %Check if etaOri exists (the center of the bath)
    cd(final_folders(ii));

    if force_sweep == true || isempty(dir("recorded*.mat")) == true
        
        numerical_parameters = struct("harmonics_qtt", simulations_cgs.modes(ii),...
            "simulation_time", inf, "version", 1);

        physical_parameters = struct("undisturbed_radius", simulations_cgs.R(ii), "initial_height", nan, ...
            "initial_velocity", simulations_cgs.U(ii), "initial_amplitudes", zeros(1, simulations_cgs.modes(ii)), ...
            "pressure_amplitudes", zeros(1, simulations_cgs.modes(ii)+1), "initial_contact_radius", 0, ...
            "rhoS", simulations_cgs.RhoS(ii), "sigmaS", simulations_cgs.SigmaS(ii));
        
        solve_motion_v2(physical_parameters, numerical_parameters);

        try
            solve_motion_v2(physical_parameters, numerical_parameters);
        catch ME
            cd(final_folders(ii))
            fprintf("Couldn't run simulation with the following parameters: \n Velocity: %g \n Modes: %g \n", ...
                simulations_cgs.U(ii), simulations_cgs.modes(ii)); 
            a = datetime('now'); a.Format = 'yyyyMMddmmss';
            save(sprintf("error_logU0=%g-%s.mat", simulations_cgs.U(ii), a), ME);
        end
    else
        fprintf("Not running simulation with the following parameters (already done): \n Velocity: %g \n Modes: %g \n", ...
                simulations_cgs.U(ii), simulations_cgs.modes(ii)); 
    end
    

end

cd(root);
delete(gcp); % Deleting current parallel workers

% Load Python3 in MACOS based on https://www.mathworks.com/matlabcentral/answers/359408-trouble-with-a-command-in-matlab-s-system
if ~ispc && system('python3 --version') ~= 0; setenv('PATH', [getenv('PATH') ':/usr/local/bin/']); end

system('python3 sending_email.py'); % Sending email to notify that's finished
diary off % turning logger off

% function final_folder = create_folder_stucture(entry)
%     base =  pwd;
%     safe_folder = fullfile(base, "D50Quant100");
% 
%     % Defining folder structure
%     physical_space = sprintf("D%gQuant%g", entry.D, entry.Quant);
%     fluid_parameters = sprintf("rho%gsigma%gnu%.0fmuair%g", entry.rho*1000, entry.sigma*100, entry.nu*10000, entry.muair);
%     sphere_parameters = sprintf("rhoS%gsigmaS%g", entry.RhoS*1000, entry.SigmaS*100);
%     radius_folder = sprintf("R%04.4gmm", entry.R*10000);
%     velocity_folder = sprintf("ImpDefCornerAng%gU%.4g", entry.Ang, entry.U);
%     modes_folder = sprintf("N=%dtol=%0.2e", entry.modes, entry.convergence_tol);
% 
%     if ~isfolder(physical_space)
%         mkdir(physical_space);
%     end
% 
%     cd(physical_space);
%     if ~isfile("zdrop.mat")
%         copyfile(fullfile(safe_folder, "ParRadDTNStops.m"), pwd);
% 
%         s = fileread(fullfile(safe_folder, "DomainMaker.m"));
%         s = regexprep(s, "D = [^\s]+;", sprintf("D = %g;", entry.D));
%         s = regexprep(s, "quant = [^\s]+;", sprintf("quant = %g;", entry.Quant));
% 
%         writeID = fopen("DomainMaker.m", 'w+');
%         fprintf(writeID, "%s", s);
%         fclose(writeID);
%         DomainMaker;
%         if ~isfile("DTN*.mat")
%             error("Integration matrix for fluid contact not found. Maually run ParRadDTNStops.m");
%         end
%     end
% 
% 
%     safe_folder = fullfile(safe_folder, "rho1000sigma7220nu98muair0");
%     if ~isfolder(fluid_parameters)
%         mkdir(fluid_parameters);
%     end
%     cd(fluid_parameters);
%     if ~isfile("muair.mat")
%         s = fileread(fullfile(safe_folder, "BathMaker.m"));
%         s = regexprep(s, "rho = [^\s]+;", sprintf("rho = %g;", entry.rho));
%         s = regexprep(s, "sigma = [^\s]+;", sprintf("sigma = %g;", entry.sigma));
%         s = regexprep(s, "nu = [^\s]+;", sprintf("nu = %.2e;", entry.nu));
%         s = regexprep(s, "muair = [^\s]+;", sprintf("muair = %g;", entry.muair));
% 
%         writeID = fopen("BathMaker.m", 'w+');
%         fprintf(writeID, "%s", s);
%         fclose(writeID);
%         BathMaker;
%     end
% 
% 
%     safe_folder = fullfile(safe_folder, "RhoS1000SigmaS7220");
%     if ~isfolder(sphere_parameters)
%         mkdir(sphere_parameters);
%        %  copyfile(fullfile(safe_folder, "DropFluidMaker.m"), pwd);
%     end
%     cd(sphere_parameters);
%     if ~isfile("sigmaS.mat")
%         s = fileread(fullfile(safe_folder, "DropFluidMaker.m"));
%         s = regexprep(s, "rhoS = [^\s]+;", sprintf("rhoS = %g;", entry.RhoS));
%         s = regexprep(s, "sigmaS = [^\s]+;", sprintf("sigmaS = %g;", entry.SigmaS));
%         writeID = fopen("DropFluidMaker.m", 'w+');
%         fprintf(writeID, "%s", s);
%         fclose(writeID);
%         DropFluidMaker;
%     end
% 
%     %safe_folder = fullfile(safe_folder, "R0350mm");
%     if ~isfolder(radius_folder)
%         mkdir(radius_folder);
%     end
%     cd(radius_folder);
%     if ~isfile("Ro.mat")
%         s = sprintf("Ro = %d; save('Ro.mat','Ro') % Ball radius in cm. This will be our unit length", entry.R);
%         writeID = fopen("RoMaker.m", "w+");
%         fprintf(writeID, "%s", s);
%         fclose(writeID);
%         RoMaker;
%     end
% 
%     %safe_folder = fullfile(safe_folder, "ImpDefCornerAng180U38");
%     if ~isfolder(velocity_folder)
%         mkdir(velocity_folder);    
%     end
%     cd(velocity_folder);
% 
%     if ~isfolder(modes_folder)
%         mkdir(modes_folder);    
%     end
%     cd(modes_folder);
% 
%     % Modify This file accordingly
%     %s = fileread(fullfile(safe_folder, "VertPolarExactSH.m"));
%     %s = regexprep(s, "U0 = 38;", sprintf("U0 = %g;", entry.U));
%     %s = regexprep(s, "N = \d+;", sprintf("N = %g;", entry.modes));
%     %writeID = fopen("VertPolarExactSH.m", 'w+');
%     %fprintf(writeID, "%s", s);
%     %fclose(writeID);
%     final_folder = pwd;
% 
%     cd(base);
% end

function parsave(fname, errormsg)
  save(fname, 'errormsg')
end

