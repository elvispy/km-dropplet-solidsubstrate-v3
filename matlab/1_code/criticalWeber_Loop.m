Ohs = linspace(0,.7,19);
Bo = 0.0189312135804878;

guesses = table();
matrix = table();
for Oh = Ohs
    %try
    fprintf("Finding value for Oh=%g, Bo=%g.\n", Oh, Bo);
        [We, guesses, matrix] = criticalWeber(Oh, Bo, guesses, matrix);
        fprintf("Found value for Oh=%g, Bo=%g. We_c=%g\n-------------------\n", Oh, Bo, We);
    %catch me
    %    disp(me);
    %end
end

notifyFinish(); % Sending me an email