
Ohs = 0.015188352468351 * [2, 5, 20, 50];
Bo = 0.0189;

options = struct('prefix', 'criticalWeOh', 'save', true);
guesses = table();
matrix = table();
for Oh = Ohs
    %try
    fprintf("Finding value for Oh=%g, Bo=%g.\n", Oh, Bo);
        [We, guesses, matrix] = criticalWeber(Oh, Bo, guesses, matrix, options);
        fprintf("Found value for Oh=%g, Bo=%g. We_c=%g\n-------------------\n", Oh, Bo, We);
    %catch me
    %    disp(me);
    %end
end

notifyFinish(); % Sending me an email