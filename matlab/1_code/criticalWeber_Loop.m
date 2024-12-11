Ohs = linspace(0,7,10);
Bo = 0.0189312135804878;

for Oh = Ohs
    try
        criticalWeber(Oh, Bo);
    catch me
        nan
    end
end