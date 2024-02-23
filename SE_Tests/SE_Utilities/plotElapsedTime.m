function h = plotElapsedTime(elapsed_time, SE_type_Group)
MC = size(elapsed_time, 2);

if MC < 2
    return;
end

%The time of the first iteration is not considered because the impedence
%and admitance matrices are built in the first iteration only.
h = figure;
plot(2 : MC, elapsed_time(1, 2 : end),  2 : MC, elapsed_time(2, 2 : end));
title('elapsed time');
legend (SE_type_Group);
ylabel('s');
xlabel('iterations');

end
