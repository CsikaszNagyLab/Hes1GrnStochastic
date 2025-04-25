close all

N = 200;
results = cell(1, N);

bmal_level = '100';

for iter = 1:N
    file = sprintf("./F_bmal"+bmal_level+"_200/simulation_results_cell_%d.csv", iter);
    results{iter} = readmatrix(file);
end

threshold = 200;
minLength = 300;

intervalStarts = [];
intervalLengths = [];
intervalFileSources = [];

periods = [];

onSum = 0;

for iter = 1:N
    T = results{iter}(:, 1);
    H = results{iter}(:, 2);

    on = H' >= threshold;

    if iter == 1
        firstT = T;
        firstSignal = H;
    end

    [~, peaks] = findpeaks(movmean(H,10), T, "MinPeakDistance", 30, "MinPeakProminence", 50);

    starts = find([on(1) on(2:end) & (~on(1:(end - 1)))]);
    ends = find([on(1:(end - 1)) & (~on(2:end)) on(end)]);

    %starts = T(starts);
    %ends = T(ends);
    if numel(starts) * numel(ends) == 0
        continue;
    end

    if starts(1) > ends(1)
        ends = ends(2:end);
    end
    if ends(end) < starts(end)
        starts = starts(1:(end - 1));
    end

    lengths = ends - starts + 1;
    mask = lengths > minLength;
    lengths = lengths(mask);
    starts = starts(mask);

    %xplot(T, H)

    intervalStarts = [intervalStarts starts];
    intervalLengths = [intervalLengths lengths];
    intervalFileSources = [intervalFileSources (lengths * 0 + 1) * iter];

    periods = [periods diff(peaks)'];

    onSampled = interp1(T, movmean(H,10), 0:0.1:72) >= threshold;
    onSum = onSum + onSampled;

    if numel(intervalLengths) >= 1054
        1;
    end
end

figure(1)
histogram(intervalLengths / 10 / 60, 'BinEdges',1:22)
title("Distribution of activation lengths (Bmal1 level: " + bmal_level + '%' + ")")
xlabel("Length (h)")
ylabel("Number of activations")

[sortedStarts, order] = sort(intervalStarts);
sortedLengths = intervalLengths(order) / 10;
figure(2)
hold on
for iter = 1:numel(sortedLengths)
    rectangle(Position=[sortedStarts(iter), iter, sortedLengths(iter), 1]);
end
hold off
title("Activation intervals")
xlabel("Time (min)")
set(gca,'YTick',[])

figure(3)
histogram(periods / 60)
title("Periods")
xlabel("Period length (h)")

figure(4)
plot(firstT / 60, firstSignal)
yline(threshold, "r--", "LineWidth", 3)
xlabel ("Time [h]")
ylabel("Hes1 molecule count")
title("Hes1 time series data (Bmal1 level: " + bmal_level + "%" + ")")
legend("", "Threshold")

figure(5)
plot(0:0.1:72, onSum)
title("Number of active cells (Bmal1 level: " + bmal_level + '%' + ")")
xlabel("Time [h]")
ylabel("Cell count")