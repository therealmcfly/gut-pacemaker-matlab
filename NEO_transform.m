function output = NEO_transform(input)
% NEO (Non-linear Energy Operator) transformation, as used by Bull et al.
% 
% INPUT: Downsampled and filtered time series signal.
% OUTPUT: Signal with NEO transofmration applied.
% Author: Larissa Marr

output = ones(1, (length(input) - 1));
for i = 2:(length(input) - 1)
    output(i) = input(i) * input(i) - input(i - 1) * input(i + 1);
end
end

