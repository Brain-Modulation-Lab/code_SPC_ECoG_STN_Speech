function M = unbal_cell2mat(C)
% C is cell
% M is matrix

% Find the length of the longest vector
max_len = max(cellfun(@length, C));

% Initialize a matrix filled with NaN
n_vectors = length(C);  % Number of vectors
M = NaN(max_len, n_vectors);  % Preallocate matrix

% Fill the matrix with vectors, padding with NaN where necessary
for i = 1:n_vectors
    current_vector = C{i};
    M(1:length(current_vector), i) = current_vector;
end