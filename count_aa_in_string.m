function amino_acid_matrix = count_aa_in_string(sequence)
    % Define the order of amino acids
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY';
    
    % Initialize the matrix to store amino acid counts
    amino_acid_matrix = zeros(1, length(amino_acids));
    
    % Count each amino acid in the sequence
    for i = 1:length(sequence)
        amino_acid = sequence(i);
        idx = find(amino_acids == amino_acid, 1);
        if ~isempty(idx)
            amino_acid_matrix(idx) = amino_acid_matrix(idx) + 1;
        else
            warning('Invalid character in sequence: %s', amino_acid);
        end
    end
end