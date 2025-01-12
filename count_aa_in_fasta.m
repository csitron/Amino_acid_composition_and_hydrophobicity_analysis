function [aa_counts,AAs,aa_counts_matrix] = count_aa_in_fasta(fasta_struct)

aa_counts = struct('A', 0, 'C', 0, 'D', 0, 'E', 0, 'F', 0, ...
    'G', 0, 'H', 0, 'I', 0, 'K', 0, 'L', 0, ...
    'M', 0, 'N', 0, 'P', 0, 'Q', 0, 'R', 0, ...
    'S', 0, 'T', 0, 'V', 0, 'W', 0, 'Y', 0);

% Loop through each sequence in the FASTA file
for i = 1:length(fasta_struct)
    sequence = fasta_struct(i).Sequence;
    % Count each amino acid in the sequence
    for j = 1:length(sequence)
        amino_acid = sequence(j);
        if isfield(aa_counts, amino_acid)
            aa_counts.(amino_acid) = aa_counts.(amino_acid) + 1;
        else
            warning('Invalid character in sequence: %s', amino_acid);
        end
    end
end

AAs = 'ACDEFGHIKLMNPQRSTVWY';
aa_counts_matrix = cell2mat(struct2cell(aa_counts));

end