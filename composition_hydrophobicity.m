function [composition_hydrophobicity_struct] = composition_hydrophobicity(ref_proteome,input_proteome)
tic
%Read in the proteome test sequences and count amino acid frequencies
proteomeStruct = fastaread(ref_proteome);
inputTest = fastaread(input_proteome);
[~,aa_string,proteome_counts] = count_aa_in_fasta(proteomeStruct)
proteome_perc = proteome_counts/sum(proteome_counts)*100;

input_counts = [];
input_names = {};

for j = 1:length(inputTest)
    aa_input = count_aa_in_string(inputTest(j).Sequence);
    input_counts(:,j) = aa_input/sum(aa_input)*100';
    input_names{j} = inputTest(j).Header;
end

%Calculate the ratio of frequencies
proteome_perc_expand = repmat(proteome_perc,1,cols(input_counts));
composition_ratio = input_counts ./ proteome_perc_expand;

%Reorder based on average values and plot
[~,ind] = sort(mean(composition_ratio,2));
composition_ratio2 = composition_ratio(ind,:);
aa_string2 = aa_string(ind);
imagesc(log2(composition_ratio2'))
set(gca, 'XTick', 1:length(aa_string2));
set(gca, 'XTickLabel', cellstr(aa_string2(:)));
set(gca, 'YTick', 1:length(inputTest));
set(gca, 'YTickLabel', input_names);
colormap(cool);colorbar
title('Amino acid enrichment log2(a2 frequency/proteome frequency)')

%Make a structure to store data
composition_hydrophobicity_struct = struct;
composition_hydrophobicity_struct.aaFrequency = composition_ratio2';
composition_hydrophobicity_struct.names = input_names';
composition_hydrophobicity_struct.aaFrequencyHeader = cellstr(num2cell(aa_string2));



% Hydrophobicity scale (Kyte-Doolittle)
hydrophobicityScale = struct('A', 1.8, 'R', -4.5, 'N', -3.5, 'D', -3.5, ...
    'C', 2.5, 'Q', -3.5, 'E', -3.5, 'G', -0.4, 'H', -3.2, 'I', 4.5, ...
    'L', 3.8, 'K', -3.9, 'M', 1.9, 'F', 2.8, 'P', -1.6, ...
    'S', -0.8, 'T', -0.7, 'W', -0.9, 'Y', -1.3, 'V', 4.2);

% Initialize hydrophobicity and charge profiles
inputCell = struct2cell(inputTest);
sequences = inputCell(2,:);
numPeptides = length(sequences);
maxLength = max(cellfun(@length, sequences)); 
hydrophobicityMatrix = zeros(numPeptides, maxLength);
chargeMatrix = zeros(numPeptides, maxLength);

% Calculate hydrophobicity and charge profiles
for i = 1:numPeptides
    seq = sequences{i};
    hydrophobicityValues = arrayfun(@(aa) hydrophobicityScale.(aa), seq, 'UniformOutput', true);
    hydrophobicityMatrix(i, 1:length(hydrophobicityValues)) = hydrophobicityValues;
end

% Calculate composite profiles (average across peptides, ignoring NaN)
compositeHydrophobicity = nanmean(hydrophobicityMatrix, 1);

% Plot individual hydrophobicity and charge profiles
figure;
hold on;
for i = 1:numPeptides
    seqLength = length(sequences{i});
    plot(1:seqLength, hydrophobicityMatrix(i, 1:seqLength), '-o', 'LineWidth', 1.5);
end
plot(1:maxLength, compositeHydrophobicity, '-k', 'LineWidth', 2.5, 'DisplayName', 'Composite Profile');
xlabel('Residue Position');
ylabel('Hydrophobicity (Kyte-Doolittle)');
title('Hydrophobicity Profiles of Peptides');
legend([input_names, {'composite'}], 'Location', 'Best');


%Make a structure to store data
composition_hydrophobicity_struct.hydrophobicityValues = hydrophobicityMatrix;
composition_hydrophobicity_struct.hydrophobicityComposite = compositeHydrophobicity;
toc
end