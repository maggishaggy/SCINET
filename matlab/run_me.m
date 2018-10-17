clear

ACTION_path = '~/Dropbox/Projects/SingleCell/repositories/ACTION/ACTION/matlab';
addpath(genpath(ACTION_path));
addpath(genpath('./utils'));

subsample_no = 100;
subsample_size = 30;
%%
expression = mmread('../input/datasets/HumanBrain/Lake_Fcx_expression.mm');

gene_table = readtable('../input/datasets/HumanBrain/Lake_Fcx_gene_annotations.txt', 'Delimiter', '\t');
gene_names = gene_table.HGNC;

sample_annotations = readtable('../input/datasets/HumanBrain/Lake_Fcx_sample_annotations.txt', 'Delimiter', '\t');

filter_mask = ismember(sample_annotations.Labels, {'End', 'Per'});
expression = expression(:, ~filter_mask);
sample_annotations = sample_annotations(~filter_mask, :);

Labels = sample_annotations.Labels;
UL = unique(Labels);


ss = sum(expression, 1);
expression = median(ss)*normalize(expression, 'pnorm', 1);
expression = spfun(@(x) log(1+x), expression);

%% Run ACTION and impute missing values
[S_r, V] = reduce_GeneExpression(expression, 30);
[H, C] = run_ACTION(S_r, 2, 15, 8);

C_stacked = horzcat(C{:});
H_stacked = vertcat(H{:});

W0 = S_r*C_stacked;
[refined_C, refined_H] = refine_solution(S_r, W0);
metacell_expression = V*S_r*refined_C;
imputed_expression = metacell_expression*refined_H;

%% Read PCNet
PCNet = mmread('../input/networks/PCNet/PCNet.mm');
PCNet_genes = readtable('../input/networks/PCNet/PCNet_genes.txt');

%% Match genes
blacklist = {'UBC', 'SUMO1', 'SUMO2', 'SUMO3'};
common_genes = setdiff(intersect(gene_names, PCNet_genes.Genes), blacklist);

[~, order1] = ismember(common_genes, gene_names);
[~, order2] = ismember(common_genes, PCNet_genes.Genes);

sample_count = sum(imputed_expression(order1, :), 2);
missed_genes = find(sample_count == 0);

degs = sum(PCNet(order2, order2), 2);
isolated_vertices = find(degs == 0);

filtered_genes = union(missed_genes, isolated_vertices);
common_genes = setdiff(common_genes, common_genes(filtered_genes));


[~, order1] = ismember(common_genes, gene_names);
[~, order2] = ismember(common_genes, PCNet_genes.Genes);
A = imputed_expression(order1, :);
net = spones(PCNet(order2, order2));

% save('../results/Lake_Fcx_mapped_datasets.mat', 'A', 'net')


%% Construct networks
mu = cell(numel(UL), 1);
sigma = cell(numel(UL), 1);
for j = 1:numel(UL)
    fprintf('\t\tCelltype: %s\n', UL{j});

    mask = strcmp(Labels, UL{j});
    cols = find(mask);    

    tic; [mu{j}, sigma{j}] = mex_constructNet(A, net, cols, subsample_no, subsample_size, 8); toc
     
    mu{j}(isnan(mu{j})) = 0;
    sigma{j}(isnan(sigma{j})) = 0;
    
    mu{j} = sparse(mu{j});
    sigma{j} = sparse(sigma{j});        
end  
