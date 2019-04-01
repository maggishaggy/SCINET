clear

ACTION_path = '../../ACTION/ACTION/matlab';
addpath(genpath(ACTION_path));
addpath(genpath('./utils'));

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


%% Run ACTION 

[S_r, V, ~,  ev] = reduceGeneExpression(expression, 30, 1, 10);
[C_trace, H_trace] = runACTION(S_r, 2, 20, 8);

[selected_archs, landmark_cells, C_stacked, H_stacked, archetype_profile, backbone] = reconstructArchetypes(expression, C_trace, H_trace, -1);

%% Interpolate cells based on archetype profile
[reduced_backbone, representatives] = reduceBackbone(backbone, C_stacked, H_stacked, 3);
core_archetypes = unique(representatives);

H_final = runSimplexRegression(S_r * C_stacked(:, core_archetypes), S_r);

imputed_archetypes = V*S_r*C_stacked(:, core_archetypes);
imputed_expression = imputed_archetypes*H_final;

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


%% Construct summary networks
subsample_no = 30; % Rather small for testing. Best selected as a proportion of samples in the group
subsample_size = 10; % Rather small for test. Suggested: 30 <

UL = {'Mic'}; % Only construct network for microglia

mu = cell(numel(UL), 1);
sigma = cell(numel(UL), 1);
for j = 1:numel(UL)
    fprintf('\t\tCelltype: %s\n', UL{j});

    mask = strcmp(Labels, UL{j});
    cols = find(mask);    

    tic; mu{j} = constructNet_summary(A, net, cols, subsample_no, subsample_size, 8); toc
     
    mu{j}(mu{j} < 2 + log10(numel(mu{j}))) = 0; % Only keep very significant edges    
    mu{j} = sparse(mu{j});
end  

% Sort edges
[ii, jj, vv] = find(tril(mu{1}));
[~, perm] = sort(vv, 'descend');
edges1 = [common_genes(ii(perm)), common_genes(jj(perm)), num2cell(vv(perm))];


%% Construct network per cell
idx = find(strcmp(Labels,'Mic'));

[subs, w] = constructNet(A, net, randsample(idx, 30), 8);

w_mean = mean(w, 2);
edges2 = [common_genes(subs(:, 1)), common_genes(subs(:, 2)), num2cell(w_mean)];

[~, perm] = sort(w_mean, 'descend');
edges2 = edges2(perm, :);
