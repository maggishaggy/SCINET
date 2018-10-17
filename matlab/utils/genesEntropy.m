function [H_norm] = genesEntropy(X)
    P = normalize(X, 'pnorm', 1, 'dim', 2);    
    I = spfun(@(x) -log2(x), P); % Information content of the relative expression of each gene       
    H = sum(P .* I, 2); % Shannon entropy of the relative expression -- overall tissue-specificity of genes. It has unit of "bits" and is between 0 -> completely selective) and log2(k) -> Uniform/HK gene    
    H_norm = H/log2(size(X, 2));
end

