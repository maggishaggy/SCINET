#include <scinet.h>

double r8_normal_01_cdf_inverse ( double p );

mat sampleUnif(int l, int m, double a, double b, int seed) {
	std::default_random_engine gen (seed);	
    std::uniform_real_distribution<double> unif(a, b);
    
	
	mat R(l, m);
	for (register int j = 0; j < m; j++) {
		for(register int i = 0; i < l; i++) {
			R(i, j) = unif(gen);
		}
	}
	return R;
}

double Edgington(vec pvals) {	
	double FACT[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800.00000, 87178291200.0000, 1307674368000.00, 20922789888000.0, 355687428096000, 6.40237370572800e15, 1.21645100408832e17, 2.43290200817664e18, 5.10909421717094e19, 1.12400072777761e21, 2.58520167388850e22, 6.20448401733239e23, 1.55112100433310e25, 4.03291461126606e26, 1.08888694504184e28, 3.04888344611714e29, 8.84176199373970e30, 2.65252859812191e32, 8.22283865417792e33, 2.63130836933694e35, 8.68331761881189e36, 2.95232799039604e38, 1.03331479663861e40, 3.71993326789901e41, 1.37637530912263e43, 5.23022617466601e44, 2.03978820811974e46, 8.15915283247898e47, 3.34525266131638e49, 1.40500611775288e51, 6.04152630633738e52, 2.65827157478845e54, 1.19622220865480e56, 5.50262215981209e57, 2.58623241511168e59, 1.24139155925361e61, 6.08281864034268e62, 3.04140932017134e64};
	
	register int j;
	double combined_p = 0;
	
    int K = pvals.n_elem;
    if(K == 0)
		return 1;
		
	else if(K == 1)
		return pvals(0);
    
    double S = sum(pvals);
		
    if(S <= 1)
        combined_p = pow(S, K) / FACT[K];
        
    else {
        for(j = 0; j <= floor(S); j++) {		
			combined_p += pow(-1, j) * r8_choose(K, j) * ( pow(S-j, K) / FACT[K] );
		}
	}
	
    combined_p = std::min(1.0, combined_p);	
    
    return combined_p;
}

vec FDR(vec pvals) {	
	pvals.replace(datum::nan, 1.0);

	unsigned int V = pvals.n_elem;
	double cV = 0, prev;
	for(int k = 1; k <= V; k++)
		cV += (1.0/k);

	uvec oidx;
	vec pval_adj;
	
	oidx = sort_index(pvals);
	pvals = sort(pvals); 

	pval_adj = zeros( V, 1 );		
	prev = 1.0;
	for (unsigned int kk = V-1; kk >= 1; kk--) { 
		pval_adj(oidx(kk)) = std::min(prev, pvals(kk)*V*cV/(kk+1));
		prev = pval_adj(oidx(kk));
	}
	pval_adj(oidx(0)) = std::min(prev, pvals(0)*V*cV);
			
	return pval_adj;
}

mat identifyDiffEdges(mat mu1, mat sigma1, int n1, mat mu2, mat sigma2, int n2) {
	field<mat> diffNets(2);
	
	uvec idx = find( mu1 > 0 );
	vec delta = mu1(idx) - mu2(idx);
	vec sigma_bar = sqrt( (square(sigma1(idx)) / n1) +  (square(sigma2(idx)) / n2) );
	
	vec t_stat = delta / sigma_bar; // Welch t-test ( for "feasible" edges, in which mu1 is positive )

	mat diffNet = zeros(size(mu1));
	diffNet(idx) = t_stat;

	return diffNet;
}

field<mat> constructNet(sp_mat A, mat net, uvec samples, int rand_sample_no, int rand_sample_size, int thread_no) {			
	thread_no = min(thread_no, rand_sample_size);
	int gene_no = A.n_rows;
	int sample_no = A.n_cols;
	
	field<mat> net_stats(2);		
	
	samples = sort(unique(samples));	
	
	mat Af(size(A));
	for(int i = 0; i < samples.n_elem; i++) {
		Af.col(i) = vec(A.col(samples(i)));
	}
		
	vec row_means = mean(Af, 1);
	uvec row_perm_forward = stable_sort_index(row_means);
	uvec row_perm = stable_sort_index(row_perm_forward);	
	vec p = (row_perm + ones(size(row_perm))) / (row_perm.n_elem + 1);		
	vec row_RINT(size(p));
	for (int j = 0; j < row_RINT.n_elem; j++) {
		double norm_inv = r8_normal_01_cdf_inverse ( p(j) );
		row_RINT(j) = norm_inv;
	}
			
	A = A.t();
	mat Z(samples.n_elem, gene_no);	
	
	#pragma omp parallel for num_threads(thread_no) 
	for(int i = 0; i < gene_no; i++) {
		vec v = vec(A.col(i)); // gene i (A is already transposed)				
		
		uvec v_perm_forward = stable_sort_index(v);
		uvec v_perm = stable_sort_index(v_perm_forward);	
		vec p = (v_perm + ones(size(v_perm))) / (v_perm.n_elem + 1);		
		
		vec v_RINT(sample_no);
		for (int j = 0; j < v.n_elem; j++) {
			double norm_inv = r8_normal_01_cdf_inverse ( p(j) );
			v_RINT(j) = norm_inv;
		}
		Z.col(i) = v_RINT(samples) + row_RINT(i);
	}
	Z = Z.t();

	uvec ind = find(net);
	vec  vv = net(ind);	
	umat subs = ind2sub(size(net), ind);
	uvec ii = trans(subs.row(0));
	uvec jj = trans(subs.row(1));


	mat subsample_weights = ones(vv.n_elem, rand_sample_no);
	//mat rand_samples = round(stats::runif<arma::mat>(rand_sample_no, rand_sample_size, 0.0, (double)samples.n_elem-1, 0));			
	mat rand_samples = round(sampleUnif(rand_sample_no, rand_sample_size, 0.0, (double)samples.n_elem-1, 0));
	
	#pragma omp parallel for num_threads(thread_no) 
	for( int k = 0; k < rand_sample_no; k++) {			
		uvec idx;	
		
		mat subsample_scores(vv.n_elem, rand_sample_size);		
		for (int i = 0; i < rand_sample_size; i++) {
			vec c = Z.col(rand_samples(k, i));		
			
			vec src = c(ii);						
			vec dst = c(jj);				
			vec min_val = arma::min(src, dst);

			vec v = ones(min_val.n_elem);			
			uvec idx = find(0 != min_val);

			// upper tail of min stat, under indepenence assumption
			if(0 < idx.n_elem) {
				vec tail_prob = square(1.0 - normcdf(min_val(idx)));						
				v(idx) = tail_prob;
			}
			
			subsample_scores.col(i) = v;						
		}	
		
		mat chi2_stat = -2.0*arma::sum(arma::log(subsample_scores), 1);				
		//mat meta_pvals = ones(size(chi2_stat)) - stats::pchisq(chi2_stat, 2.0*rand_sample_size);
					
		vec meta_pvals(vv.n_elem);
		double cum, ccum, threshold, dof = 2*rand_sample_size;
		for(int j = 0; j < vv.n_elem; j++) {			
			threshold = chi2_stat(j);
			
			cumchi( &threshold, &dof, &cum, &ccum );
			meta_pvals(j) = ccum;
		}
		
		vec pvals_corr = FDR(meta_pvals);	
		pvals_corr.transform( [](double val) { return (val < 1e-300?1e-300:val); } );
		
		vec weights = -arma::log10(pvals_corr);
		
		subsample_weights.col(k) = weights;
	} 
	
	subsample_weights.replace(datum::nan, 0);  // replace each NaN with 0
	subsample_weights.replace(datum::inf, 0);  // replace each inf with 0

	int nV = net.n_cols;
	if(1 < rand_sample_no) {
		vec mu = mean(subsample_weights, 1);
		vec sigma = stddev(subsample_weights, 0, 1);
		
		net_stats(0) = mat(sp_mat(subs, mu, nV, nV));
		net_stats(1) = mat(sp_mat(subs, sigma, nV, nV));
	}
	else {
		net_stats(0) = mat(sp_mat(subs, subsample_weights, nV, nV));		
		net_stats(1) = zeros(nV, nV);
	}
	
	return net_stats;
}


