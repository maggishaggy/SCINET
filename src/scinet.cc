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

namespace SCINET {
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

	field<mat> constructNet(mat A, mat net, uvec samples, int thread_no) {
		int gene_no = A.n_rows;
		int sample_no = A.n_cols;
		
		field<mat> res(2);
		uvec ind = find(trimatl(net) > 0);
		umat subs = trans(ind2sub(size(net), ind));
		
		uvec ii = subs.col(0);
		uvec jj = subs.col(1);

		res(0) = conv_to<mat>::from(subs + ones(size(subs)));
		


		samples = sort(unique(samples));	
		
		mat Af = zeros(A.n_rows, samples.n_elem);
		for(int i = 0; i < samples.n_elem; i++) {
			Af.col(i) = vec(A.col(samples(i)));
		}
		
		printf("Computing gene-specificity factor\n");	
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
		
		printf("Computing Rank-Based Inverse Normal Transformation\n");	
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
			Z.col(i) = (v_RINT(samples) + row_RINT(i)) / sqrt(2);
		}
		Z = Z.t();


		mat weights = zeros(ii.n_elem, samples.n_elem);
		
		int perc = 0, total_idx = 0;
		#pragma omp parallel for num_threads(thread_no) 
		for( int i = 0; i < samples.n_elem; i++) {			
			if(round(100*(double)total_idx/samples.n_elem) > perc) {
				perc = round(100*(double)total_idx/samples.n_elem);
				printf("%d %% done\n", perc); fflush(stdout);
				fflush(stdout);
			}	
			vec c = Z.col(i);		
			
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
			weights.col(i) = -arma::log10(v);
			
			total_idx++;
		}

		res(1) = weights;
		
		return res;
	}


	field<mat> constructNet_summary(mat A, mat net, uvec samples, int rand_sample_no, int rand_sample_size, int thread_no) {			
		thread_no = min(thread_no, rand_sample_size);
		int gene_no = A.n_rows;
		int sample_no = A.n_cols;
		
		field<mat> net_stats(2);		
		
		samples = sort(unique(samples));	
		
		mat Af = zeros(A.n_rows, samples.n_elem);
		for(int i = 0; i < samples.n_elem; i++) {
			Af.col(i) = vec(A.col(samples(i)));
		}
		
		printf("Computing gene-specificity factor\n");	
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
		
		printf("Computing Rank-Based Inverse Normal Transformation\n");	
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
			Z.col(i) = (v_RINT(samples) + row_RINT(i)) / sqrt(2);
		}
		Z = Z.t();

		uvec ind = find(net);
		vec  vv = net(ind);	
		umat subs = ind2sub(size(net), ind);
		uvec ii = trans(subs.row(0));
		uvec jj = trans(subs.row(1));


		printf("Obtaining workspace\n");
		// Instead of the whole distribution, only keep data needed to computed mean/std
		//mat subsample_weights = ones(vv.n_elem, rand_sample_no);
		vec subsample_sums = zeros(vv.n_elem, 1);
		vec subsample_sums_sq = zeros(vv.n_elem, 1);
		
		//mat rand_samples = round(stats::runif<arma::mat>(rand_sample_no, rand_sample_size, 0.0, (double)samples.n_elem-1, 0));			
		mat rand_samples = round(sampleUnif(rand_sample_no, rand_sample_size, 0.0, (double)samples.n_elem-1, 0));	
		
		int perc = 0, total_idx = 0;
		for( int k = 0; k < rand_sample_no; k++) {			
			if(round(100*(double)total_idx/rand_sample_no) > perc) {
				perc = round(100*(double)total_idx/rand_sample_no);
				printf("%d %% done\n", perc); fflush(stdout);
				fflush(stdout);
			}		
			
			//mat subsample_scores(vv.n_elem, rand_sample_size);		
			
			vec chi2_stat = zeros(vv.n_elem, 1);
			#pragma omp parallel for num_threads(thread_no) 
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
				
				chi2_stat += -2.0*arma::log(v);
			}	
						
			vec meta_pvals(vv.n_elem);
			double cum, ccum, threshold, dof = 2*rand_sample_size;
			for(int j = 0; j < vv.n_elem; j++) {			
				threshold = chi2_stat(j);
				
				cumchi( &threshold, &dof, &cum, &ccum );
				meta_pvals(j) = ccum;
			}
			
			vec pvals_corr = meta_pvals;//FDR(meta_pvals);	
			pvals_corr.transform( [](double val) { return (val < 1e-300?1e-300:val); } );
			
			vec weights = -arma::log10(pvals_corr);
			weights.replace(datum::nan, 0);  // replace each NaN with 0
			weights.replace(datum::inf, 0);  // replace each inf with 0

			
			subsample_sums += weights;
			subsample_sums_sq += square(weights);
			
			total_idx++;
		} 
		
		vec sigma;
		vec mu = subsample_sums / rand_sample_no;	
		if(rand_sample_no > 1) {
			sigma = (subsample_sums_sq - square(subsample_sums)/(rand_sample_no)) / (rand_sample_no-1);
		} 
		else {
			sigma = zeros(size(mu));
		}

		// Reconstruct networks
		int nV = net.n_cols;
		net_stats(0) = mat(sp_mat(subs, mu, nV, nV)); 
		net_stats(0) = ( net_stats(0) + trans(net_stats(0)) ) / 2; // Deal with numerical issues that violate symmetry
		
		net_stats(1) = mat(sp_mat(subs, sigma, nV, nV));
		net_stats(1) = ( net_stats(1) + trans(net_stats(1)) ) / 2; // Deal with numerical issues that violate symmetry
				
		return net_stats;
	}


	field<mat> constructNet_crossing(mat A, mat net, uvec rows1, uvec samples1, uvec rows2, uvec samples2, int rand_sample_no, int rand_sample_size, int thread_no) {
		int gene_no = A.n_rows;
		int sample_no = A.n_cols;
		
		field<mat> res(3);
		uvec ind = find(trimatl(net(rows1, rows2) > 0));
		umat subs = trans(ind2sub(size(net(rows1, rows2)), ind));
		
		uvec ii = rows1(subs.col(0));
		uvec jj = rows2(subs.col(1));
		
		umat global_subs = join_horiz(ii, jj);
		
		res(0) = conv_to<mat>::from(global_subs + ones(size(subs)));
		


		printf("Computing gene-specificity factor (for group 1)\n");	
		samples1 = sort(unique(samples1));		
		mat Af1 = zeros(A.n_rows, samples1.n_elem);
		for(int i = 0; i < samples1.n_elem; i++) {
			Af1.col(i) = vec(A.col(samples1(i)));
		}
		vec row_means1 = mean(Af1, 1);
		uvec row_perm_forward1 = stable_sort_index(row_means1);
		uvec row_perm1 = stable_sort_index(row_perm_forward1);	
		vec p1 = (row_perm1 + ones(size(row_perm1))) / (row_perm1.n_elem + 1);		
		vec row_RINT1 = zeros(size(p1));
		for (int j = 0; j < row_RINT1.n_elem; j++) {
			double norm_inv = r8_normal_01_cdf_inverse ( p1(j) );
			row_RINT1(j) = norm_inv;
		}
		


		printf("Computing gene-specificity factor (for group 1)\n");	
		samples2 = sort(unique(samples2));		
		mat Af2 = zeros(A.n_rows, samples2.n_elem);
		for(int i = 0; i < samples2.n_elem; i++) {
			Af2.col(i) = vec(A.col(samples2(i)));
		}
		vec row_means2 = mean(Af2, 1);
		uvec row_perm_forward2 = stable_sort_index(row_means2);
		uvec row_perm2 = stable_sort_index(row_perm_forward2);	
		vec p2 = (row_perm2 + ones(size(row_perm2))) / (row_perm2.n_elem + 1);		
		vec row_RINT2 = zeros(size(p2));
		for (int j = 0; j < row_RINT2.n_elem; j++) {
			double norm_inv = r8_normal_01_cdf_inverse ( p2(j) );
			row_RINT2(j) = norm_inv;
		}


				
		A = A.t();
		mat Z1(samples1.n_elem, gene_no);	
		mat Z2(samples2.n_elem, gene_no);	
		
		printf("Computing Rank-Based Inverse Normal Transformation\n");	
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
			Z1.col(i) = (v_RINT(samples1) + row_RINT1(i))/ sqrt(2);
			Z2.col(i) = (v_RINT(samples2) + row_RINT2(i))/ sqrt(2);
		}
		Z1 = Z1.t();
		Z2 = Z2.t();


		mat X = zeros(rows1.n_elem, rand_sample_no);
		umat rand_samples1 = conv_to<umat>::from(round(sampleUnif(rand_sample_no, rand_sample_size, 0.0, (double)Z1.n_cols-1, 0)));	
		for(int i = 0; i < rand_sample_no; i++) {
			X.col(i) = sum(Z1.cols(rand_samples1.row(i)), 1);
			
		}
		X /= sqrt(rand_sample_size);
		
		mat Y = zeros(rows2.n_elem, rand_sample_no);
		umat rand_samples2 = conv_to<umat>::from(round(sampleUnif(rand_sample_no, rand_sample_size, 0.0, (double)Z2.n_cols-1, 0)));	
		for(int i = 0; i < rand_sample_no; i++) {
			Y.col(i) = sum(Z2.cols(rand_samples2.row(i)), 1);
		}
		Y /= sqrt(rand_sample_size);
			
		
		
		vec subsample_sums = zeros(ii.n_elem, 1);
		vec subsample_sums_sq = zeros(ii.n_elem, 1);
		
		
		int perc = 0, total_idx = 0;
		#pragma omp parallel for num_threads(thread_no) 
		for( int k = 0; k < rand_sample_no; k++) {			
			if(round(100*(double)total_idx/rand_sample_no) > perc) {
				perc = round(100*(double)total_idx/rand_sample_no);
				printf("%d %% done\n", perc); fflush(stdout);
				fflush(stdout);
			}		
				
			vec x = X.col(k);		
			vec y = Y.col(k);
			vec min_val = arma::min(x(ii), y(jj));

			vec v = ones(min_val.n_elem);			
			uvec idx = find(0 != min_val);

			// upper tail of min stat, under indepenence assumption
			if(0 < idx.n_elem) {
				vec tail_prob = square(1.0 - normcdf(min_val(idx)));						
				v(idx) = tail_prob;
			}
			v = clamp(v, 1e-300, 1);		
			vec weights = -log10(weights);
			
			subsample_sums += weights;
			subsample_sums_sq += square(weights);
			
			total_idx++;
		} 
		
		vec sigma;
		vec mu = subsample_sums / rand_sample_no;	
		if(rand_sample_no > 1) {
			sigma = (subsample_sums_sq - square(subsample_sums)/(rand_sample_no)) / (rand_sample_no-1);
		} 
		else {
			sigma = zeros(size(mu));
		}
		
		res(1) = mu;
		res(2) = sigma;
		
		return res;
	}
}
