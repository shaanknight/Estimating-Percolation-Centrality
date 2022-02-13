#include<bits/stdc++.h>
#include<chrono>
using namespace std;

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

// compile : g++ -O3 -static-libstdc++ approximatePercolationCentrality.cpp

const int INF = 1e9;

int N,M;
vector<vector<int> > adj;
vector<double> closeness,percolation,test_percolation;
vector<double> x;

// scale : 1 for Linear Scaling, 0 for No Scaling
// dir : 0 for Forward Brandes, 1 for Backward Brandes
void brandes(int node,vector<double> x, vector<vector<int> > &adj,double *ptr,bool scale,bool dir,double factor)
{	
	int N = (int)x.size()-1;
	queue<int> q;
	stack<int> st;
	vector<int> dist(N+1,-1);
	vector<long double> sig(N+1,0.0),delta(N+1,0.0);
	vector<vector<int> > pr(N+1);

	int u = node;
	q.push(u);
	dist[u] = 0;
	sig[u] = 1.0;

	while(!q.empty())
	{
		u = q.front();
		q.pop();
		st.push(u);

		for(auto v:adj[u])
		{
			if(dist[v] < 0)
			{
				dist[v] = dist[u]+1;
				q.push(v);
			}
			if(dist[v] == dist[u]+1)
			{
				pr[v].push_back(u);
				sig[v] = sig[u]+sig[v];
			}
		}
	}

	while(!(st.empty()))
	{
		u = st.top();
		st.pop();
		for(auto p:pr[u])
		{
			double g;
			g = sig[p]/sig[u];
			if(scale)
				g = (g*dist[p])/(double)(dist[u]);
			if(dir)
				g = g*(max(x[u]-x[node],(double)(0.0))+delta[u]);
			else
				g = g*(max(x[node]-x[u],(double)(0.0))+delta[u]);
			delta[p] = delta[p]+g;
		}
		if(u != node)
			ptr[u] += delta[u]/factor;
		pr[u].clear();
		delta[u] = 0;
		sig[u] = 0;
		dist[u] = -1;
	}
}

// Dual Brandes has Linear Scaling by default
void dual_brandes(int node,vector<double> x, vector<vector<int> > &adj,double *ptr,vector<bool> is_sampled,double factor)
{	
	int N = (int)x.size()-1;
	queue<int> q;
	stack<int> st;
	vector<int> dist(N+1,-1);
	vector<long double> sig(N+1,0.0),delta(N+1,0.0);
	vector<vector<int> > pr(N+1);

	int u = node;
	q.push(u);
	dist[u] = 0;
	sig[u] = 1.0;

	while(!q.empty())
	{
		u = q.front();
		q.pop();
		st.push(u);

		for(auto v:adj[u])
		{
			if(dist[v] < 0)
			{
				dist[v] = dist[u]+1;
				q.push(v);
			}
			if(dist[v] == dist[u]+1)
			{
				pr[v].push_back(u);
				sig[v] = sig[u]+sig[v];
			}
		}
	}

	while(!(st.empty()))
	{
		u = st.top();
		st.pop();
		for(auto p:pr[u])
		{
			double g = (sig[p]/sig[u]);
			g = (g*dist[p])/(double)(dist[u]);
			if(!is_sampled[u] || u>node)
				g *= abs(x[u]-x[node])+delta[u];
			else
				g *= delta[u]; 
			delta[p] += g;
		}
		if(u != node)
			ptr[u] += delta[u]/factor;
	}
}

void compute_closeness_centrality(int src,vector<vector<int> > &adj,double *ptr)
{
	int N = (int)x.size()-1;
	queue<int> q;
	vector<int> dist(N+1,-1);

	int u = src;
	q.push(u);
	dist[u] = 0;
	double sum = 0.0;
	int itr = 0;

	while(!q.empty())
	{
		u = q.front();
		q.pop();
		if(u != src)
			sum += (double)(dist[u]);
		itr++;
		for(auto v:adj[u])
		{
			if(dist[v] < 0)
			{
				dist[v] = dist[u]+1;
				q.push(v);
			}
		}
	} 
	if(sum)
		ptr[src] = (itr-1)/sum;
}

pair<double,vector<double> > compute_constants()
{
	vector<pair<double,int> > perc(N+1);
	vector<double> contrib(N+1,0.0);
	for(int i=1;i<=N;++i)
		perc[i] = {x[i],i};
	sort(perc.begin(),perc.end());
	long double carry = 0,sum_x = 0;
	for(int i=1;i<=N;++i)
	{
		contrib[perc[i].second] = (long double)(i-1)*perc[i].first-carry;
		carry += perc[i].first;
		sum_x += contrib[perc[i].second];
	}
	carry = 0;
	for(int i=N;i>=1;i--)
	{
		contrib[perc[i].second] += carry-(long double)(N-i)*perc[i].first;
		carry += perc[i].first;
	}
	return make_pair(sum_x,contrib);
}
// Algorithm A
void test_uniformSampling(int samples,double *ptr)
{
	auto res = compute_constants();
    double sum_x = res.first;
    vector<double> contrib = res.second;

	vector<int> permutation(N);
    for (int i = 1; i <= N; i++)
        permutation[i-1] = i;
    shuffle(permutation.begin(), permutation.end(), rng);

	for(int i=1;i<=samples;++i)
		brandes(permutation[i-1],x,adj,ptr,0,0,1.0);

	for(int i=1;i<=N;++i)
	{
		test_percolation[i] /= samples;
		test_percolation[i] /= sum_x-contrib[i];
		test_percolation[i] *= N;
	}
}
// ALgorithm B
void test_closenessBasedSampling(int samples,double *ptr)
{
	auto res = compute_constants();
    double sum_x = res.first;
    vector<double> contrib = res.second;

    fill(closeness.begin(),closeness.end(),0);

    double *ctr = &closeness[0];
	for(int i=1;i<=N;++i)
		compute_closeness_centrality(i,adj,ctr);

	double denom = 0.0;
	vector<double> prefix_sum_prob(N+1,0);
	vector<int> sampled_node;
	for(int i=1;i<=N;++i)
	{
		if(!(adj[i].empty()))
		{
			denom += x[i]*closeness[i];
			prefix_sum_prob[i] = denom;
			sampled_node.push_back(i);
		}
		else
			prefix_sum_prob[i] = prefix_sum_prob[i-1];
	}

	for(int i=1;i<=samples;++i)
	{
		double p = (0.000001*(rng()%1000000))*denom;
		int l = 0;
		int r = (int)(sampled_node.size())-1;
		int query_node = 0;
		while(l<=r)
		{
			int mid = (l+r)/2;
			query_node = sampled_node[(l+r)/2];
			if(prefix_sum_prob[query_node]>=p && prefix_sum_prob[query_node-1]<p)
				break;
			else if(prefix_sum_prob[query_node]<p)
				l = mid+1;
			else
				r = mid-1;
		}
		brandes(query_node,x,adj,ptr,0,0,(prefix_sum_prob[query_node]-prefix_sum_prob[query_node-1])/denom);
	}

	for(int i=1;i<=N;++i)
	{
		test_percolation[i] /= samples;
		test_percolation[i] /= sum_x-contrib[i];
	}
}
// Algorithm C
void test_betterApproximation(int samples,double *ptr)
{
	auto res = compute_constants();
    double sum_x = res.first;
    vector<double> contrib = res.second;

	vector<int> permutation(N);
    for (int i = 1; i <= N; i++)
        permutation[i-1] = i;
    shuffle(permutation.begin(), permutation.end(), rng);

	for(int i=1;i<=samples;++i)
	{
		int r = rng()%2;
		if(r) brandes(permutation[i-1],x,adj,ptr,1,0,0.5);
		else brandes(permutation[i-1],x,adj,ptr,1,1,0.5);
	}

	for(int i=1;i<=N;++i)
	{
		test_percolation[i] /= samples;
		test_percolation[i] /= sum_x-contrib[i];
		test_percolation[i] *= N;
	}
}
// Dual-Brandes
void test_dualBrandesNonUniform(int samples,double *ptr)
{
	auto res = compute_constants();
    double sum_x = res.first;
    vector<double> contrib = res.second;

    vector<int> set_samples;
    double denom = 0.0;
	vector<double> prefix_sum_prob(N+1,0);
	vector<int> sampled_node;
	for(int i=1;i<=N;++i)
	{
		if(!(adj[i].empty()))
		{
			denom += contrib[i];
			prefix_sum_prob[i] = denom;
			sampled_node.push_back(i);
		}
		else
			prefix_sum_prob[i] = prefix_sum_prob[i-1];
	}
	vector<bool> is_sampled(N+1,0);
	vector<double> scaling_ratio;
	for(int i=1;i<=samples;++i)
	{
		while(1)
		{
			double p = (0.000001*(rng()%1000000))*denom;
			int l = 0;
			int r = (int)(sampled_node.size())-1;
			int query_node = 0;
			while(l<=r)
			{
				int mid = (l+r)/2;
				query_node = sampled_node[(l+r)/2];
				if(prefix_sum_prob[query_node]>=p && prefix_sum_prob[query_node-1]<p)
					break;
				else if(prefix_sum_prob[query_node]<p)
					l = mid+1;
				else
					r = mid-1;
			}
			if(!is_sampled[query_node])
			{
				set_samples.push_back(query_node);
				is_sampled[query_node] = 1;
				break; 
			}
		}
	}

	for(int i=1;i<=samples;++i)
		scaling_ratio.push_back((prefix_sum_prob[set_samples[i-1]]-prefix_sum_prob[set_samples[i-1]-1])/denom);
	
	ptr = &test_percolation[0];
	for(int i=1;i<=samples;++i)
		dual_brandes(set_samples[i-1],x,adj,ptr,is_sampled,scaling_ratio[i-1]);

	for(int i=1;i<=N;++i)
	{
		test_percolation[i] /= 0.5*(samples+(1.0*(N-samples)*samples)/(1.0*N));
		test_percolation[i] /= sum_x-contrib[i];
	}
}

int calculate_diameter(vector<vector<int> > &adj)
{
	int diam = 0;
	for(int i=1;i<=10;++i)
	{
		int N = (int)adj.size()-1;
		int u = rand()%N + 1;
		int dmax = 0;

		for(int j=1;j<=20;++j)
		{
			queue<int> q;
			q.push(u);
			vector<int> dist(N+1,-1);
			dist[u] = 0;

			while(!q.empty())
			{
				u = q.front();
				q.pop();

				for(auto v:adj[u])
				{
					if(dist[v] < 0)
					{
						dist[v] = dist[u]+1;
						q.push(v);
					}
				}
			}
			if(dist[u] >= dmax)
				dmax = dist[u];
			else
				break;
		}
		diam = max(diam,dmax);
	}
	return diam;
}

int main( int argc, char **argv ) {
	ios::sync_with_stdio(false);
	cin.tie(0);
	cout.tie(0);

    string input = argv[1];
    int sampling_pcnt = atoi(argv[2]);
    string algo_choice = argv[3];

    ifstream fin(input);
	fin >> N >> M;
	int u,v;
	adj.resize(N+1);
	x.push_back(0);
	for(int i=0;i<N;++i)
	{
		double prc = 0.01*(rng()%100);
		x.push_back(prc);
	}
	for(int i=0;i<M;++i)
	{
		fin >> u >> v;
		adj[u].push_back(v);
		adj[v].push_back(u);
	}

	percolation.resize(N+1,0.0);
	test_percolation.resize(N+1,0.0);
	closeness.resize(N+1,0.0);
	// int diameter = calculate_diameter(g);
	// cerr << diameter << endl;

    auto t1 = std::chrono::high_resolution_clock::now();

    // exact PC computation
    auto res = compute_constants();
    double sum_x = res.first;
    vector<double> contrib = res.second;

	double *ptr = &percolation[0];
	for(int i=1;i<=N;++i)
		brandes(i,x,adj,ptr,0,0,1.0);

	for(int i=1;i<=N;++i)
		percolation[i] /= sum_x-contrib[i];

	auto t2 = std::chrono::high_resolution_clock::now();

	// Benchmarking Approx PC computation
    int T = 5;
    int samples = (N*sampling_pcnt)/100;
	double mse = 0.0,approx_runtime = 0.0,top_identification_ratio = 0.0;
	long long int inversions = 0;
	vector<double> top_sum_errors(101,0.0),approx_ratio(101,0.0);
	for(int tries=1;tries<=T;++tries)
	{
		auto t3 = std::chrono::high_resolution_clock::now();

		ptr = &test_percolation[0];
		if(algo_choice == "Algorithm-A") test_uniformSampling(samples,ptr);
		else if(algo_choice == "Algorithm-B") test_closenessBasedSampling(samples,ptr);
		else if(algo_choice == "Algorithm-C") test_betterApproximation(samples,ptr);
		else if(algo_choice == "Dual-Brandes") test_dualBrandesNonUniform(samples,ptr);
		else
		{
			cerr << "Non Recognisable algorithm choice, Use : " << endl;
			cerr << "Algorithm-A" << endl;
			cerr << "Algorithm-B" << endl;
			cerr << "Algorithm-C" << endl;
			cerr << "Dual-Brandes" << endl;
			return 0;
		}	    

		auto t4 = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t4 - t3 ).count();
		approx_runtime += (double)(duration);

		vector<pair<double,int> > pc_order;
		vector<int> pos_pc_order(N+1);
		for(int i=1;i<=N;++i)
			pc_order.push_back(make_pair(percolation[i],i));
		sort(pc_order.begin(),pc_order.end());
		reverse(pc_order.begin(),pc_order.end());
		for(int i=0;i<N;++i)
			pos_pc_order[pc_order[i].second] = i+1;

		vector<pair<double,int> > test_pc_order;
		vector<int> pos_test_pc_order(N+1);
		for(int i=1;i<=N;++i)
			test_pc_order.push_back(make_pair(test_percolation[i],i));
		sort(test_pc_order.begin(),test_pc_order.end());
		reverse(test_pc_order.begin(),test_pc_order.end());
		for(int i=0;i<N;++i)
			pos_test_pc_order[test_pc_order[i].second] = i+1;

		for(int i=0;i<N;++i)
		{
			int s = pc_order[i].second;
			if(percolation[s])
				mse = mse + (percolation[s]-test_percolation[s])*(percolation[s]-test_percolation[s]);
			for(int j=1;j<=100;++j)
			{
				if(i<=(N*j)/100)
				{
					if(percolation[s])
					{
						approx_ratio[j] = approx_ratio[j] + test_percolation[s]/percolation[s];
						top_sum_errors[j] = top_sum_errors[j] + abs(percolation[s]-test_percolation[s])/percolation[s];
					}
					else
						approx_ratio[j] = approx_ratio[j] + 1;
				}
			}
		}
		
		for(int i=0;i<N/100;++i)
		{
			int s = pc_order[i].second;
			if(pos_test_pc_order[s]*100<=N)
				top_identification_ratio += 100.0/N;
		}

		for(int i=0;i<N;++i)
		{
			int s = pc_order[i].second;
			for(int j=0;j<i;++j)
			{
				int t = pc_order[j].second;
				if(test_percolation[t]<test_percolation[s])
					inversions++;
			}
		}
		fill(test_percolation.begin(),test_percolation.end(),0);
	}

	cerr << "Relative Error : " << top_sum_errors[100]/(1.0*T*N) << endl;
	cerr << "Top 1% identification ratio : " << top_identification_ratio/T << endl;
	cerr << "Inversion Ratio : " << (long double)(2.0*inversions)/(long double)(1.0*T*N*(N-1)) << endl;
	cerr << "MSE : " << mse/(1.0*T*N) << endl;
	cerr << "Approx Ratio : " << approx_ratio[100]/(T*N) << endl;
	
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	cerr << "Runtime for Exact PC Computation : " << duration << " mu.s." <<endl;
	cerr << "Runtime for Approx PC Computation : " << approx_runtime/T << " mu.s." <<endl;

	return 0;
}