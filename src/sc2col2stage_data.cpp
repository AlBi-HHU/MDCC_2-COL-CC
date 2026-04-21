#include "sc2col_feasible.cpp"
#include <iostream>

using namespace Rcpp;


//Convert distance matrix to sorted distance list
static void
  convert_to_sorted_distances(const Rcpp::NumericMatrix& datas,std::vector<Distance>& sorted_distances,int N) {
    sorted_distances.clear();
    sorted_distances.reserve(1LL * N * (N - 1) / 2);
    int M = datas.ncol(); 
    for (int i = 0; i < N; i++)
    {
      for (int j = i + 1; j < N; j++)
      {
        double d = 0.0;
        for (int k = 0; k < M; k++)
        {
          double diff = datas(i, k) - datas(j, k);
          d += diff * diff;
        }
        d = std::sqrt(d);
        Distance e{i, j, d};
        sorted_distances.push_back(e);
      }
    }
    std::sort(sorted_distances.begin(), sorted_distances.end());
  }


// Convert distance matrix to sorted distance list,but only keep the smallest J distances
static void
  convert_to_smallest_j_distances(const Rcpp::NumericMatrix& datas,std::vector<Distance>& sorted_distances,int N,size_t J) {
    
    sorted_distances.clear();
    sorted_distances.reserve(J);
    
    // We use a max heap to keep track of the smallest J distances.
    std::priority_queue<Distance> heap;
    const double* ptr = datas.begin();
    int M = datas.ncol(); 

    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
          double d = 0.0;
          for (int k = 0; k < M; k++)
          {
            double diff = ptr[i + k * N] - ptr[j + k * N];
            d += diff * diff;
          }
          d = std::sqrt(d);
          Distance e{i, j, d};
          if (heap.size() < (size_t)J) {
            heap.push(e);
          } else if (d < heap.top().dist) {
            heap.pop();
            heap.push(e);
          }
        }
    
      }

    sorted_distances.resize(J);
    for (int i = J - 1; i >= 0; --i) {
      sorted_distances[i] = heap.top();
      heap.pop();
    }
  }




static void stage(std::vector<Distance>& sorted_distances,
                  const int N,
                  const int A,
                  const int B,
                  std::vector<std::vector<int>>& adj,
                  std::vector<std::vector<int>>& last_adj,
                  std::vector<int>& nodes,
                  std::vector<char>& node_seen,
                  std::vector<int>& color_fixated,
                  std::vector<int>& last_color_fixated,
                  std::vector<double>& dispersions_considered,
                  bool& dispersion_found,
                  FeasibleResult& res) {
  size_t index = 0;
  while (index < sorted_distances.size() && !dispersion_found) {
    //step 3
    auto dispersion = sorted_distances[index].dist;
    //step 4 & 5
    while (index < sorted_distances.size() && dispersion >= sorted_distances[index].dist) {
      add_edge_to_adj(adj, nodes, node_seen, sorted_distances, index);
    }
    //step 6
    res = sc2col_feasible(adj, nodes, color_fixated, N, A, B);
    if (res.feasible == true) {
      last_adj = adj;
      last_color_fixated = color_fixated;
      dispersions_considered.push_back(dispersion);
    }else{
      dispersion_found = true;
      dispersions_considered.push_back(dispersion);
    }
  }
}

// [[Rcpp::export]]
List solver_sc2col2stage_data(const NumericMatrix& datas,
                             const int N,
                             const IntegerVector& target_groups){
  
  int A = target_groups[0];
  int B = target_groups[1];
  bool dispersion_found = false;
  //step 1 in KAD
  std::vector<std::vector<int>> adj(N);
  std::vector<std::vector<int>> last_adj(N);
  
  std::vector<int> nodes;
  std::vector<char> node_seen(N, 0);
  
  std::vector<int> color_fixated(N, NA_INTEGER);
  std::vector<int> last_color_fixated(N, NA_INTEGER);
  
  std::vector<double> dispersions_considered;
  // Although in theory we could consider up to N*(N-1)/2 edges, in practice we often find a feasible solution with far fewer edges. 
  //If N is very large, preallocating for too many edges could lead to immediate memory issues. 
  //Therefore, we preallocate a smaller number and will dynamically expand if we exceed it.
  dispersions_considered.reserve(2*N); 
  nodes.reserve(N);

  
  FeasibleResult res;
  //step 2 in KAD
  std::vector<Distance> sorted_distances;
  convert_to_smallest_j_distances(datas, sorted_distances, N, N);

  //2Stage optimization: First use a smaller J to quickly find a feasible solution, 
  //then use the dispersion of that solution as an upper bound to filter out larger distances. 
  //Then only sort and iterate over the remaining distances. 
  //This can significantly reduce the number of distances we need to sort and iterate over, improving efficiency.
  stage(sorted_distances, N, A, B, adj, last_adj, nodes, node_seen, color_fixated, last_color_fixated, dispersions_considered, dispersion_found, res);
  if(!dispersion_found){
    //Stage 2
    adj.clear();
    adj.resize(N);
    last_adj.clear();
    last_adj.resize(N);
    nodes.clear();
    std::fill(node_seen.begin(), node_seen.end(), 0);
    std::fill(color_fixated.begin(), color_fixated.end(), NA_INTEGER);
    dispersions_considered.clear();
    sorted_distances.clear();
    convert_to_sorted_distances(datas, sorted_distances, N);
    stage(sorted_distances, N, A, B, adj, last_adj, nodes, node_seen, color_fixated, last_color_fixated, dispersions_considered, dispersion_found, res);
  }
  if (dispersions_considered.size()==1) {
    return List::create(
      Named("dispersion") = R_PosInf,
      Named("groups") = res.reason,
      Named("groups_fixated") = R_NilValue,
      Named("edges") = R_NilValue,
      Named("dispersions_considered") = NumericVector::create(R_PosInf)
    );
  } else {
    std::vector<int> final_groups = fill_isolated_nodes(last_color_fixated, A);
    Output out = convert_to_output(
      dispersions_considered.back(),
      final_groups,
      last_color_fixated,
      last_adj,
      dispersions_considered
    );
    return List::create(
      Named("dispersion") = out.dispersion,
      Named("groups") = out.groups,
      Named("groups_fixated") = out.groups_fixated,
      Named("edges") = out.edges,
      Named("dispersions_considered") = out.dispersions_considered
    );
  }
  
}








