#include "sc2col_feasible.cpp"
#include <iostream>

using namespace Rcpp;





// Convert distance matrix to sorted distance list
static void
  convert_to_sorted_distances(const Rcpp::NumericMatrix& datas,std::vector<Distance>& sorted_distances,int N) {
    int M = datas.ncol(); 
    const double* ptr = datas.begin();
    sorted_distances.reserve(1LL * N * (N - 1) / 2);
    
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
        sorted_distances.push_back(e);
      }
    }
    std::sort(sorted_distances.begin(), sorted_distances.end());
  }




// [[Rcpp::export]]
List solver_sc2col_data(const NumericMatrix& datas,
                               const int N,
                               const IntegerVector& target_groups){
  const int A = target_groups[0];
  const int B = target_groups[1];


  bool dispersion_found = false;
  //step 1 in KAD
  std::vector<std::vector<int>> adj(N);
  std::vector<std::vector<int>> last_adj(N);

  std::vector<int> nodes;
  std::vector<char> node_seen(N, 0);

  std::vector<int> color_fixated(N, NA_INTEGER);
  std::vector<int> last_color_fixated(N, NA_INTEGER);

  std::vector<double> dispersions_considered;
  dispersions_considered.reserve(N * (N - 1) / 2); 
  nodes.reserve(N);

  FeasibleResult res;
  //step 2 in KAD
  std::vector<Distance> sorted_distances;
  size_t index = 0;
  convert_to_sorted_distances(datas,sorted_distances,N);
  
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



