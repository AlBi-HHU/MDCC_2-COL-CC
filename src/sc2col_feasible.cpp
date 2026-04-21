#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <chrono>
#include <unordered_map>
#include <queue>


using namespace Rcpp;

// Graph and connected component information structure
struct BipartiteComponent {
  std::vector<int> nodes;
  int a_count = 0;
  int b_count = 0;
};

// Feasibility result structure
struct FeasibleResult {
  bool feasible=false;
  std::string reason;
};

// Dynamic programming state structure
struct DpState {
  int component_index=-1;
  int a_sum=-1;
  int last_state_sum=-1;
  bool fliped=false;
};
//


// Distance structure
struct Distance {
  int u, v;
  double dist;
  
  bool operator<(const Distance& other) const {
    return dist < other.dist;  
  }
};

// Output result structure
struct Output {
  double dispersion;
  IntegerVector groups;
  IntegerVector groups_fixated;
  IntegerMatrix edges;
  NumericVector dispersions_considered;
};




// BFS bipartite coloring and component identification
static FeasibleResult
bfs_bipartite_and_components(const std::vector<std::vector<int>>& adj,
                             const std::vector<int>& nodes,
                             int N,
                             std::vector<int>& color_fixated,
                             std::vector<BipartiteComponent>& components) {
  FeasibleResult res;
  std::fill(color_fixated.begin(), color_fixated.end(), NA_INTEGER);
  std::queue<int> q;
  
  for (int s : nodes) {
    if (color_fixated[s] != NA_INTEGER) continue;
    
    BipartiteComponent comp;
    
    color_fixated[s] = 1;
    
    comp.nodes.push_back(s);
    comp.a_count+=1;
    q.push(s);
    while (!q.empty()) {
      int u = q.front(); 
      q.pop();
      for (int v : adj[u]) {
        if (color_fixated[v] == NA_INTEGER) {
          color_fixated[v] = color_fixated[u] == 1 ? 2 : 1;
          color_fixated[v] == 1 ? comp.a_count++ : comp.b_count++;
          comp.nodes.push_back(v);
          q.push(v);
        } else {
          if (color_fixated[v] == color_fixated[u]) {
            res.reason = "Graph is not bipartite (odd cycle detected).";
            return res;
          }
        }
      }
    }
    
    components.push_back(std::move(comp));
  }
  res.feasible = true;
  return res;
}



//Component pruning, removing the balanced nodes in each component

static FeasibleResult component_pruning(
    std::vector<BipartiteComponent>& comps,
    int& goal_A,
    int& goal_B)
{
  FeasibleResult res;
  for (auto& comp : comps) {
    int balance = std::min(comp.a_count, comp.b_count);
    comp.a_count -= balance;
    comp.b_count -= balance;
    goal_A -= balance;
    goal_B -= balance;
    if (std::min(goal_A, goal_B) < 0) {
      res.feasible = false;
      res.reason = "Infeasible: size constraint cannot be satisfied after pruning.";
      return res;
    }
  }
  res.feasible = true;
  return res;
}



//Use dynamic programming to decide whether to flip each component
static FeasibleResult dp(const std::vector<BipartiteComponent>& comps,int N,int target_a_sum, int isolated_count, std::vector<bool>& flips){
  FeasibleResult result;
  const size_t comps_count=  comps.size();
  flips.assign(comps_count, false);
  
  std::map<int, DpState> states, next_states;

  // Initial state: sum is 0, no components used, no flips
  states[0] = { -1, 0, -1, false };
  bool found = false;
  int min_needed = std::max(0, target_a_sum - isolated_count);
  
  //DP
  for (size_t i = 0; i < comps_count; ++i) {
    next_states.clear();
    //Early check if the target sum is already found
    auto it  = states.lower_bound(min_needed);
    if (it != states.end()) {
      found = true;
      break;
    }
    const int a = comps[i].a_count;
    const int b = comps[i].b_count;
    if (std::max(a,b)==0)
    {
      continue;
    }
    
    for (const auto& [sum, state] : states) {

      // not flip
      if (states.find(sum + a) == states.end()){
        if (sum + a <= target_a_sum){
          next_states[sum + a] = { (int)i, sum + a, sum, false };
        }
      }else next_states[sum + a] = states[sum + a];

      // flip
      if (states.find(sum + b) == states.end()){
        if (sum + b <= target_a_sum){
          next_states[sum + b] = { (int)i, sum + b, sum, true };
        }
      }else next_states[sum + b] = states[sum + b];
    }
    
    states.swap(next_states);

    
  }
  
  if(found){
    result.feasible = true;

    // backtrack chain construction
    std::map<int,bool> backtrack_chain;
    DpState current_state = states.lower_bound(min_needed)->second;
    while (current_state.last_state_sum != -1){ 
      backtrack_chain[current_state.component_index] = current_state.fliped;
      current_state = states[current_state.last_state_sum];
    }
    
    

    // backtracking
    for (size_t i = 0; i < comps_count; i++)
    {
      if (backtrack_chain.find(i) != backtrack_chain.end()){
        flips[i] = backtrack_chain[i];
      }else {
        auto component= comps[i];
        component.a_count == 0 ? flips[i] = false : flips[i] = true;
      }
    }
  }else{
    result.feasible = false;
    result.reason = "Dp failed to find a feasible assignment.";
  }
  return result;
}



//Construct fixated groups based on flip information and isolated nodes
static void construct_groups_fixated(
    std::vector<int>& color_fixated,
    std::vector<BipartiteComponent>& comps,
    const std::vector<bool>& flips) {
  const size_t comps_count = comps.size();
  for (size_t i = 0; i < comps_count; ++i) {
    auto& comp = comps[i];
    const bool do_flip = flips[i];
    if (do_flip)
    {
      int temp = comp.a_count;
      comp.a_count = comp.b_count;
      comp.b_count = temp;
      for (int u : comp.nodes) {
        color_fixated[u] == 1 ? color_fixated[u] = 2 : color_fixated[u] = 1;
      }
    }
  }
  
}

//Add edge to adjacency list
static void
add_edge_to_adj(std::vector<std::vector<int>>& adj,
                std::vector<int>& nodes,
                std::vector<char>& node_seen,
                std::vector<Distance>& sorted_distances,
                size_t& index) {
  auto d = sorted_distances[index];
  adj[d.u].push_back(d.v);
  adj[d.v].push_back(d.u);
  
  if (!node_seen[d.u]) {
    node_seen[d.u] = 1;
    nodes.push_back(d.u);
  }
  if (!node_seen[d.v]) {
    node_seen[d.v] = 1;
    nodes.push_back(d.v);
  }
  index++;
}

//Package and return the final result
static Output convert_to_output(double dispersion,
                                const std::vector<int>& groups,
                                const std::vector<int>& groups_fixated,
                                const std::vector<std::vector<int>>& adj,
                                const std::vector<double>& dispersions_considered){
  Output output;
  output.dispersion = dispersion;
  output.groups = IntegerVector(groups.begin(), groups.end());
  output.groups_fixated = IntegerVector(groups_fixated.begin(), groups_fixated.end());
  std::vector<std::pair<int,int>> edges_list;
  int N = adj.size();
  for (int u = 0; u < N; ++u) {
    for (int v : adj[u]) {
      if (u < v) { 
        edges_list.emplace_back(u + 1, v + 1); 
      }
    }
  }
  IntegerMatrix edges_mat(edges_list.size(), 2);
  for (size_t i = 0; i < edges_list.size(); ++i) {
    edges_mat(i, 0) = edges_list[i].first;
    edges_mat(i, 1) = edges_list[i].second;
  }
  output.edges = edges_mat;
  output.dispersions_considered = NumericVector(dispersions_considered.begin(), dispersions_considered.end());
  return output;
}
//Fill isolated nodes into groups based on remaining capacity
static std::vector<int> fill_isolated_nodes(const std::vector<int>& color_fixated,
                                            int A)  {
  
  int current_a_count = 0;
  for (size_t i = 0; i < color_fixated.size(); ++i) {
    if (color_fixated[i] == 1) {
      current_a_count++;
    }
  }
  std::vector<int> final_groups = color_fixated;
  for (size_t i = 0; i < color_fixated.size(); ++i) {
    if (color_fixated[i] == NA_INTEGER) {
      if (current_a_count < A) {
        final_groups[i] = 1;
        current_a_count++;
      } else {
        final_groups[i] = 2;
      }
    }
  }
  return final_groups;
}



// Main solver function to determine if the current graph's bipartite coloring is feasible
static FeasibleResult sc2col_feasible(std::vector<std::vector<int>>& adj,
                                      std::vector<int>& nodes,
                                      std::vector<int>& color_fixated,
                                      const int N,
                                      int A,
                                      int B) {
  
  int isolated_count = N - (int)nodes.size();
  std::vector<BipartiteComponent> comps;
  
  FeasibleResult res;
  
  res = bfs_bipartite_and_components(
    adj,nodes, N, color_fixated, comps);  
  
  if (res.feasible == false) {
    return res;
  }
  
  res = component_pruning(comps,A,B);
  
  if (res.feasible == false) {
    return res;
  }
  
  std::vector<bool> flips;

  res = dp(comps,N,A, isolated_count, flips);

  if (res.feasible == false) {
    return res;
  }
  

  construct_groups_fixated(
    color_fixated, comps, flips);
  
  return res;
}

