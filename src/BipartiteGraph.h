//
// Created by Danial Dervovic on 30/03/2017.
//

#ifndef MBEA_BIPARTITEGRAPH_H
#define MBEA_BIPARTITEGRAPH_H

#include <vector>
#include <iostream>
#include <string>
#include "../utilities/Defines.h"
#include "../utilities/Utility.h"
#include "../utilities/Timer.h"
#include <unordered_map>
using namespace std;
class VertexDegree{
public:
    unsigned vertex;
    unsigned degree;
    inline bool operator < (const VertexDegree &other) const {
        return degree > other.degree || (degree == other.degree && vertex < other.vertex);
    }
    
    VertexDegree();
    VertexDegree(unsigned v, unsigned d);
    ~VertexDegree();
};


class BipartiteGraph {

public:

    //BipartiteGraph(const std::vector<std::vector<int>>& incidence_matrix);
    BipartiteGraph(const char *_dir,const ui alpha,const ui beta,const ui delta,const ui side,const ui strategy);
    void read_graph_binary(const char *_dir);
    void one_side_alpha_beta_pruned();
    void one_side_alpha_beta_pruned(queue<ui> removed);
    void one_side_construct_graph();
    void one_side_colorful_alpha_beta_pruned();  
    void one_side_get_num_fair_bicliques();
    void get_num_bicliques();
    void alpha_beta_pruned();
    void alpha_beta_pruned(queue<ui> removed);  
    void two_side_alpha_beta_pruned();   
    void two_side_alpha_beta_pruned(queue<ui> removed);         
    void two_side_construct_graph();
    void two_side_colorful_alpha_beta_pruned();  
    void two_side_get_num_fair_bicliques();
    void reformat_graph();
    vector<vector<ui>> combination(const vector<ui> &R_in);    
    void fairbiclique_append(const vector<ui> L_in,const vector<ui> R_in);
    void fairbiclique_append_one_side(const vector<vector<ui>> res_R,const vector<ui> L_in);
    void fairbiclique_append(const vector<ui> L_in,const vector<vector<ui>> R_in);
    vector<ui> sort_vertices(const ui strategy);
    bool check_in_biclique_one_side(const vector<ui> L_in,const vector<ui> R_in,ui strategy);
    void fairbiclique_append_two_side(const vector<vector<ui>> l_res,const vector<vector<ui>> r_res,const vector<ui> L_in,const vector<ui> R_in);
    ui get_common_neighbour(ui v,const vector<ui> set_in); 
    vector<ui> get_common_neighbour_attr(ui v,const vector<ui> set_in);
    void baseline_two_side_biclique_find();
    void baseline_one_side_biclique_find();
    void combination_(const vector<ui> vec_in,vector<ui> &now,vector<vector<ui>> &res,ui target_size,ui now_size,ui index,ui vec_in_size);    
    void baseline_two_side_biclique_find(const vector<ui> &L_in, const vector<ui> &R_in, const vector<ui> &P_in, const vector<ui> &Q_in,
                                        const ui L_in_size, const ui R_in_size,const ui P_in_size,const ui Q_in_size);
    void baseline_one_side_biclique_find(const vector<ui> &L_in, const vector<ui> &R_in, const vector<ui> &P_in, const vector<ui> &Q_in,
                                        const ui L_in_size, const ui R_in_size,const ui P_in_size,const ui Q_in_size);                                          
    vector<ui> get_common_neighbour_set(const vector<ui> set_in);    
    ui get_n();
    ui get_m();
    ui get_l_r_num();
    ui get_ui_abs(ui x,ui y);
    string vector_to_string(const vector<ui> input);
    vector<ui> get_intersection_set(const vector<ui> plus,const vector<ui> ori);
    bool is_the_same_vec(const vector<ui> plus,const vector<ui> ori);
    void baseline_one_side_fairbiclique_find(const vector<ui> &L_in, const vector<ui> &R_in, const vector<ui> &P_in, const vector<ui> &Q_in,
                                        const ui L_in_size, const vector<ui> R_in_size,const vector<ui> P_in_size,const ui Q_in_size);    
    void baseline_one_side_fairbiclique_find();                                        
    void baseline_initialize_lrpq(vector<ui> &L, vector<ui> &R,vector<ui> &P,vector<ui> &Q);
    void baseline_two_side_fairbiclique_find();
    void baseline_two_side_fairbiclique_find(const vector<ui> &L_in, const vector<ui> &R_in, const vector<ui> &P_in, const vector<ui> &Q_in,
                                        const vector<ui> L_in_size, const vector<ui> R_in_size,const vector<ui> P_in_size,const ui Q_in_size);  
    void naive_enumeration_one_side();
    void naive_enumeration_one_side_fairbiclique();
    void naive_enumeration_one_side_fairbiclique(const vector<ui> &L_in, const vector<ui> &R_in, const vector<ui> &P_in, const vector<ui> &Q_in,
                                        const ui L_in_size, const vector<ui> R_in_size,const vector<ui> P_in_size,const ui Q_in_size);
    void naive_enumeration_two_side();
    void naive_enumeration_two_side_fairbiclique();
    void naive_enumeration_two_side_fairbiclique(const vector<ui> &L_in, const vector<ui> &R_in, const vector<ui> &P_in, const vector<ui> &Q_in,
                                        const vector<ui> L_in_size, const vector<ui> R_in_size,const vector<ui> P_in_size,const ui Q_in_size);    
    bool is_a_maximal_fair_subset(const vector<ui> subset,const vector<ui> set);

    void change_mem(long long tmp,int flag);
    void change_local_mem(long long tmp,int flag);
    void clear_local_mem();
    void clear_mem();
protected:

    ui n; //#nodes of the graph
    ui m; //#edges of the graph
    ui alpha;
    ui beta;
    ui delta;
    ui side;
    ui strategy;
    ui right_start; //#V part start of the graph
    ui *pstart; //start positions of neighbors of vertices in the array "edges"
    ui *edges; //concatenation of neighbors of all vertices   
    ui *attr;
    ui *vis;       
    ui *degree;
    ui **attr_degree;
    ui *colorful_min;
    ui attr_size;
    ui *color;
    long long max_mem=0;
    long long now_mem=0;
    long long max_local_mem=0;
    long long now_local_mem=0;

    ui R_n;
    ui R_m;
    vector<ui> R_edges;
    ui *R_degree;
    ui *R_pstart;    
    vector<ui> R_nodes;      
    ui L_n;
    ui L_m;
    vector<ui> L_edges;
    ui *L_degree;
    ui *L_pstart;    
    vector<ui> L_nodes;

    ui *rank;
    
    unordered_map<string,ui> fairbiclique_map;
    ui fairbiclique_count=0;
    ui redundant_count=0;    
   
};


#endif //MBEA_BIPARTITEGRAPH_H
