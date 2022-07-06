//
// Created by Danial Dervovic on 30/03/2017.
//

#include "BipartiteGraph.h"
#include <memory>
#include <string>
#include <fstream>
#include<iostream> 
using namespace std;
VertexDegree::VertexDegree(): vertex(-1), degree(-1){}

VertexDegree::VertexDegree(unsigned v, unsigned d): vertex(v), degree(d){}

VertexDegree::~VertexDegree(){}


bool cmpByDegAsc(VertexDegree &v1, VertexDegree &v2){
    return v1.degree < v2.degree || (v1.degree == v2.degree && v1.vertex < v2.vertex);
}

bool cmpByDegDes(VertexDegree &v1, VertexDegree &v2){
    return v1.degree > v2.degree || (v1.degree == v2.degree && v1.vertex > v2.vertex);
}


BipartiteGraph::BipartiteGraph(const char *_dir,const ui alpha_,const ui beta_,const ui delta_,const ui side_,const ui strategy_) {
	n = 0;
	m = 0;
  attr_size=2; //assume that there is only two attribute
  alpha=alpha_;
  beta=beta_;
  delta=delta_;
  side=side_;
  strategy=strategy_;
	pstart = nullptr;
	edges = nullptr;
  attr = nullptr;
  attr_degree=nullptr;
  degree=nullptr;
  colorful_min=nullptr;
  R_m=0;
  R_n=0;
  R_nodes.clear();
  R_pstart=nullptr;
  L_m=0;
  L_n=0;
  L_nodes.clear();
  fairbiclique_count=0; 
  L_pstart=nullptr;  
  read_graph_binary(_dir);
}
void BipartiteGraph::change_mem(long long tmp, int flag){
  if(flag){
    now_mem+=tmp;
    max_mem = now_mem>max_mem? now_mem: max_mem;
  }
  else now_mem-=tmp;
}

void BipartiteGraph::clear_mem(){
  now_mem=0;
  max_mem=0;
}

void BipartiteGraph::change_local_mem(long long tmp, int flag){
  if(flag){
    now_local_mem+=tmp;
    max_local_mem=now_local_mem>max_local_mem?now_local_mem:max_local_mem;
  }
  else now_local_mem-=tmp;
}

void BipartiteGraph::clear_local_mem(){
  now_local_mem=0;
}

ui BipartiteGraph::get_l_r_num(){
  ui left_count=0;
  for(ui i=0;i<right_start;i++){
    if(!vis[i])
      left_count++;
  }
  ui right_count=0;
  for(ui i=right_start;i<n;i++){
    if(!vis[i]){
      right_count++;
    }
  }
  printf("<Now total graph node number is %d,left graph node number is %d,right graph node number is %d,total edges is %d >\n",left_count+right_count,left_count,right_count,pstart[n]);
  return 0;
}

vector<ui> BipartiteGraph::get_common_neighbour_set(const vector<ui> vec_in){
    vector<ui> res;
    vector<ui> vec=vec_in;
    change_mem(vec_in.size()*sizeof(ui),1);
    for(ui i=pstart[vec_in[0]];i<pstart[vec_in[0]+1];i++){
        if(!vis[edges[i]]){
          res.emplace_back(edges[i]);
        }
    }
    for(ui i=1;i<vec_in.size();i++){
      vector<ui> res_prime;
      ui k=pstart[vec_in[i]];
      ui t=0;
      while(k<pstart[vec_in[i]+1]&&t<res.size()){
        if(vis[edges[k]]){
          k++;
        }else if(res[t]==edges[k]){
          res_prime.emplace_back(res[t]);
          t++;
          k++;
        }else if(res[t]>edges[k]){
          k++;
        }else if(edges[k]>res[t]){
          t++;
        }
      }
      res=res_prime;        
    }
    change_mem(res.size()*sizeof(ui),1);
    change_mem(res.size()*sizeof(ui),0);
    change_mem(vec_in.size()*sizeof(ui),0);
    return res;  
}

vector<ui> BipartiteGraph::sort_vertices(ui strategy){
    vector<ui> res;
    printf("strategy is %d\n",strategy);
    switch (strategy) {
        case 0: {// random vertex order
            srand(time(0));
            for(ui i=0;i<n;i++){
              if(!vis[i])
                res.emplace_back(i);
            }       
            break;
        }
        case 1: {//descending degree vertex order
            cout << "<Descend degree vertex order is applied...>" << endl;
            
            vector<VertexDegree> verdegs;
            for(ui i = 0; i < n; i++){
                if(!vis[i]){
                  VertexDegree vd(i, degree[i]);
                  verdegs.emplace_back(vd);
                }
            }
            sort(verdegs.begin(), verdegs.end(),cmpByDegDes);
            for(ui i=0;i<verdegs.size();i++){
              res.emplace_back(verdegs[i].vertex);
            }
            break;
        }
        case 2:{//ascending degree vertex order
            cout << "<Ascend degree vertex order is applied...>" << endl;            
            vector<VertexDegree> verdegs;
            for(ui i = 0; i < n; i++){
                if(!vis[i]){
                  VertexDegree vd(i, degree[i]);
                  verdegs.emplace_back(vd);
                }
            }
            sort(verdegs.begin(), verdegs.end(),cmpByDegAsc);
            for(ui i=0;i<verdegs.size();i++){
              res.emplace_back(verdegs[i].vertex);
            }
            break;
        }  
  }
  cout << "<sort vertices is over>" << endl;
  return res;
}


ui BipartiteGraph::get_n() {
  return n;
}

ui BipartiteGraph::get_m() {
  return m;
}

void BipartiteGraph::read_graph_binary(const char *dir) {
	printf("# Start reading graph, Require files \"b_degree.bin\" and \"b_adj.bin\"\n");
	FILE *f = Utility::open_file((dir + string("/b_degree.bin")).c_str(), "rb");

	ui tt;
	fread(&tt, sizeof(ui), 1, f);
	if(tt != sizeof(ui)) {
		printf("sizeof unsigned int is different: b_degree.bin(%u), machine(%lu)\n", tt, sizeof(ui));
		return ;
	}

	fread(&n, sizeof(ui), 1, f);
	fread(&m, sizeof(ui), 1, f);
	fread(&right_start, sizeof(ui), 1, f);
	printf("*\tn = %s; m = %s (undirected)\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str());

	degree = new ui[n]();
	vis = new ui[n](); 
  color= new ui[n]();
  change_mem(sizeof(ui)*n*3, 1);
  for(ui i=0;i< n;i++){
    vis[0]=0;
    color[i]=n;
  }
	fread(degree, sizeof(ui), n, f);

	fclose(f);

	long long sum = 0;
	for(ui i = 0;i < n;i ++) sum += degree[i];
	assert(sum == m);

	f = Utility::open_file((dir + string("/b_adj.bin")).c_str(), "rb");

	if(pstart != nullptr) delete[] pstart;
	pstart = new ui[n+1]();
	if(edges != nullptr) delete[] edges;
	edges = new ui[m]();
  change_mem(sizeof(ui)*(n+1), 1);
	change_mem(sizeof(ui)*(m), 1);
  pstart[0] = 0;
  int count=0;
	for(ui i = 0;i < n;i ++) {
		if(degree[i] > 0) fread(edges+pstart[i], sizeof(ui), degree[i], f);
		pstart[i+1] = pstart[i] + degree[i];
    if(degree[i] > 5){
      count++;
    }
	}
	fclose(f);
	f = Utility::open_file((dir + string("/attr.bin")).c_str(), "rb");
  attr= new ui[n]();
  change_mem(sizeof(ui)*(n), 1);
  fread(attr, sizeof(ui), n, f);
  //const ui *attr=attr_in;
  attr_degree=new ui*[n];
  for(ui i=0;i<n;i++){
    attr_degree[i]=new ui[attr_size]();
  }
  change_mem(sizeof(ui)*(n*attr_size), 1);
  for(ui i=0; i<n; i++){
    for(ui j=pstart[i];j<pstart[i+1];j++){
      ui attr_=attr[edges[j]];
      attr_degree[i][attr_]++;
    }
  }

  for(ui i=0;i<n;i++){
    if(attr[i]>attr_size||attr[i]<0){
      printf("attribute error\n");
      exit(1);
    }
  }
  get_l_r_num();
  printf("mem of graph = %lld\n", now_mem);
	printf("<here we end reading graph>\n\n\n");  
}

void BipartiteGraph::one_side_alpha_beta_pruned(){
    Timer t;
    printf("< here we start alpha beta pruned >\n");  
    queue<ui> removed;    
    for(ui i=0;i<n;i++){
      if((i>=right_start&&degree[i]<alpha)){
        vis[i]=true; 
        removed.push(i);
        continue;          
      }else if(i<right_start){
          for(ui k=0;k<attr_size;k++){
            if(attr_degree[i][k]<beta){
              vis[i]=true;
              removed.push(i);
              break;
            }
          }
        }
    }
    change_mem(sizeof(removed)+sizeof(ui)*removed.size(),1);
    one_side_alpha_beta_pruned(removed);
    printf("<time for alpha beta pruned: %s (microseconds)>\n", Utility::integer_to_string(t.elapsed()).c_str());     
    printf("< here we end alpha beta pruned >\n\n\n");      
    return ;
}

void BipartiteGraph::one_side_alpha_beta_pruned(queue<ui> removed){
    ui max_removed=removed.size();
    while(!removed.empty()){
      ui u=removed.front();
      removed.pop();
      for(ui i=pstart[u];i<pstart[u+1];i++){
        ui j=edges[i];
        if(vis[j])
          continue;
        if(u<right_start){
          if(--degree[j]<alpha){
            vis[j]=true; 
            removed.push(j);
            continue;
          }          
        }else{
          if(--attr_degree[j][attr[u]]<beta){                      
              vis[j]=true;
              removed.push(j);
              continue;
            }
        } 
      }
  }
  change_mem(sizeof(removed)+sizeof(ui)*max_removed,0);  
  return ;

}

void BipartiteGraph::one_side_colorful_alpha_beta_pruned(){
    Timer t;
    one_side_construct_graph();
    printf("< here we start colorful alpha beta pruned >\n");  
    ui* cvis=new ui[n]();
    ui* head = new ui[n]();
    ui* nxt = new ui[n]();    
    change_mem(sizeof(ui)*n*3,1);
    ui max_degree=0;
    
    for(ui i=0; i<n; i++){    
      cvis[i]=0;          
      head[i] =n;      
      color[i]=n;      
    }
    ui cut_thre=beta*attr_size-1;
    ui first_reduction=0;
    for(ui i=0; i<R_n; i++){
        ui tmp_degree=R_degree[R_nodes[i]];
        if(tmp_degree<cut_thre){
            vis[R_nodes[i]]=1;
            first_reduction++;
            continue;
        }
        nxt[R_nodes[i]]=head[tmp_degree];
    		head[tmp_degree]=R_nodes[i];
        if(tmp_degree > max_degree) 
          max_degree = tmp_degree;
    }
   	ui max_color = 0;
    for(ui ii=max_degree; ii>=1; ii--){
      for(ui jj=head[ii]; jj!=n; jj=nxt[jj]){
        for(ui j = R_pstart[jj];j < R_pstart[jj]+R_degree[jj];j ++){
          ui c = color[R_edges[j]];
          if(c != n) {
            cvis[c] = 1;
          }
        }
   			for(ui j = 0;;j ++){
         if(!cvis[j]) {
            color[jj] = j;
            if(j > max_color) 
              max_color = j;
            break;     
			    }
        }
   			for(ui j = R_pstart[jj];j < R_pstart[jj]+R_degree[jj];j ++) {
            ui c = color[R_edges[j]];
            if(c != n) 
              cvis[c] = 0;
			  }
      }
    }    
    max_color++;

    delete[] cvis;
    delete[] head;
    delete[] nxt;
    change_mem(sizeof(ui)*n*3,0);
    ui*** colorful_r =  new ui **[n];   
    ui** colorful_d = new ui*[n];   
    change_mem(sizeof(ui)*n*attr_size,1);
    change_mem(sizeof(ui)*n*attr_size*max_color,1);
    for(ui i=0; i<n; i++){
        if(vis[i])
          continue;
        colorful_d[i] = new ui[attr_size]();
        for(ui j=0; j<attr_size; j++)
            colorful_d[i][j]=0;
    }
    for(ui i=0; i<n; i++){
        if(vis[i])
          continue;
        colorful_r[i] = new ui*[attr_size];
        for(ui j=0; j<attr_size; j++){
            colorful_r[i][j] = new ui[max_color]();
            for(ui k=0; k<max_color; k++)
                colorful_r[i][j][k]=0;
        }
    }
   
    for(ui i=0; i<R_n; i++){
        if(vis[R_nodes[i]]){
          continue;    
        }
        for(ui j=R_pstart[R_nodes[i]]; j<R_pstart[R_nodes[i]]+R_degree[R_nodes[i]]; j++){
            ui neighbour=R_edges[j];
            if(vis[neighbour]==1) continue;
            if((colorful_r[R_nodes[i]][attr[neighbour]][color[neighbour]]++)==0){
                colorful_d[R_nodes[i]][attr[neighbour]]++;
            }
        }
    }
    ui deletenode=0;
    queue<ui> Q;
    ui second_reduction=0;
    for(ui i=0; i<R_n; i++){
        deletenode=0;
        ui u=R_nodes[i];
        if(vis[u]==1){
          continue;
        }
        if(colorful_d[u][attr[u]]<beta-1){
            deletenode=1; 
        }else{
          for(ui j=0; j<attr_size; j++){
            if(attr[u]==j)
              continue;
            if(colorful_d[u][j]<beta){
              deletenode=1; 
              break;
            }
          }
        }
        if(deletenode == 1){
            Q.push(u);
            vis[u]=1;
            second_reduction++;           
        }
    }
    ui maxQ=Q.size();
    change_mem(sizeof(Q)+maxQ*sizeof(ui),1);
    while(!Q.empty()){
        if(Q.size()>maxQ) maxQ=Q.size();
        ui cur=Q.front();        

        Q.pop(); 
        for(ui i=R_pstart[cur]; i<R_pstart[cur]+R_degree[cur]; i++){
            ui neighbour=R_edges[i];
            deletenode=0;            
            if(!vis[neighbour]){
                if(--colorful_r[neighbour][attr[cur]][color[cur]] <= 0){
                    colorful_d[neighbour][attr[cur]]--;
                }
                if(colorful_d[neighbour][attr[neighbour]]<beta-1){
                  deletenode=1; 
                }else{                
                  for(ui j=0; j<attr_size; j++){
                  if(j==attr[neighbour])
                    continue;
                  else if(colorful_d[neighbour][j]<beta){
                    deletenode=1; 
                    break;
                    }
                  }
                }
            if(deletenode==1){
              Q.push(neighbour);
              vis[neighbour]=1;
            }
          }
        }
    }
    change_mem(sizeof(Q)+maxQ*sizeof(ui),0);
    for(ui i=0;i<right_start;i++){
      if(!vis[i]){
        for(ui j=pstart[i];j<pstart[i+1];j++){
            ui neighbour=edges[j];            
            if(!vis[neighbour]){
              if((colorful_r[i][attr[neighbour]][color[neighbour]]++)==0){
                  colorful_d[i][attr[neighbour]]++;
              }
            }
          }
        }
    }
    queue<ui> removed; 
    for(ui i=0;i<right_start;i++){
        if(!vis[i]){
          for(ui j=0;j<attr_size;j++){
            if(colorful_d[i][j]<beta){
              removed.push(i);
              vis[i]=true;       
              break;
            }
          }
        }
    }  
    ui max_removed=removed.size();
    change_mem(sizeof(removed)+max_removed*sizeof(ui),1);
    while(!removed.empty()){
        ui u=removed.front();
        removed.pop();
        for(ui i=pstart[u];i<pstart[u+1];i++){
          ui j=edges[i];
          if(vis[j])
            continue;
          if(j>=right_start){
            degree[j]--;          
            if(degree[j]<alpha){
              vis[j]=true; 
              removed.push(j);
              continue;
            }
          }else{
            if(--colorful_r[j][attr[u]][color[u]]<=0){
                colorful_d[j][attr[u]]--;            
            }
            for(ui k=0;k<attr_size;k++){
              if(colorful_d[j][k]<beta){
                vis[j]=true; 
                removed.push(j);
                break;
              } 
            }
          }
      }
   }
   change_mem(sizeof(removed)+max_removed*sizeof(ui),0);  
   printf("< time for colorful alpha beta pruned: %s (microseconds) >\n", Utility::integer_to_string(t.elapsed()).c_str());     
   printf("< here we end colorful alpha beta pruned >\n\n\n");     
   return ;
}

void BipartiteGraph::one_side_construct_graph(){
    Timer t;
    printf("< here we start construct graph >\n");   
    R_degree= new ui[n]();
    R_pstart= new ui[n+1]();    
    
    ui *umap=new ui[n]();
    ui idx=0;
    ui *common_neigh=new ui[n]();  
    change_mem(sizeof(ui)*n*4+1,1);    
    for(ui i=right_start;i<n;i++){
      idx=0;
      if(vis[i])
        continue;
      for(ui j=pstart[i];j<pstart[i+1];j++){
        ui u=edges[j];
        if(vis[u])
          continue;
        else{
          for(ui k=pstart[u];k<pstart[u+1];k++){       
            ui v=edges[k];
            if(vis[v])
              continue;
            if(v!=i){
              umap[v]=umap[v]+1;
              if(umap[v]==1){
                common_neigh[idx++]=v;
              }
            }
          }            
        }
      }
      bool flag=true;
      for(ui j=0;j<idx;j++){
        ui k=common_neigh[j];
        if(umap[k]>=alpha){
          if(flag){
            flag=false;
            R_nodes.emplace_back(i);
            R_n++;
            R_pstart[i]=R_m;
            R_edges.emplace_back(k);
            R_m++;              
          }else{
            R_edges.emplace_back(k); 
            R_m++;
          }
        }
        umap[k]=0;
      }          
    } 
    R_nodes.emplace_back(n);
    R_pstart[n]=R_m;
    printf("R_m number is %d, m number is %d\n",R_m,m);
    change_mem(sizeof(R_nodes)+R_nodes.size()*sizeof(ui),1);
    change_mem(sizeof(R_edges)+R_edges.size()*sizeof(ui),1);
    delete umap;
    delete common_neigh;
    
    for(ui i=0;i<R_n;i++){
      R_degree[R_nodes[i]]=R_pstart[R_nodes[i+1]]-R_pstart[R_nodes[i]];
    }

    for(ui i=right_start;i<n;i++){
      if(R_degree[i]==0){
        vis[i]=1;
      }
    }
    printf("< time for construct graph: %s (microseconds) >\n", Utility::integer_to_string(t.elapsed()).c_str());      
    printf("< here we end construct graph> \n\n\n");       
    return ;
}

void BipartiteGraph::reformat_graph(){
  ui j=0;
  Timer t;
  printf("< here we start reformat graph>\n");    
  for(ui i=0;i<n;i++){
    ui tmp=pstart[i];
    pstart[i]=j;      
    for(ui k=tmp;k<pstart[i+1];k++){
      if(vis[edges[k]]){
          continue;
      }else{
          edges[j++]=edges[k];
      }
    }
    degree[i]=j-pstart[i];
  }
  pstart[n]=j;
  printf("< time for reformat graph: %s (microseconds) >\n", Utility::integer_to_string(t.elapsed()).c_str());       
  get_l_r_num();  
  printf("< here we end reformat graph> \n\n\n");       
}
void BipartiteGraph::baseline_initialize_lrpq(vector<ui> &L, vector<ui> &R,vector<ui> &P,vector<ui> &Q){
  vector<ui> res=sort_vertices(strategy);
  for(ui i=0;i<right_start;i++){
    if(vis[i])
      continue;
    else{
      L.emplace_back(i);
    }
  }

  for(ui i=0;i<res.size();i++){
    if(res[i]>=right_start){
      P.emplace_back(res[i]);
    }
  }
}


ui BipartiteGraph::get_common_neighbour(ui v,const vector<ui> set_in){
    vector<ui> input_vector=set_in;
    ui i=pstart[v],j=0,common_number=0;
    ui input_vector_size=input_vector.size();
    while(i<pstart[v+1]&&j<input_vector_size){
        if(vis[edges[i]]){
            i++;
            continue;
        }else if(edges[i]==input_vector[j]){
            i++;
            j++;
            common_number++;
        }else if(edges[i]<input_vector[j]){
            i++;
            continue;
        }else if(input_vector[j]<edges[i]){
            j++;
            continue;
        }
    }
    return common_number;
}

vector<ui> BipartiteGraph::get_common_neighbour_attr(ui v,const vector<ui> set_in){
    vector<ui> input_vector=set_in;
    ui i=pstart[v],j=0;
    vector<ui> common_number(2,0);
    ui input_vector_size=input_vector.size(); 
    while(i<pstart[v+1]&&j<input_vector_size){
        if(vis[edges[i]]){
            i++;
            continue;
        }else if(edges[i]==input_vector[j]){
            common_number[attr[edges[i]]]++;
            i++;
            j++;
        }else if(edges[i]<input_vector[j]){
            i++;
            continue;
        }else if(input_vector[j]<edges[i]){
            j++;
            continue;
        }
    }
    return common_number;
}

bool BipartiteGraph::check_in_biclique_one_side(const vector<ui> L_in,const vector<ui> R_in,ui strategy){
  vector<ui> res,generate,bound;
  if(strategy==0){
    generate=L_in;
    bound=R_in;
  }else{
    generate=R_in;
    bound=L_in;
  }
  ui generate_size=generate.size();
  ui bound_size=bound.size();
  for(ui i=pstart[generate[0]];i<pstart[generate[0]+1];i++){
    if(!vis[edges[i]]){
      res.emplace_back(edges[i]);
    }
  }
  for(ui i=1;i<generate_size;i++){
    if(res.size()<=bound_size){
      return true;
    }
    vector<ui> res_prime;
    ui k=pstart[generate[i]];
    ui t=0;
    while(k<pstart[generate[i]+1]&&t<res.size()){
      if(vis[edges[k]]){
        k++;
      }else if(res[t]==edges[k]){
        res_prime.emplace_back(res[t]);
        t++;
        k++;
      }else if(res[t]>edges[k]){
        k++;
      }else if(edges[k]>res[t]){
        t++;
      }
    }
    res=res_prime;        
  }     
  if(res.size()<=bound_size){
    return true;
  }else{
    return false;
  }  
}
void BipartiteGraph::fairbiclique_append_one_side(const vector<vector<ui>> res,const vector<ui> L_in){
  for(ui i=0;i<res.size();i++){
    vector<ui> r_res=res[i];
    if(check_in_biclique_one_side(L_in,r_res,1)){
      fairbiclique_append(L_in,r_res);
    }
  }
  
  return ;
}
void BipartiteGraph::fairbiclique_append_two_side(const vector<vector<ui>> l_res,const vector<vector<ui>> r_res,const vector<ui> L_in,const vector<ui> R_in){
  for(ui i=0;i<l_res.size();i++){
    vector<ui> l=l_res[i];

    if(check_in_biclique_one_side(l,R_in,0)){
      for(ui j=0;j<r_res.size();j++){
        vector<ui> r=r_res[j];
        if(check_in_biclique_one_side(L_in,r,1)){
          fairbiclique_append(l,r);
        }
      }
    }
  }
  
  return ;
}

void BipartiteGraph::combination_(const vector<ui> vec_in,vector<ui> &now,vector<vector<ui>> &res,ui target_size,ui now_size,ui index,ui vec_in_size){
  if(now_size==target_size){
    res.emplace_back(now);
    return ;
  }
  while(index<vec_in_size){
    if(vec_in_size-index+now_size<target_size){
      return ;
    }
    now.emplace_back(vec_in[index]);
    index++;
    now_size++;
    combination_(vec_in,now,res,target_size,now_size,index,vec_in_size);
    now.pop_back();
    now_size--;
  }
  return ;
}
string BipartiteGraph::vector_to_string(const vector<ui> input){
  string out_string="";
  vector<ui> input_number=input;
  sort(input_number.begin(),input_number.end());
  for(ui i=0;i<input_number.size();i++){
    out_string=out_string+"_"+to_string(input_number[i]); 
  }  
  return out_string;
}

void BipartiteGraph::fairbiclique_append(const vector<ui> L_in,const vector<ui> R_in){
    string out=vector_to_string(L_in);
    out=out+vector_to_string(R_in);
    ofstream fout;      
    fout.open("output.txt",ios::app);
    if(!fairbiclique_map[out]){
        fairbiclique_map[out]=1;
        fairbiclique_count++;
        fout <<out<<endl;
    }else{
      redundant_count++;
    }
    fout.close();
    return ;
}

void BipartiteGraph::fairbiclique_append(const vector<ui> L_in,const vector<vector<ui>> R_in){
    vector<ui> L=L_in;
    vector<ui> R=R_in[0];
    for(ui i=1;i<R_in.size();i++){
      R.insert(R.end(),R_in[i].begin(),R_in[i].end());
    }
    sort(L.begin(),L.end());
    sort(R.begin(),R.end());
    string out="";
    for(ui i=0;i<L.size();i++){
        out=out+"_"+to_string(L[i]); 
    }
    for(ui i=0;i<R.size();i++){
        out=out+"_"+to_string(R[i]);
    }
    if(!fairbiclique_map[out]){
        fairbiclique_map[out]=1;
        fairbiclique_count++;
        cout<<out<<endl;
    }else{
      redundant_count++;
    }
    return ;
}

vector<vector<ui>> BipartiteGraph::combination(const vector<ui> &R_in){   
    vector<vector<ui>> R_attr_vec(2);
    vector<vector<ui>> res;
    vector<ui> attr_number(2,0);
    for(ui i=0;i<R_in.size();i++){
        R_attr_vec[attr[R_in[i]]].emplace_back(R_in[i]);
        attr_number[attr[R_in[i]]]++;
    }
    if(attr_number[0]>attr_number[1]){   
      vector<ui> tmp;           
      combination_(R_attr_vec[0],tmp,res,min(attr_number[1]+delta,attr_number[0]),0,0,attr_number[0]);         
      for(ui i=0;i<res.size();i++){
        res[i].insert(res[i].end(),R_attr_vec[1].begin(),R_attr_vec[1].end());    
      }          
    }else{
      vector<ui> tmp;           
      combination_(R_attr_vec[1],tmp,res,min(attr_number[0]+delta,attr_number[1]),0,0,attr_number[1]);
      for(ui i=0;i<res.size();i++){
        res[i].insert(res[i].end(),R_attr_vec[0].begin(),R_attr_vec[0].end());    
      }               
    }
    if(res.size()==0){
      printf("error\n");
      exit(1);
    }

    return res;
}

void BipartiteGraph::baseline_one_side_fairbiclique_find(){
  vector<ui> L,R,P,Q;
  baseline_initialize_lrpq(L,R,P,Q);
  vector<ui> P_size(2,0);
  vector<ui> R_size(2,0);
  for(ui i=0;i<P.size();i++){
    P_size[attr[P[i]]]++;
  }

  baseline_one_side_fairbiclique_find(L,R,P,Q,L.size(),R_size,P_size,0);
  return ;
}
void BipartiteGraph::baseline_one_side_biclique_find(){
  vector<ui> L,R,P,Q;
  baseline_initialize_lrpq(L,R,P,Q);
  ui P_size=P.size();
  ui R_size=R.size();
  baseline_one_side_biclique_find(L,R,P,Q,L.size(),R_size,P_size,0);
  return ;
}
void BipartiteGraph::baseline_two_side_fairbiclique_find(){
  vector<ui> L,R,P,Q;
  baseline_initialize_lrpq(L,R,P,Q);
  vector<ui> P_size(2,0);
  vector<ui> R_size(2,0);
  vector<ui> L_size(2,0);  
  for(ui i=0;i<P.size();i++){
    P_size[attr[P[i]]]++;
  }
  for(ui i=0;i<L.size();i++){
    L_size[attr[L[i]]]++;
  }
  baseline_two_side_fairbiclique_find(L,R,P,Q,L_size,R_size,P_size,0);
  return ;
}

void BipartiteGraph::baseline_two_side_biclique_find(){
  vector<ui> L,R,P,Q;
  baseline_initialize_lrpq(L,R,P,Q);
  ui P_size=P.size();
  ui R_size=0;
  ui L_size=L.size();  
  baseline_two_side_biclique_find(L,R,P,Q,L_size,R_size,P_size,0);
  return ;
}

void BipartiteGraph::baseline_one_side_fairbiclique_find(const vector<ui> &L_in, const vector<ui> &R_in, const vector<ui> &P_in, const vector<ui> &Q_in,
                                        const ui L_in_size, const vector<ui> R_in_size,const vector<ui> P_in_size,const ui Q_in_size) 
{
    vector<ui> L = L_in;
    vector<ui> R = R_in;
    vector<ui> P = P_in;
    vector<ui> Q = Q_in;

    ui L_size=L_in_size;
    vector<ui> R_size=R_in_size;
    vector<ui> P_size=P_in_size;
    ui Q_size=Q_in_size;
    
    ui total_P_size=P_size[0]+P_size[1];
    
    change_mem(L.size()*sizeof(ui),1);
    change_mem(R.size()*sizeof(ui),1);
    change_mem(P.size()*sizeof(ui),1);
    change_mem(Q.size()*sizeof(ui),1);
    while(total_P_size) {
        ui x=P.back();

        vector<ui> R_prime=R;
        vector<ui> L_prime;
        vector<ui> L_overlap;
        vector<ui> R_prime_size=R_size;
        vector<ui> P_prime_size(2,0);
        ui Q_prime_size=0;      

        R_prime.emplace_back(x);
        R_prime_size[attr[x]]++;
        ui L_prime_size=0;
        vector<ui> Q_Candidate;      

        for (ui i = 0; i < L_size; i++) {
            ui u = L[i];
            bool neighbour_flag=false;
            for(ui j=pstart[u];j<pstart[u+1];j++){
                if(vis[edges[j]])
                    continue;
                if(edges[j]==x){
                    L_prime.emplace_back(u);
                    L_prime_size++;
                    neighbour_flag=true;
                    break;
                }
            }
            if(!neighbour_flag){
              L_overlap.emplace_back(u);
            }
        }
        change_mem(L_prime.size()*sizeof(ui),1);
        change_mem(R_prime.size()*sizeof(ui),1);        
        Q_Candidate.emplace_back(x);
        vector<ui> P_prime;
        vector<ui> Q_prime;
        bool is_maximal = true;
        
        for(ui j=0;j<Q_size;j++){
            ui v=Q[j];
            ui num_L_prime_neighbours=get_common_neighbour(v,L_prime);
            if (num_L_prime_neighbours == L_prime_size) {
                is_maximal = false;
                break;
            }
            else if (num_L_prime_neighbours >=alpha) {
                Q_prime.emplace_back(v);
                Q_prime_size++;
            }                        
        }
        change_mem(Q_prime.size()*sizeof(ui),1);         
        if(is_maximal){
          for (ui j = 0; j < total_P_size; j++) {
              ui v = P[j];
              if(v == x) {
                continue;
              }
              ui num_L_prime_neighbours = get_common_neighbour(v,L_prime);
              if (num_L_prime_neighbours == L_prime_size) {
                ui num_L_overlap_neighbours=get_common_neighbour(v,L_overlap);
                if(num_L_overlap_neighbours==0){
                  Q_Candidate.emplace_back(v);
                }
                R_prime.emplace_back(v);
                R_prime_size[attr[v]]++;
              }
              else if (num_L_prime_neighbours >= alpha) {
                P_prime.emplace_back(v);
                P_prime_size[attr[v]]++;                  
              }
            }
            change_mem(P_prime.size()*sizeof(ui),1);
            if(L_prime_size>=alpha){
              if(R_prime_size[0]>=beta&&R_prime_size[1]>=beta){
                ui gap=get_ui_abs(R_prime_size[0],R_prime_size[1]);
                if(gap<=delta){
                  fairbiclique_append(L_prime,R_prime);
                }else{
                  vector<vector<ui>> R_candidate=combination(R_prime);
                  for(ui i=0;i<R_candidate.size();i++){
                    vector<ui> common_L=get_common_neighbour_set(R_candidate[i]);
                    if(is_the_same_vec(common_L,L_prime)){
                      fairbiclique_append(L_prime,R_candidate[i]);
                    }
                  }
                }
              }
              if(R_prime_size[0]+P_prime_size[0]>=beta&&R_prime_size[1]+P_prime_size[1]>=beta){
                if(!P_prime.empty()){
                  baseline_one_side_fairbiclique_find(L_prime, R_prime, P_prime, Q_prime,L_prime_size,R_prime_size,P_prime_size,Q_prime_size);
                }
              }          
          }
        }
        ui P_index=0,Q_Candidate_index=0;
        vector<ui> new_P;
        while(P_index<total_P_size){
          if(P[P_index]==Q_Candidate[Q_Candidate_index]){
            P_size[attr[P[P_index]]]--;
            P_index++;
            Q_Candidate_index++;
          }else{
            new_P.emplace_back(P[P_index]);
            P_index++;
          }
        }
        P=new_P;
        total_P_size=P_size[0]+P_size[1];
        Q.insert(Q.end(),Q_Candidate.begin(),Q_Candidate.end());
        Q_size=Q.size();
        change_mem(L_prime.size()*sizeof(ui),0);
        change_mem(R_prime.size()*sizeof(ui),0);
        change_mem(P_prime.size()*sizeof(ui),0);
        change_mem(Q_prime.size()*sizeof(ui),0);
    }
    change_mem(L.size()*sizeof(ui),0);
    change_mem(R.size()*sizeof(ui),0);
    change_mem(P.size()*sizeof(ui),0);
    change_mem(Q.size()*sizeof(ui),0);
    return ;
}

void BipartiteGraph::baseline_two_side_fairbiclique_find(const vector<ui> &L_in, const vector<ui> &R_in, const vector<ui> &P_in, const vector<ui> &Q_in,
                                        const vector<ui> L_in_size, const vector<ui> R_in_size,const vector<ui> P_in_size,const ui Q_in_size) 
{
    vector<ui> L = L_in;
    vector<ui> R = R_in;
    vector<ui> P = P_in;
    vector<ui> Q = Q_in;

    vector<ui> L_size=L_in_size;
    vector<ui> R_size=R_in_size;
    vector<ui> P_size=P_in_size;
    ui Q_size=Q_in_size;
    
    ui total_P_size=P_size[0]+P_size[1];
    ui total_L_size=L_size[0]+L_size[1];    
    change_mem(L.size()*sizeof(ui),1);
    change_mem(R.size()*sizeof(ui),1);
    change_mem(P.size()*sizeof(ui),1);
    change_mem(Q.size()*sizeof(ui),1);
    while(total_P_size) {
        ui x=P.back();
        vector<ui> R_prime=R;
        vector<ui> L_prime;
        vector<ui> R_prime_size=R_size;
        vector<ui> P_prime_size(2,0);
        ui Q_prime_size=0;      
        vector<ui> Q_Candidate;   
        vector<ui> L_overlap;
        R_prime.emplace_back(x);
        R_prime_size[attr[x]]++;
        vector<ui> L_prime_size(2,0);
        for (ui i = 0; i < total_L_size; i++) {
            ui u = L[i];
            bool neighbour_flag=false;
            for(ui j=pstart[u];j<pstart[u+1];j++){
                if(vis[edges[j]])
                    continue;
                if(edges[j]==x){
                    L_prime.emplace_back(u);
                    L_prime_size[attr[u]]++;
                    neighbour_flag=true;
                    break;
                }
            }
            if(!neighbour_flag){
              L_overlap.emplace_back(u);
            }            
        }
        Q_Candidate.emplace_back(x);
        vector<ui> P_prime;
        vector<ui> Q_prime;
        bool is_maximal = true;
        ui total_L_prime_size=L_prime_size[0]+L_prime_size[1];
        for(ui j=0;j<Q_size;j++){
            ui v=Q[j];
            vector<ui> num_L_prime_neighbours=get_common_neighbour_attr(v,L_prime);
            if (num_L_prime_neighbours[0]+num_L_prime_neighbours[1] == total_L_prime_size) {
                is_maximal = false;
                break;
            }
            else if (num_L_prime_neighbours[0]>=alpha&&num_L_prime_neighbours[1] >=alpha) {
                Q_prime.emplace_back(v);
                Q_prime_size++;
            }                        
        }
        change_mem(L_prime.size()*sizeof(ui),1);
        change_mem(R_prime.size()*sizeof(ui),1);         
        if(is_maximal){
          for (ui j = 0; j < total_P_size; j++) {
              ui v = P[j];
              if(v == x) {
                continue;
              }
              vector<ui> num_L_prime_neighbours = get_common_neighbour_attr(v,L_prime);          
              if (num_L_prime_neighbours[0]+num_L_prime_neighbours[1] == total_L_prime_size) {
                ui num_L_overlap_neighbours=get_common_neighbour(v,L_overlap);
                if(num_L_overlap_neighbours==0){
                  Q_Candidate.emplace_back(v);
                }               
                R_prime.emplace_back(v);
                R_prime_size[attr[v]]++;
              }
              else if (num_L_prime_neighbours[0] >= alpha&&num_L_prime_neighbours[1]>=alpha) {                
                P_prime.emplace_back(v);
                P_prime_size[attr[v]]++;                 
              }
            }
            
            if(L_prime_size[0]>=alpha&&L_prime_size[1]>=alpha){
              if(R_prime_size[0]>=beta&&R_prime_size[1]>=beta){
                  vector<vector<ui>> R_candidate;   
                  if(get_ui_abs(R_prime_size[0],R_prime_size[1])<=delta){
                    R_candidate.emplace_back(R_prime);
                  }else{
                    R_candidate=combination(R_prime);                    
                  }                                      
                  for(ui i=0;i<R_candidate.size();i++){
                    vector<ui> common_L=get_common_neighbour_set(R_candidate[i]);                                               
                    if(is_the_same_vec(common_L,L_prime)){                   
                      vector<vector<ui>> L_candidate=combination(L_prime);
                      sort(R_candidate[i].begin(),R_candidate[i].end()); 
                      for(ui j=0;j<L_candidate.size();j++){                      
                        vector<ui> L_candidate_neighbour=get_common_neighbour_set(L_candidate[j]);                     
                        if(is_a_maximal_fair_subset(R_candidate[i],L_candidate_neighbour)){
                          fairbiclique_append(L_candidate[j],R_candidate[i]);  
                        }
                      }
                    }
                  }
              }
              change_mem(P_prime.size()*sizeof(ui),1);
              change_mem(Q_prime.size()*sizeof(ui),1);
              if(R_prime_size[0]+P_prime_size[0]>=beta&&R_prime_size[1]+P_prime_size[1]>=beta){
                if(!P_prime.empty()){
                  baseline_two_side_fairbiclique_find(L_prime, R_prime, P_prime, Q_prime,L_prime_size,R_prime_size,P_prime_size,Q_prime_size);
                }
              }
            }          
            change_mem(L_prime.size()*sizeof(ui),0);
            change_mem(R_prime.size()*sizeof(ui),0);
            change_mem(P_prime.size()*sizeof(ui),0);
            change_mem(Q_prime.size()*sizeof(ui),0);
        }
        ui P_index=0,Q_Candidate_index=0;
        vector<ui> new_P;
        while(P_index<total_P_size){
          if(P[P_index]==Q_Candidate[Q_Candidate_index]){
            P_size[attr[P[P_index]]]--;
            P_index++;
            Q_Candidate_index++;
          }else{
            new_P.emplace_back(P[P_index]);
            P_index++;
          }
        }
        P=new_P;
        total_P_size=P_size[0]+P_size[1];
        Q.insert(Q.end(),Q_Candidate.begin(),Q_Candidate.end());
        Q_size=Q.size();
    }
    change_mem(L.size()*sizeof(ui),0);
    change_mem(R.size()*sizeof(ui),0);
    change_mem(P.size()*sizeof(ui),0);
    change_mem(Q.size()*sizeof(ui),0);
    return ;
}


void BipartiteGraph::one_side_get_num_fair_bicliques(){
  one_side_alpha_beta_pruned();
  reformat_graph();
  one_side_colorful_alpha_beta_pruned();
  reformat_graph();
  Timer t;  
  baseline_one_side_fairbiclique_find();
  printf("mem of graph = %lld\n", max_mem);
  printf("< time for find fairbiclique: %s (microseconds) >\n", Utility::integer_to_string(t.elapsed()).c_str());
  printf("< find maximal fair biclique: %d, find redundant maximal fair biclique: %d >\n", fairbiclique_count,redundant_count);     
  return ;
}

void BipartiteGraph::two_side_get_num_fair_bicliques(){
  two_side_alpha_beta_pruned();
  reformat_graph();
  two_side_colorful_alpha_beta_pruned();
  reformat_graph();
  Timer t;
  baseline_two_side_fairbiclique_find();
  printf("mem of graph = %lld\n", max_mem);  
  printf("< time for find fairbiclique: %s (microseconds) >\n", Utility::integer_to_string(t.elapsed()).c_str());
  printf("< find maximal fair biclique: %d, find redundant maximal fair biclique: %d >\n", fairbiclique_count,redundant_count);     
  return ;
}

ui BipartiteGraph::get_ui_abs(ui x,ui y){
  if(x>y){
    return x-y;
  }else{
    return y-x;
  }
}

void BipartiteGraph::naive_enumeration_one_side_fairbiclique(const vector<ui> &L_in, const vector<ui> &R_in, const vector<ui> &P_in, const vector<ui> &Q_in,
                                        const ui L_in_size, const vector<ui> R_in_size,const vector<ui> P_in_size,const ui Q_in_size) 
{
    vector<ui> L = L_in;
    vector<ui> R = R_in;
    vector<ui> P = P_in;
    vector<ui> Q = Q_in;

    ui L_size=L_in_size;
    vector<ui> R_size=R_in_size;
    vector<ui> P_size=P_in_size;
    ui Q_size=Q_in_size;
    change_mem(L.size()*sizeof(ui),1);
    change_mem(R.size()*sizeof(ui),1);
    change_mem(P.size()*sizeof(ui),1);
    change_mem(Q.size()*sizeof(ui),1);
    ui total_P_size=P_size[0]+P_size[1];
    
    while(total_P_size) {
        ui x=P[0];
        vector<ui> R_prime=R;
        vector<ui> L_prime;
        vector<ui> P_prime;
        vector<ui> Q_prime;  
        vector<ui> R_prime_size=R_size;
        ui L_prime_size=0;
        vector<ui> P_prime_size(2,0);      
        ui Q_prime_size=0;

        R_prime.emplace_back(x);
        R_prime_size[attr[x]]++;        

        bool is_maximal = true;      
        for (ui i = 0; i < L_size; i++) {
            ui u = L[i];
            for(ui j=pstart[u];j<pstart[u+1];j++){
                if(vis[edges[j]])
                    continue;
                if(edges[j]==x){
                    L_prime.emplace_back(u);
                    L_prime_size++;
                    break;
                }
            }
        }
        change_mem(L_prime.size()*sizeof(ui),1);
        change_mem(R_prime.size()*sizeof(ui),1);        
        vector<ui> Q_attr_size(2,0);
        for(ui i=0;i<Q_size;i++){
          ui v=Q[i];
          ui num_L_prime_neighbours=get_common_neighbour(v,L_prime);
          if (num_L_prime_neighbours == L_prime_size) {
            Q_attr_size[attr[v]]++;
          } 
          if(num_L_prime_neighbours>=alpha){
            Q_prime.push_back(Q[i]);
            Q_prime_size++;
          }
        }
        if(L_prime_size<alpha){
          is_maximal=false;          
        }
        if(Q_attr_size[0]>=1&&Q_attr_size[1]>=1){
          is_maximal=false;
        }  
        if(is_maximal){
          vector<ui> R_prime_plus;
          vector<ui> R_prime_plus_size(2,0);
          for (ui i = 0; i < total_P_size; i++) {
            ui v = P[i];
            if(v == x) {
              continue;
            }
            ui num_L_prime_neighbours = get_common_neighbour(v,L_prime);
            if (num_L_prime_neighbours == L_prime_size) {
              R_prime_plus.emplace_back(v);
              R_prime_plus_size[attr[v]]++;
              P_prime.emplace_back(v);
              P_prime_size[attr[v]]++;
            }
            else if (num_L_prime_neighbours >=alpha) {
              P_prime.emplace_back(v);
              P_prime_size[attr[v]]++;          
            }
          }
          if(R_prime_plus_size[0]+R_prime_plus_size[1]==P_prime_size[0]+P_prime_size[1]){
            ui gap=get_ui_abs(R_prime_size[0]+R_prime_plus_size[0],R_prime_size[1]+R_prime_plus_size[1]);
            if(gap<=delta){
              R_prime_size[0]=R_prime_size[0]+R_prime_plus_size[0];
              R_prime_size[1]=R_prime_size[1]+R_prime_plus_size[1];
              R_prime.insert(R_prime.end(),R_prime_plus.begin(),R_prime_plus.end());
              R_prime_plus_size[0]=0;
              R_prime_plus_size[1]=0;
              vector <ui>().swap(P_prime);
              P_prime_size[0]=0;
              P_prime_size[1]=0;
            }
          }
                    
          ui gap=get_ui_abs(R_prime_size[0],R_prime_size[1]);

          if(gap>delta){
            is_maximal=false;
          }      
          if((R_prime_size[0]<beta)||(R_prime_size[1]<beta)){
            is_maximal=false;            
          }         
          if(is_maximal){    
                if(R_prime_plus_size[0]+Q_attr_size[0]>=1&&R_prime_plus_size[1]+Q_attr_size[1]>=1){
                  is_maximal=false;                
                }else if(R_prime_plus_size[0]+Q_attr_size[0]==0&&R_prime_plus_size[1]+Q_attr_size[1]>=1){
                  gap=get_ui_abs(R_prime_size[0],R_prime_size[1]+1);
                  if(gap<=delta){
                    is_maximal=false;
                  }
                }else if(R_prime_plus_size[1]+Q_attr_size[1]==0&&R_prime_plus_size[0]+Q_attr_size[0]>=1){
                  gap=get_ui_abs(R_prime_size[0]+1,R_prime_size[1]);                
                  if(gap<=delta){
                    is_maximal=false;
                  }                
                }                                 
                if(is_maximal){                       
                  fairbiclique_append(L_prime,R_prime);
                }                
          }
          
        bool continue_maximal_flag=true;
        if((R_prime_size[0]+P_prime_size[0]<beta)||(R_prime_size[1]+P_prime_size[1]<beta)){
          continue_maximal_flag=false;
        }
        change_mem(P_prime.size()*sizeof(ui),1);
        change_mem(Q_prime.size()*sizeof(ui),1);
        if ((P_prime_size[0]+P_prime_size[1])&&(continue_maximal_flag)) {
          naive_enumeration_one_side_fairbiclique(L_prime, R_prime, P_prime, Q_prime,L_prime_size,R_prime_size,P_prime_size,Q_prime_size);
        }

      }
      P_size[attr[x]]--;      
      P.erase(P.begin());
      total_P_size--;
      Q.emplace_back(x);
      Q_size++;    
      change_mem(L_prime.size()*sizeof(ui),0);
      change_mem(R_prime.size()*sizeof(ui),0);
      change_mem(P_prime.size()*sizeof(ui),0);
      change_mem(Q_prime.size()*sizeof(ui),0);      
    }
    change_mem(L.size()*sizeof(ui),0);
    change_mem(R.size()*sizeof(ui),0);
    change_mem(P.size()*sizeof(ui),0);
    change_mem(Q.size()*sizeof(ui),0);
    return ;
}

void BipartiteGraph::naive_enumeration_one_side_fairbiclique(){
  vector<ui> L,R,P,Q;
  baseline_initialize_lrpq(L,R,P,Q);
  vector<ui> R_size(2,0);
  vector<ui> P_size(2,0);
  for(ui i=0;i<P.size();i++){
    P_size[attr[P[i]]]++;
  }
  naive_enumeration_one_side_fairbiclique(L,R,P,Q,L.size(),R_size,P_size,0);
}

void BipartiteGraph::naive_enumeration_one_side(){
  one_side_alpha_beta_pruned();
  reformat_graph();
  one_side_colorful_alpha_beta_pruned();
  reformat_graph();
  Timer t;  
  naive_enumeration_one_side_fairbiclique();
  printf("mem of graph = %lld\n", max_mem);    
  printf("< time for find fairbiclique: %s (microseconds) >\n", Utility::integer_to_string(t.elapsed()).c_str());
  printf("< find maximal fair biclique: %d, find redundant maximal fair biclique: %d >\n", fairbiclique_count,redundant_count);
  return ;
}

void BipartiteGraph::two_side_colorful_alpha_beta_pruned(){
    Timer t;
    two_side_construct_graph();    
    printf("< here we start colorful alpha beta pruned >\n");  
    ui* cvis=new ui[n]();
    ui* head = new ui[n]();
    ui* nxt = new ui[n]();    
    colorful_min=new ui[n]();
    change_mem(sizeof(ui)*n*4,1);
    ui totol_max_color=0;
    ui max_degree=0;
    
    for(ui i=0; i<n; i++){    
      cvis[i]=0;          
      head[i]=n;
      color[i]=n;            
    }
    ui cut_thre=beta*attr_size-1;
    for(ui i=0; i<R_n; i++){
        ui tmp_degree=R_degree[R_nodes[i]];
        if(tmp_degree<cut_thre){
            vis[R_nodes[i]]=1;
            continue;
        }
        nxt[R_nodes[i]]=head[tmp_degree];
    		head[tmp_degree]=R_nodes[i];
        if(tmp_degree > max_degree) 
          max_degree = tmp_degree;
    }

   	ui max_color = 0;
    for(ui ii=max_degree; ii>=1; ii--){
      for(ui jj=head[ii]; jj!=n; jj=nxt[jj]){
        for(ui j = R_pstart[jj];j < R_pstart[jj]+R_degree[jj];j ++){
          ui c = color[R_edges[j]];
          if(c != n) {
            cvis[c] = 1;
          }
        }
   			for(ui j = 0;;j ++){
         if(!cvis[j]) {
            color[jj] = j;
            if(j > max_color) 
              max_color = j;
            break;     
			    }
        }
   			for(ui j = R_pstart[jj];j < R_pstart[jj]+R_degree[jj];j ++) {
            ui c = color[R_edges[j]];
            if(c != n) 
              cvis[c] = 0;
			  }
      }
    }
    max_color++;
    totol_max_color=max(max_color,totol_max_color);
    ui*** colorful_r =  new ui **[n];   
    ui** colorful_d = new ui*[n];   
    change_mem(sizeof(ui)*n*attr_size,1);
    change_mem(sizeof(ui)*n*attr_size*max_color,1);
    for(ui i=0; i<n; i++){
        if(vis[i])
          continue;
        colorful_d[i] = new ui[attr_size]();
        for(ui j=0; j<attr_size; j++)
            colorful_d[i][j]=0;
    }
    for(ui i=right_start; i<n; i++){
        if(vis[i])
          continue;
        colorful_r[i] = new ui*[attr_size];
        for(ui j=0; j<attr_size; j++){
            colorful_r[i][j] = new ui[max_color]();
            for(ui k=0; k<max_color; k++)
                colorful_r[i][j][k]=0;
        }
    }
   
    for(ui i=0; i<R_n; i++){
        if(vis[R_nodes[i]]){
          continue;    
        }
        for(ui j=R_pstart[R_nodes[i]]; j<R_pstart[R_nodes[i]]+R_degree[R_nodes[i]]; j++){
            ui neighbour=R_edges[j];

            if(vis[neighbour]==1) continue;
            if((colorful_r[R_nodes[i]][attr[neighbour]][color[neighbour]]++)==0){
                colorful_d[R_nodes[i]][attr[neighbour]]++;
            }
        }
    }

    //Here we start construct left 2-hop graph
    max_degree=0;
    for(ui i=0; i<n; i++){    
      cvis[i]=0;          
      head[i] =n;      
    }    
    cut_thre=alpha*attr_size-1;
    for(ui i=0; i<L_n; i++){
        ui tmp_degree=L_degree[L_nodes[i]];
        if(tmp_degree<cut_thre){
            vis[L_nodes[i]]=1;
            continue;
        }
        nxt[L_nodes[i]]=head[tmp_degree];
    		head[tmp_degree]=L_nodes[i];
        if(tmp_degree > max_degree) 
          max_degree = tmp_degree;
    }

   	max_color = 0;
    for(ui ii=max_degree; ii>=1; ii--){
      for(ui jj=head[ii]; jj!=n; jj=nxt[jj]){
        for(ui j = L_pstart[jj];j < L_pstart[jj]+L_degree[jj];j ++){
          ui c = color[L_edges[j]];
          if(c != n) {
            cvis[c] = 1;
          }
        }
   			for(ui j = 0;;j ++){
         if(!cvis[j]) {
            color[jj] = j;
            if(j > max_color) 
              max_color = j;
            break;     
			    }
        }
   			for(ui j = L_pstart[jj];j < L_pstart[jj]+L_degree[jj];j ++) {
            ui c = color[L_edges[j]];
            if(c != n) 
              cvis[c] = 0;
			  }
      }
    }
    max_color++;
    totol_max_color=max(max_color,totol_max_color);    
    delete[] cvis;
    delete[] head;
    delete[] nxt;
    change_mem(sizeof(ui)*n*3,0);    
    for(ui i=0; i<right_start; i++){
        if(vis[i])
          continue;
        colorful_r[i] = new ui*[attr_size];
        for(ui j=0; j<attr_size; j++){
            colorful_r[i][j] = new ui[max_color]();
            for(ui k=0; k<max_color; k++)
                colorful_r[i][j][k]=0;
        }
    }

    for(ui i=0; i<L_n; i++){
        if(vis[L_nodes[i]]){
          continue;    
        }
        for(ui j=L_pstart[L_nodes[i]]; j<L_pstart[L_nodes[i]]+L_degree[L_nodes[i]]; j++){
            ui neighbour=L_edges[j];

            if(vis[neighbour]==1) continue;
            if((colorful_r[L_nodes[i]][attr[neighbour]][color[neighbour]]++)==0){
                colorful_d[L_nodes[i]][attr[neighbour]]++;
            }
        }
    }    



    ui deletenode=0;
    queue<ui> Q;
    for(ui i=0; i<R_n; i++){
        deletenode=0;
        ui u=R_nodes[i];
        if(vis[u]==1){
          continue;
        }
        if(colorful_d[u][attr[u]]<beta-1){
            deletenode=1; 
        }else{
          for(ui j=0; j<attr_size; j++){
            if(attr[u]==j)
              continue;
            if(colorful_d[u][j]<beta){
              deletenode=1; 
              break;
            }
          }
        }
        if(deletenode == 1){
            Q.push(u);
            vis[u]=1;
        }
    }

    for(ui i=0; i<L_n; i++){
        deletenode=0;
        ui u=L_nodes[i];
        if(vis[u]==1){
          continue;
        }
        if(colorful_d[u][attr[u]]<alpha-1){
            deletenode=1; 
        }else{
          for(ui j=0; j<attr_size; j++){
            if(attr[u]==j)
              continue;
            if(colorful_d[u][j]<alpha){
              deletenode=1; 
              break;
            }
          }
        }
        if(deletenode == 1){
            Q.push(u);
            vis[u]=1;
        }
    }

    ui maxQ=Q.size();
    while(!Q.empty()){
        if(Q.size()>maxQ) maxQ=Q.size();
        ui cur=Q.front();        
        Q.pop(); 
        if(cur>=right_start){
          for(ui i=R_pstart[cur]; i<R_pstart[cur]+R_degree[cur]; i++){
              ui neighbour=R_edges[i];
              deletenode=0;            
              if(!vis[neighbour]){                                  
                  if(--colorful_r[neighbour][attr[cur]][color[cur]] <= 0){
                      colorful_d[neighbour][attr[cur]]--;
                  }
                  if(colorful_d[neighbour][attr[neighbour]]<beta-1){
                    deletenode=1; 
                  }else{                
                    for(ui j=0; j<attr_size; j++){
                      if(j==attr[neighbour])
                        continue;
                    else if(colorful_d[neighbour][j]<beta){
                      deletenode=1; 
                      break;
                      }
                    }
                  }
              if(deletenode==1){
                Q.push(neighbour);
                vis[neighbour]=1;
              }
            }
          }          
        }else{
          for(ui i=L_pstart[cur]; i<L_pstart[cur]+L_degree[cur]; i++){
              ui neighbour=L_edges[i];
              deletenode=0;            
              if(!vis[neighbour]){
                  if(--colorful_r[neighbour][attr[cur]][color[cur]] <= 0){
                      colorful_d[neighbour][attr[cur]]--;
                  }
                  if(colorful_d[neighbour][attr[neighbour]]<alpha-1){
                    deletenode=1; 
                  }else{                
                    for(ui j=0; j<attr_size; j++){
                    if(j==attr[neighbour])
                      continue;
                    else if(colorful_d[neighbour][j]<alpha){
                      deletenode=1; 
                      break;
                      }
                    }
                  }
              if(deletenode==1){
                Q.push(neighbour);
                vis[neighbour]=1;
              }
            }
          } 
        }
    }   

    for(ui i=0; i<n; i++){
        if(vis[i])
          continue;
        colorful_d[i] = new ui[attr_size]();
        for(ui j=0; j<attr_size; j++)
            colorful_d[i][j]=0;
    }
    for(ui i=0; i<n; i++){
        if(vis[i])
          continue;
        delete colorful_r[i];
        colorful_r[i] = new ui*[attr_size];        
        for(ui j=0; j<attr_size; j++){
            colorful_r[i][j] = new ui[totol_max_color]();
            for(ui k=0; k<totol_max_color; k++)
                colorful_r[i][j][k]=0;
        }
    }

    for(ui i=0;i<n;i++){
      if(!vis[i]){
        for(ui j=pstart[i];j<pstart[i+1];j++){
            ui neighbour=edges[j];            
            if(!vis[neighbour]){

              if((colorful_r[i][attr[neighbour]][color[neighbour]]++)==0){
                  colorful_d[i][attr[neighbour]]++;
              }
            }
          }
        }
    }

    queue<ui> removed; 
    for(ui i=0;i<n;i++){
        if(!vis[i]){
          if(i>right_start){
            for(ui j=0;j<attr_size;j++){
              if(colorful_d[i][j]<alpha){
                removed.push(i);
                vis[i]=true;
                break;
              }
            }
          }else{
            for(ui j=0;j<attr_size;j++){
              if(colorful_d[i][j]<beta){
                removed.push(i);
                vis[i]=true;
                break;
              }
            }
          }
        }
    }

    while(!removed.empty()){
        ui u=removed.front();
        removed.pop();
        for(ui i=pstart[u];i<pstart[u+1];i++){
          ui j=edges[i];
          if(vis[j])
            continue;
          else{
            if(--colorful_r[j][attr[u]][color[u]]<=0){
                colorful_d[j][attr[u]]--; 
                if(j<right_start){
                  if(colorful_d[j][attr[u]]<beta){
                    vis[j]=true; 
                    removed.push(j);
                  }
                }else{
                  if(colorful_d[j][attr[u]]<alpha){
                    vis[j]=true; 
                    removed.push(j);
                  }                  
                }           

            }
          }    
    }
  }
  for(ui i=0;i<n;i++){
    if(!vis[i]){
      colorful_min[i]=max(colorful_d[i][0]+1,colorful_d[i][1]);
    }
  }
  printf("< time for colorful alpha beta pruned: %s (microseconds) >\n", Utility::integer_to_string(t.elapsed()).c_str());     
  printf("< here we end colorful alpha beta pruned >\n\n\n");     
  return ;
}

void BipartiteGraph::naive_enumeration_two_side_fairbiclique(){
  vector<ui> L,R,P,Q;
  baseline_initialize_lrpq(L,R,P,Q);
  vector<ui> R_size(2,0);
  vector<ui> P_size(2,0);
  vector<ui> L_size(2,0);
  for(ui i=0;i<P.size();i++){
    P_size[attr[P[i]]]++;
  }
  for(ui i=0;i<L.size();i++){
    L_size[attr[L[i]]]++;
  }
  naive_enumeration_two_side_fairbiclique(L,R,P,Q,L_size,R_size,P_size,0);
  return ;
}
vector<ui> BipartiteGraph::get_intersection_set(const vector<ui> plus,const vector<ui> ori){
  ui i=0,j=0;
  ui plus_size=plus.size(),ori_size=ori.size();
  vector<ui> res;
  while(i<plus_size&&j<ori_size){
    if(plus[i]==ori[j]){
      i++;
      j++;
      continue;
    }else if(ori[j]>plus[i]){
      res.emplace_back(plus[i]);
      i++;
    }else{
      printf("dont's exits\n");
      exit(1);
    }
  }
  vector<ui> attr_res(2,0);
  for(ui i=0;i<res.size();i++){
    attr_res[attr[res[i]]]++;
  }
  return attr_res;
}


bool BipartiteGraph::is_the_same_vec(const vector<ui> left,const vector<ui> right){
  ui left_size=left.size();
  ui right_size=right.size();
  vector<ui> input_left=left;
  vector<ui> input_right=right;
  sort(input_left.begin(),input_left.end());
  sort(input_right.begin(),input_right.end());
  if(left_size!=right_size){
    return false;
  }else{
    for(ui i=0;i<left_size;i++){
      if(input_left[i]==input_right[i]){
        continue;
      }else{
        return false;
      }
    }
    return true;
  }
  return true;
}
bool BipartiteGraph::is_a_maximal_fair_subset(const vector<ui> subset,const vector<ui> set){

  vector<ui> set_attr_size(2,0);
  vector<ui> subset_attr_size(2,0);
  ui i=0,j=0;
  while(j<set.size()){
    if(i==subset.size()){
      set_attr_size[attr[set[j]]]++;      
      j++;
    }
    else if(set[j]==subset[i]){
      subset_attr_size[attr[set[j]]]++;
      i++;
      j++;
    }else{
      set_attr_size[attr[set[j]]]++;
      j++;
    }
  }
  if(set_attr_size[0]>0&&set_attr_size[1]>0){
    return false;
  }else if(set_attr_size[0]>0){
    ui gap=get_ui_abs(subset_attr_size[0]+1,subset_attr_size[1]);
    if(gap<=delta){
      return false;
    }
  }else if(set_attr_size[1]>0){
    ui gap=get_ui_abs(subset_attr_size[0],subset_attr_size[1]+1);
    if(gap<=delta){
      return false;
    }    
  }else{
    return true;
  }
  return true;
}
void BipartiteGraph::naive_enumeration_two_side_fairbiclique(const vector<ui> &L_in, const vector<ui> &R_in, const vector<ui> &P_in, const vector<ui> &Q_in,
                                        const vector<ui> L_in_size, const vector<ui> R_in_size,const vector<ui> P_in_size,const ui Q_in_size) 
{
    //here we only consider 2 attribute situation
    vector<ui> L = L_in;
    vector<ui> R = R_in;
    vector<ui> P = P_in;
    vector<ui> Q = Q_in;
    vector<ui> L_size=L_in_size;
    vector<ui> R_size=R_in_size;
    vector<ui> P_size=P_in_size;
    ui Q_size=Q_in_size;

    ui total_P_size=P_size[0]+P_size[1];
    change_mem(L.size()*sizeof(ui),1);
    change_mem(R.size()*sizeof(ui),1);
    change_mem(P.size()*sizeof(ui),1);
    change_mem(Q.size()*sizeof(ui),1);
    while(total_P_size) {
        ui x=P[0];
        vector<ui> R_prime=R;
        vector<ui> L_prime;
        vector<ui> P_prime;
        vector<ui> Q_prime;  
        vector<ui> R_prime_size=R_size;
        vector<ui> L_prime_size(2,0);
        vector<ui> P_prime_size(2,0);      
        ui Q_prime_size=0;

        R_prime.emplace_back(x);
        R_prime_size[attr[x]]++;        

        bool is_maximal = true;      
        for (ui i = 0; i < L_size[0]+L_size[1]; i++) {
            ui u = L[i];
            for(ui j=pstart[u];j<pstart[u+1];j++){
                if(vis[edges[j]])
                    continue;
                if(edges[j]==x){
                    L_prime.emplace_back(u);
                    L_prime_size[attr[u]]++;
                    break;
                }
            }
        }
        ui total_L_prime_size=L_prime_size[0]+L_prime_size[1];
        vector<ui> Q_attr_size(2,0);
        for(ui i=0;i<Q_size;i++){
          ui v=Q[i];
          vector<ui> num_L_prime_neighbours=get_common_neighbour_attr(v,L_prime);
          if (num_L_prime_neighbours[0]+num_L_prime_neighbours[1]== total_L_prime_size) {
            Q_attr_size[attr[v]]++;
          } 
          if(num_L_prime_neighbours[0]>=alpha&&num_L_prime_neighbours[1]>=alpha){
            Q_prime.push_back(Q[i]);
            Q_prime_size++;
          }
        }

        if(L_prime_size[0]<alpha||L_prime_size[1]<alpha){
          is_maximal=false;          
        }

        if(Q_attr_size[0]>=1&&Q_attr_size[1]>=1){
          is_maximal=false;
        }

        if(is_maximal){
          vector<ui> R_prime_plus;
          vector<ui> R_prime_plus_size(2,0);
          for (ui i = 0; i < total_P_size; i++) {
            ui v = P[i];
            if(v == x) {
              continue;
            }
            vector<ui> num_L_prime_neighbours = get_common_neighbour_attr(v,L_prime);
            if (num_L_prime_neighbours[0]+num_L_prime_neighbours[1] == total_L_prime_size) {
              R_prime_plus.emplace_back(v);
              R_prime_plus_size[attr[v]]++;
              P_prime.emplace_back(v);
              P_prime_size[attr[v]]++;
            }
            else if (num_L_prime_neighbours[0] >= alpha&&num_L_prime_neighbours[1]>=alpha) {
              P_prime.emplace_back(v);
              P_prime_size[attr[v]]++;   
            }
          }
          if(R_prime_plus_size[0]+R_prime_plus_size[1]==P_prime_size[0]+P_prime_size[1]){
            ui gap=get_ui_abs(R_prime_size[0]+R_prime_plus_size[0],R_prime_size[1]+R_prime_plus_size[1]);
            if(gap<=delta){
              R_prime_size[0]=R_prime_size[0]+R_prime_plus_size[0];
              R_prime_size[1]=R_prime_size[1]+R_prime_plus_size[1];
              R_prime.insert(R_prime.end(),R_prime_plus.begin(),R_prime_plus.end());
              R_prime_plus_size[0]=0;
              R_prime_plus_size[1]=0;
              vector <ui>().swap(P_prime);
              P_prime_size[0]=0;
              P_prime_size[1]=0;
            }
          }

          ui gap=get_ui_abs(R_prime_size[0],R_prime_size[1]);

          if(gap>delta){
            is_maximal=false;
          }      
          if((R_prime_size[0]<beta)||(R_prime_size[1]<beta)){
            is_maximal=false;            
          }
          change_mem(L_prime.size()*sizeof(ui),1);
          change_mem(R_prime.size()*sizeof(ui),1);
          change_mem(P_prime.size()*sizeof(ui),1);
          change_mem(Q_prime.size()*sizeof(ui),1);
          if(is_maximal){    
                if(R_prime_plus_size[0]+Q_attr_size[0]>=1&&R_prime_plus_size[1]+Q_attr_size[1]>=1){
                  is_maximal=false;                
                }else if(R_prime_plus_size[0]+Q_attr_size[0]==0&&R_prime_plus_size[1]+Q_attr_size[1]>=1){
                  gap=get_ui_abs(R_prime_size[0],R_prime_size[1]+1);
                  if(gap<=delta){
                    is_maximal=false;
                  }
                }else if(R_prime_plus_size[1]+Q_attr_size[1]==0&&R_prime_plus_size[0]+Q_attr_size[0]>=1){
                  gap=get_ui_abs(R_prime_size[0]+1,R_prime_size[1]);                
                  if(gap<=delta){
                    is_maximal=false;
                  }                
                }                                 
                if(is_maximal){  
                  //here we got one side fair, however it's not double side fair/                     
                  gap=get_ui_abs(L_prime_size[0],L_prime_size[1]);
                  if(gap<=delta){
                    fairbiclique_append(L_prime,R_prime);
                  }else{
                   vector<vector<ui>> L_candidate=combination(L_prime);
                   sort(R_prime.begin(),R_prime.end());
                   
                   for(ui i=0;i<L_candidate.size();i++){
                     vector<ui> L_candidate_neighbour=get_common_neighbour_set(L_candidate[i]);
                      if(is_a_maximal_fair_subset(R_prime,L_candidate_neighbour)){
                          fairbiclique_append(L_candidate[i],R_prime);
                      } 
                    }
                  }
                }                
          }
          
        bool continue_maximal_flag=true;
        if((R_prime_size[0]+P_prime_size[0]<beta)||(R_prime_size[1]+P_prime_size[1]<beta)){
          continue_maximal_flag=false;
        }

        if ((P_prime_size[0]+P_prime_size[1])&&(continue_maximal_flag)) {
          naive_enumeration_two_side_fairbiclique(L_prime, R_prime, P_prime, Q_prime,L_prime_size,R_prime_size,P_prime_size,Q_prime_size);
        }
        change_mem(L_prime.size()*sizeof(ui),0);
        change_mem(R_prime.size()*sizeof(ui),0);
        change_mem(P_prime.size()*sizeof(ui),0);
        change_mem(Q_prime.size()*sizeof(ui),0);        
      }
      P_size[attr[x]]--;      
      P.erase(P.begin());
      total_P_size--;
      Q.emplace_back(x);
      Q_size++;    
    }
    change_mem(L.size()*sizeof(ui),0);
    change_mem(R.size()*sizeof(ui),0);
    change_mem(P.size()*sizeof(ui),0);
    change_mem(Q.size()*sizeof(ui),0);    
    return ;
}

void BipartiteGraph::two_side_construct_graph(){
    Timer t;
    printf("< here we start construct graph >\n");   
    R_degree= new ui[n]();
    R_pstart= new ui[n+1]();    
    ui **umap=new ui*[n];
    for(ui i=0;i<n;i++){
      umap[i]=new ui[attr_size]();
    }    
    ui idx=0;
    ui *common_neigh=new ui[n]();  
    
    for(ui i=right_start;i<n;i++){
      idx=0;
      if(vis[i])
        continue;
      for(ui j=pstart[i];j<pstart[i+1];j++){
        ui u=edges[j];
        if(vis[u])
          continue;
        else{
          for(ui k=pstart[u];k<pstart[u+1];k++){       
            ui v=edges[k];
            if(vis[v])
              continue;
            if(v!=i){
              umap[v][attr[u]]=umap[v][attr[u]]+1;
              if(umap[v][0]+umap[v][1]==1){
                common_neigh[idx++]=v;
              }
            }
          }            
        }
      }
      bool flag=true;
      for(ui j=0;j<idx;j++){
        ui k=common_neigh[j];
        if(umap[k][0]>=alpha&&umap[k][1]>=alpha){
          if(flag){
            flag=false;
            R_nodes.emplace_back(i);
            R_n++;
            R_pstart[i]=R_m;
            R_edges.emplace_back(k);
            R_m++;              
          }else{
            R_edges.emplace_back(k); 
            R_m++;
          }
        }
        umap[k][0]=0;
        umap[k][1]=0;
      }          
    } 
    R_nodes.emplace_back(n);
    R_pstart[n]=R_m;
    printf("R_m number is %d, R_n number is %d\n",R_m,R_n);  
    for(ui i=0;i<R_n;i++){
      R_degree[R_nodes[i]]=R_pstart[R_nodes[i+1]]-R_pstart[R_nodes[i]];
    }
    for(ui i=right_start;i<n;i++){
      if(R_degree[i]==0){
        vis[i]=1;
      }
    }    
    L_degree= new ui[n]();
    L_pstart= new ui[n+1]();    
    for(ui i=0;i<n;i++){
      common_neigh[i]=0;
    }
    idx=0;
    for(ui i=0;i<right_start;i++){
      idx=0;
      if(vis[i])
        continue;
      for(ui j=pstart[i];j<pstart[i+1];j++){
        ui u=edges[j];
        if(vis[u])
          continue;
        else{
          for(ui k=pstart[u];k<pstart[u+1];k++){       
            ui v=edges[k];
            if(vis[v])
              continue;
            if(v!=i){
              umap[v][attr[u]]=umap[v][attr[u]]+1;
              if(umap[v][0]+umap[v][1]==1){
                common_neigh[idx++]=v;
              }
            }
          }            
        }
      }
      bool flag=true;
      for(ui j=0;j<idx;j++){
        ui k=common_neigh[j];
        if(umap[k][0]>=beta&&umap[k][1]>=beta){
          if(flag){
            flag=false;
            L_nodes.emplace_back(i);
            L_n++;
            L_pstart[i]=L_m;
            L_edges.emplace_back(k);
            L_m++;              
          }else{
            L_edges.emplace_back(k); 
            L_m++;
          }
        }
        umap[k][0]=0;
        umap[k][1]=0;
      }          
    } 
    L_nodes.emplace_back(right_start);
    L_pstart[right_start]=L_m;
    printf("L_m number is %d, L_n number is %d\n",L_m,L_n);  
    for(ui i=0;i<L_n;i++){
      L_degree[L_nodes[i]]=L_pstart[L_nodes[i+1]]-L_pstart[L_nodes[i]];
    }        
    for(ui i=0;i<right_start;i++){
      if(L_degree[i]==0){
        vis[i]=1;
      }
    }     
    printf("< time for construct graph: %s (microseconds) >\n", Utility::integer_to_string(t.elapsed()).c_str());      
    printf("< here we end construct graph> \n\n\n");       
    return ;
}

void BipartiteGraph::naive_enumeration_two_side(){  
  two_side_alpha_beta_pruned();
  reformat_graph();
  two_side_colorful_alpha_beta_pruned();
  reformat_graph();
  Timer t;    
  naive_enumeration_two_side_fairbiclique();
  printf("mem of graph = %lld\n", max_mem);  
  printf("< time for find fairbiclique: %s (microseconds) >\n", Utility::integer_to_string(t.elapsed()).c_str());
  printf("< find maximal fair biclique: %d, find redundant maximal fair biclique: %d >\n", fairbiclique_count,redundant_count);     
  return ;
}
void BipartiteGraph::alpha_beta_pruned(queue<ui> removed){
    ui max_removed=removed.size();
    while(!removed.empty()){
      ui u=removed.front();
      removed.pop();
      for(ui i=pstart[u];i<pstart[u+1];i++){
        ui j=edges[i];
        if(vis[j])
          continue;
        if(j<right_start){
          if(--degree[j]<beta){                      
              vis[j]=true;
              removed.push(j);
              continue;
          }
        }else{
          if(--degree[j]<alpha){                      
              vis[j]=true;
              removed.push(j);
              continue;
          }          
        }
      } 
    }
    change_mem(sizeof(removed)+ sizeof(int)*max_removed, 0);
    return ;
}

void BipartiteGraph::alpha_beta_pruned(){
  Timer t;
  printf("< here we start alpha beta pruned >\n");  
  queue<ui> removed;    
  for(ui i=0;i<right_start;i++){
      if(degree[i]<beta){
        vis[i]=true;
        removed.push(i);
        break;
      }
  }  
  for(ui i=right_start;i<n;i++){
      if(degree[i]<alpha){
        vis[i]=true;
        removed.push(i);
        break;
      }
  }
  change_mem(sizeof(removed)+ sizeof(int)*removed.size(), 1);
  alpha_beta_pruned(removed);
  printf("<time for alpha beta pruned: %s (microseconds)>\n", Utility::integer_to_string(t.elapsed()).c_str());     
  printf("< here we end alpha beta pruned >\n\n\n");      
  return ;
}

void BipartiteGraph::two_side_alpha_beta_pruned(queue<ui> removed){
    ui max_removed=removed.size();
    while(!removed.empty()){
      ui u=removed.front();
      removed.pop();
      for(ui i=pstart[u];i<pstart[u+1];i++){
        ui j=edges[i];
        if(vis[j])
          continue;
        if(j<right_start){
          if(--attr_degree[j][attr[u]]<beta){                      
              vis[j]=true;
              removed.push(j);
              continue;
          }
        }else{
          if(--attr_degree[j][attr[u]]<alpha){                      
              vis[j]=true;
              removed.push(j);
              continue;
          }          
        }
      } 
    }
    change_mem(sizeof(removed)+ sizeof(int)*max_removed, 0);
    return ;
}

void BipartiteGraph::two_side_alpha_beta_pruned(){
  Timer t;
  printf("< here we start alpha beta pruned >\n");  
  queue<ui> removed;    
  for(ui i=0;i<right_start;i++){
    for(ui k=0;k<attr_size;k++){
      if(attr_degree[i][k]<beta){
        vis[i]=true;
        removed.push(i);
        break;
      }
    }
  }  
  for(ui i=right_start;i<n;i++){
    for(ui k=0;k<attr_size;k++){
      if(attr_degree[i][k]<alpha){
        vis[i]=true;
        removed.push(i);
        break;
      }
    }
  }
  change_mem(sizeof(removed)+ sizeof(int)*removed.size(), 1);
  two_side_alpha_beta_pruned(removed);
  printf("<time for alpha beta pruned: %s (microseconds)>\n", Utility::integer_to_string(t.elapsed()).c_str());     
  printf("< here we end alpha beta pruned >\n\n\n");      
  return ;
}
void BipartiteGraph::baseline_one_side_biclique_find(const vector<ui> &L_in, const vector<ui> &R_in, const vector<ui> &P_in, const vector<ui> &Q_in,
                                        const ui L_in_size, const ui R_in_size,const ui P_in_size,const ui Q_in_size) 
{
    vector<ui> L = L_in;
    vector<ui> R = R_in;
    vector<ui> P = P_in;
    vector<ui> Q = Q_in;
    ui L_size=L_in_size;
    ui R_size=R_in_size;
    ui P_size=P_in_size;
    ui Q_size=Q_in_size;
    
    while(P_size) {
        ui x=P.back();
        vector<ui> R_prime=R;
        vector<ui> L_prime;
        ui R_prime_size=R_size;
        ui P_prime_size=0;
        ui Q_prime_size=0;      
        R_prime.emplace_back(x);
        R_prime_size++;
        ui L_prime_size=0;
        for (ui i = 0; i < L_size; i++) {
            ui u = L[i];
            for(ui j=pstart[u];j<pstart[u+1];j++){
                if(vis[edges[j]])
                    continue;
                if(edges[j]==x){
                    L_prime.emplace_back(u);
                    L_prime_size++;
                    break;
                }
            }
        }      
        vector<ui> P_prime;
        vector<ui> Q_prime;
        bool is_maximal = true;
        
        for(ui j=0;j<Q_size;j++){
            ui v=Q[j];
            ui num_L_prime_neighbours=get_common_neighbour(v,L_prime);
            if (num_L_prime_neighbours == L_prime_size) {
                is_maximal = false;
                break;
            }
            else if (num_L_prime_neighbours >0) {
                Q_prime.emplace_back(v);
                Q_prime_size++;
            }                        
        }      
        if(is_maximal){
          for (ui j = 0; j < P_size; j++) {
              ui v = P[j];
              if(v == x) {
                continue;
              }
              ui num_L_prime_neighbours = get_common_neighbour(v,L_prime);
              if (num_L_prime_neighbours == L_prime_size) {
                R_prime.emplace_back(v);
                R_prime_size++;
              }
              else if (num_L_prime_neighbours >0) {
                P_prime.emplace_back(v);
                P_prime_size++;                  
              }
            }
//            printf("%d %d\n",L_prime_size,R_prime_size);
            if(L_prime_size>=alpha&&R_prime_size>=2*beta){
              fairbiclique_append(L_prime,R_prime);
            }
            if(!P_prime.empty()&&L_prime_size>=alpha){
                baseline_one_side_biclique_find(L_prime, R_prime, P_prime, Q_prime,L_prime_size,R_prime_size,P_prime_size,Q_prime_size);
            }          
        }
        P.pop_back();
        P_size--;
        Q.push_back(x);
        Q_size++;
    }
    return ;
}

void BipartiteGraph::baseline_two_side_biclique_find(const vector<ui> &L_in, const vector<ui> &R_in, const vector<ui> &P_in, const vector<ui> &Q_in,
                                        const ui L_in_size, const ui R_in_size,const ui P_in_size,const ui Q_in_size) 
{
    vector<ui> L = L_in;
    vector<ui> R = R_in;
    vector<ui> P = P_in;
    vector<ui> Q = Q_in;

    ui L_size=L_in_size;
    ui R_size=R_in_size;
    ui P_size=P_in_size;
    ui Q_size=Q_in_size;
    
    while(P_size) {
        ui x=P.back();
        vector<ui> R_prime=R;
        vector<ui> L_prime;
        ui R_prime_size=R_size;
        ui P_prime_size=0;
        ui Q_prime_size=0;      
        R_prime.emplace_back(x);
        R_prime_size++;
        ui L_prime_size=0;    

        for (ui i = 0; i < L_size; i++) {
            ui u = L[i];
            for(ui j=pstart[u];j<pstart[u+1];j++){
                if(vis[edges[j]])
                    continue;
                if(edges[j]==x){
                    L_prime.emplace_back(u);
                    L_prime_size++;
                    break;
                }
            }
        }      
        vector<ui> P_prime;
        vector<ui> Q_prime;
        bool is_maximal = true;
        
        for(ui j=0;j<Q_size;j++){
            ui v=Q[j];
            ui num_L_prime_neighbours=get_common_neighbour(v,L_prime);
            if (num_L_prime_neighbours == L_prime_size) {
                is_maximal = false;
                break;
            }
            else if (num_L_prime_neighbours >0) {
                Q_prime.emplace_back(v);
                Q_prime_size++;
            }                        
        }      
        if(is_maximal){
          for (ui j = 0; j < P_size; j++) {
              ui v = P[j];
              if(v == x) {
                continue;
              }
              ui num_L_prime_neighbours = get_common_neighbour(v,L_prime);
              if (num_L_prime_neighbours == L_prime_size) {
                R_prime.emplace_back(v);
                R_prime_size++;
              }
              else if (num_L_prime_neighbours >0) {
                P_prime.emplace_back(v);
                P_prime_size++;                  
              }
            }
            if(L_prime_size>=2*alpha&&R_prime_size>=2*beta){
              fairbiclique_append(L_prime,R_prime);
            }
            if(!P_prime.empty()&&L_prime_size>=2*alpha){
                baseline_one_side_biclique_find(L_prime, R_prime, P_prime, Q_prime,L_prime_size,R_prime_size,P_prime_size,Q_prime_size);
            }          
        }
        P.pop_back();
        P_size--;
        Q.push_back(x);
        Q_size++;
    }
    return ;
}

void BipartiteGraph::get_num_bicliques(){
  one_side_alpha_beta_pruned();
  reformat_graph();
  Timer t;  
  baseline_one_side_biclique_find();
  printf("mem of graph = %lld\n", max_mem);
  printf("< time for find biclique: %s (microseconds) >\n", Utility::integer_to_string(t.elapsed()).c_str());
  printf("< find maximal biclique: %d, find redundant maximal fair biclique: %d >\n", fairbiclique_count,redundant_count);     
  return ;
}