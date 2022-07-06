#include "../utilities/Defines.h"
#include "../utilities/Utility.h"
#include <unordered_map>

using namespace std;

int main(int argc, char *argv[]) {
	if(argc < 4) {
		printf("Usage: [1]exe [2]graph-dir [3]edge-list-filename [4]right_start [5]drop\n");
		printf("\tIn the edgelist file, lines starting with non-number characters are skipped!\n");
		return 0;
	}
  srand((unsigned int)time(NULL));
	vector<ui> nodes;
	vector<pair<ui,ui> > edges; 

	string dir = string(argv[1]);
	string file_name = string(argv[2]);
  ui right_start=atoi(argv[3]);
  ui droptype=atoi(argv[4]);  
	FILE *f = Utility::open_file((dir + string("/") + file_name).c_str(), "r");
	char buf[1024];
	ui a, b;
  float p=0.6;
	while(fgets(buf, 1024, f)) {
		char comment = 1;
		for(ui j = 0;buf[j] != '\0';j ++) if(buf[j] != ' '&&buf[j] != '\t') {
			if(buf[j] >= '0'&&buf[j] <= '9') comment = 0;
			break;
		}
		if(comment) continue;

		for(ui j = 0;buf[j] != '\0';j ++) if(buf[j] < '0'||buf[j] > '9') buf[j] = ' ';
		sscanf(buf, "%u%u", &a, &b);
 	  b=b+right_start;
    if(droptype==1){
   			  ui t = rand() % 100;
  			  if (t >= p*100){
           continue;
          }
    }   
		nodes.push_back(a);
		nodes.push_back(b);
		edges.push_back(make_pair(a,b));
		edges.push_back(make_pair(b,a));
	}

	fclose(f);

	sort(nodes.begin(), nodes.end());
	nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end());

	printf("min id = %u, max id = %u, n = %lu\n", nodes.front(), nodes.back(), nodes.size());
  
	sort(edges.begin(), edges.end());
	edges.erase(unique(edges.begin(), edges.end()), edges.end());

	map<ui,ui> M;
	for(ui i = 0;i < nodes.size();i ++) M[nodes[i]] = i;

	char preserved = 1;
	for(ui i = 0;i < nodes.size();i ++) if(nodes[i] != i) preserved = 0;
	if(!preserved) printf("Node ids are not preserved!\n");

  vector<pair<ui,ui>> new_edges;
  if(droptype==2){
      unordered_map<ui,ui> node_map;
      for(ui i=0;i<nodes.size();i++){
      	  ui t = rand() % 100;
  		  	if (t >= p*100){
             node_map[nodes[i]]=0;
             continue;
          }else{
            node_map[nodes[i]]=1;
          }
      }
      for(ui i=0;i<edges.size();i++){
        if(node_map[edges[i].first]==0||node_map[edges[i].second]==0){
          continue;
        }else{
          new_edges.emplace_back(edges[i]);
        }
      }
    edges=new_edges;  
  }
  
	ui n = nodes.size();
	ui m = edges.size();  
	printf("n = %s, (undirected) m = %s\n", Utility::integer_to_string(n).c_str(), Utility::integer_to_string(m/2).c_str());

	ui *pstart = new ui[n+1];
	ui *edge = new ui[m];
  ui *attr = new ui[n];
	ui j = 0;
	right_start=0;
	for(ui i = 0;i < n;i ++) {
    attr[i]=rand()%2;

		pstart[i] = j;
		while(j < m&&edges[j].first == nodes[i]) {
			edge[j] = M[edges[j].second];
			++ j;
		}
	}
 	pstart[n] = j;
  	for(ui i=0;i < n;i ++){
      if(edge[pstart[i]]>i){
        right_start++;
      }
  	}
	printf("%d\n",right_start);

	f = Utility::open_file((dir + string("/b_degree.bin")).c_str(), "wb");

	ui tt = sizeof(ui);
	fwrite(&tt, sizeof(ui), 1, f);
	fwrite(&n, sizeof(ui), 1, f);
	fwrite(&m, sizeof(ui), 1, f);
	fwrite(&right_start, sizeof(ui), 1, f);
	ui *degree = new ui[n];
	for(ui i = 0;i < n;i ++) degree[i] = pstart[i+1]-pstart[i];
	fwrite(degree, sizeof(ui), n, f);
	fclose(f);

	f = Utility::open_file((dir + string("/b_adj.bin")).c_str(), "wb");
	fwrite(edge, sizeof(ui), m, f);
	fclose(f);

	f = Utility::open_file((dir + string("/attr.bin")).c_str(), "wb");
	fwrite(attr, sizeof(ui), n, f);
	fclose(f);
 
	delete[] pstart;
	delete[] edge;
	delete[] degree;
}