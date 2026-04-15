#include "graph64.hpp"

// static DEFAULTOPTIONS(options);
// statsblk(stats);
// setword nauty_workspace[160*MAXM];
set *nauty_gv;
static const uint64 ZRO_V_BIT = 0x00UL;
static const uint64 TWO_V_BIT = 0x01UL;
static const uint64 THR_V_BIT = 0x02UL;
static const uint64 FOU_V_BIT = 0x03UL;
static const uint64 V_BIT_MSK = 0x03UL;
static const uint64 ONE_E_BIT = 0x04UL;
static const uint64 TWO_E_BIT = 0x08UL;
static const uint64 THR_E_BIT = 0x0CUL;
static const uint64 E_BIT_MSK = 0x0CUL;

void init_graph(graph64 &g, short size, unsigned short num_vcolors, 
                unsigned short num_ecolors, bool directed) {
	for (int i = 0; i!=64 ; ++i) {
		g.matrix[i] = 0;
	}
	g.size = size;
	g.has_vertex_colors = (num_vcolors > 1);
	g.has_edge_colors = (num_ecolors > 1);
	g.num_vertex_colors = num_vcolors;
	g.num_edge_colors = num_ecolors;
	g.directed = directed;
	
	g.codestamp = ZRO_V_BIT;
	
	g.num_vertex_bits = 0;
	if (num_vcolors > 1)   {g.num_vertex_bits = 2; g.codestamp = TWO_V_BIT; }
	if (num_vcolors > 4)   {g.num_vertex_bits = 3; g.codestamp = THR_V_BIT; }	
	if (num_vcolors > 8)   {g.num_vertex_bits = 4; g.codestamp = FOU_V_BIT; }	
	   
	if (num_ecolors ==1)   {g.num_edge_bits = 1; g.codestamp |= ONE_E_BIT; }    
	if (num_ecolors > 1)   {g.num_edge_bits = 2; g.codestamp |= TWO_E_BIT; }
	if (num_ecolors > 3)   {g.num_edge_bits = 3; g.codestamp |= THR_E_BIT; }	
   
	
	// g.g_N = (g.has_edge_colors) ? size*size : size;  // Graph is uncolored at edges
	// g.g_M = (g.g_N + WORDSIZE - 1) / WORDSIZE;
	// for (int i = 0; i != g.g_N; ++i) 
	// 	EMPTYSET( ( GRAPHROW(g.nauty_g, i, g.g_M) ) , g.g_M);
	// options.writeautoms = FALSE;
	// options.getcanon = TRUE;

	// options.defaultptn = (g.has_edge_colors || g.has_vertex_colors) ? FALSE : TRUE;
	
	// if (directed) {
	// 	options.digraph = TRUE;
	// 	options.invarproc = adjacencies;
	// 	options.mininvarlevel = 1;
	// 	options.maxinvarlevel = 10;
	// }

	// nauty_check(WORDSIZE, g.g_M, g.g_N, NAUTYVERSIONID);

}

graphcode64 toHashCode(graph64 &g) {

	graphcode64 ret = 0ULL;
	
	if ((!g.has_vertex_colors) && (!g.has_edge_colors)) { //graph is not colored
	
		// Create nautypp graph from adjacency matrix
		nautypp::Graph npp_g(g.size);
		for (int i = 0; i != g.size; ++i) {
			for (int j = i+1; j != g.size; ++j) {
				if (get_element(g, i, j)) {
					npp_g.link(i, j);
				}
			}
		}
		
		// Read canonical adjacency matrix and build the hash
		for (int a = 0; a != g.size; ++a) {
			for (int b = 0; b != g.size; ++b) {
	            if (a!=b) {
				   ret <<= 1;
	   			   ret |= npp_g.has_edge(a, b) ? 1 : 0;
	            }
			}
		}

		ret <<= 4;
        ret |= g.codestamp;
		
	} else { // graph is colored

		// convert for nauty
		int index = 0;
		uint32 sortarray[MAXN];

		// Initialize sortarray with original vertices (for partition refinement)
		for (int i = 0; i != g.size; ++i) {
			sortarray[i] = (0UL | i) | (g.matrix[i*8 + i] << 16);
		}
		
		// Add intermediate vertices for colored edges
		register unsigned int edgecolor_ij;
		if (g.has_edge_colors) {
			for (int i = 0; i != g.size; ++i) {
				for (int j = 0; j != g.size; ++j) {
					if (i != j) {
						edgecolor_ij = get_element(g, i, j);
						// Create intermediate vertex for colored edges (color > 1)
						if (edgecolor_ij > 1) {
							sortarray[index] = (0UL | index) | (j<<8) | (i<<12) | (edgecolor_ij << 20);
							++index;
						}
					}
				}
			}
		}
		
		// Create the auxiliary graph with original + intermediate vertices
		nautypp::Graph npp_g(index);
		
		// Add uncolored edges (edge color = 1 or uncolored)
		for (int i = 0; i != g.size; ++i) {
			for (int j = i+1; j != g.size; ++j) {
				unsigned int ecolor = get_element(g, i, j);
				// Add if uncolored (0) or has color 1
				if (ecolor == 0 || ecolor == 1) {
					npp_g.link(i, j);
				}
			}
		}
		
		// Add edges through intermediate vertices for colored edges
		if (g.has_edge_colors) {
			int inter_count = 0;
			for (int i = 0; i != g.size; ++i) {
				for (int j = 0; j != g.size; ++j) {
					if (i != j) {
						edgecolor_ij = get_element(g, i, j);
						if (edgecolor_ij > 1) {
							int intermediate_vertex = g.size + inter_count;
							// Connect original vertices through intermediate vertex
							npp_g.link(i, intermediate_vertex);
							npp_g.link(intermediate_vertex, j);
							inter_count++;
						}
					}
				}
			}
		}
		
		// Extract hash from original vertices adjacency in canonical order
		for (int a = 0; a != g.size; ++a) {
			for (int b = 0; b != g.size; ++b) {
				ret <<= (a==b) ? g.num_vertex_bits : g.num_edge_bits;
				
				if (a == b) {
					// Vertex color from original matrix diagonal
					ret |= get_element(g, a, a);
				} else {
					// Edge color from original matrix
					ret |= get_element(g, a, b);
				}
			}
		}
		
		ret <<= 4;
		ret |= g.codestamp;
	}
	
	return ret;
}



//
// TODO
//
void readHashCode(graph64 &g, graphcode64 gc) {
     unsigned short vbits = 0, ebits = 1 , vmsk = 0, emsk = 0;
     switch (gc & V_BIT_MSK) {
            case ZRO_V_BIT : vbits = 0; vmsk =  0; break;
            case TWO_V_BIT : vbits = 2; vmsk =  3; break;
            case THR_V_BIT : vbits = 3; vmsk =  7; break;                        
            case FOU_V_BIT : vbits = 4; vmsk = 15; break;            
     }
     switch (gc & E_BIT_MSK) {
            case ONE_E_BIT : ebits = 1; emsk = 1; break;
            case TWO_E_BIT : ebits = 2; emsk = 3; break;
            case THR_E_BIT : ebits = 3; emsk = 7; break;                        
     }     
     
     gc >>= 4; // remove codestamp
     
     for (int a = g.size-1; a >= 0; --a) {
         for (int b = g.size-1; b >= 0; --b) {
             if (a==b) {
                set_element(g, a, b, gc & vmsk); 
                gc >>= vbits;         
             } else {
                set_element(g, a, b, gc & emsk);
                gc >>= ebits;                    
             }
         }
     }     
}


graphcode64 getGraphID(graph64 &g, graphcode64 gc) {
     unsigned short vbits = 0, ebits = 1 , vmsk = 0, emsk = 0;
     switch (gc & V_BIT_MSK) {
            case ZRO_V_BIT : vbits = 0; vmsk =  0; break;
            case TWO_V_BIT : vbits = 2; vmsk =  3; break;
            case THR_V_BIT : vbits = 3; vmsk =  7; break;                        
            case FOU_V_BIT : vbits = 4; vmsk = 15; break;            
     }
     switch (gc & E_BIT_MSK) {
            case ONE_E_BIT : ebits = 1; emsk = 1; break;
            case TWO_E_BIT : ebits = 2; emsk = 3; break;
            case THR_E_BIT : ebits = 3; emsk = 7; break;                        
     }     
     
     gc >>= 4; // remove codestamp
     
    // graph nau_c[MAXN * MAXM];
	// graph nau_g[MAXN * MAXM];
	// short gn = g.size;
	// short gm = (gn + WORDSIZE - 1) / WORDSIZE;
	// for (int i = 0; i != g.g_N; ++i) 
	// 	EMPTYSET( ( GRAPHROW(nau_g, i, gm) ) , gm);	

    // Reconstruct the adjacency matrix from the compressed hash code
	for (int a = g.size-1; a >= 0; --a) {
		for (int b = g.size-1; b >= 0; --b) {
			if (a == b) {
				gc >>= vbits;         
			} else {
				unsigned short ecolor = gc & emsk;
				if (ecolor > 0) {
					set_element(g, a, b, ecolor);
				}
				gc >>= ebits;
			}
		}
	}
	
	// Create nautypp graph from reconstructed adjacency matrix
	nautypp::Graph npp_g(g.size);
	
	// Populate edges from the reconstructed matrix
	for (int i = 0; i != g.size; ++i) {
		for (int j = i+1; j != g.size; ++j) {
			unsigned short edge_val = get_element(g, i, j);
			if (edge_val > 0) {
				npp_g.link(i, j);
			}
		}
	}
	
	// nautypp computes canonical form internally
	// Extract the canonical adjacency matrix and return its hash
	graphcode64 ret = 0;
	
	for (int a = 0; a != g.size; ++a) {
		for (int b = 0; b != g.size; ++b) {
			ret <<= 1;
			if (npp_g.has_edge(a, b)) {
				ret |= 1;
			}
		}
	}
	
	return ret;
}