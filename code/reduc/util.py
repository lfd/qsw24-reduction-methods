import networkx as nx
import sys
import os

def MultiGraphToLatexTikz(G, radius):
        """
        Input: networkx MultiGraph with compact node names
                radius: Node radius (typ 2, big graphs 10)
        Returns: latex string of that drawn MultiGraph
        """
        
        out = "\\resizebox{\\linewidth}{!}{%\n\\begin{tikzpicture}\n\t\\begin{scope}[circle]\n\t\t\draw\n"
        
        
        size = G.number_of_nodes()
        i = 0
     
        for n in G.nodes():
            out += "\t\t(" + str((i*360.0)/size) +  ":" + str(radius) + ") \tnode (x" + str(n) + ") {$x_{" + str(n) + " }$}\n"
            i+=1
        
        out += "\t\t;\n\t\\end{scope}\n"
        
        
    
        nG = nx.Graph()
        nG.add_nodes_from(G.nodes())
        nG.add_edges_from(G.edges())
        
        out += "\t\\begin{scope}[-]\n"
        for e in nG.edges():
            num_ed = G.number_of_edges(e[0],e[1])
            if num_ed == 1:
                out += "\t\t\\draw (x" + str(e[0]) + ") to " + "(x" + str(e[1]) + ");\n"
            else :
                out += "\t\t\\draw (x" + str(e[0]) + ") to node[above, sloped] {$" + str(num_ed) + "$} (x" + str(e[1]) + ");\n"
            
        
        out += "\t\\end{scope}\n\end{tikzpicture}\n}"
        
        return out
