from typing import List, Tuple

class Graph:
    def __init__(self):
        self._nodes:List = []     
        self._edges:List[Tuple] = []

    @property
    def roots(self) -> List:
        roots:List = []
        for node in self._nodes:
            isChild = False
            for edge in self._edges:
                if edge[-1] == node:
                    isChild = True
                    continue
            if isChild == False:
                roots.append(node)
        return roots.copy()
    
    @property
    def nodes(self) -> List:
        return self._nodes.copy()

    @property
    def edges(self) -> List:
        return self._edges.copy()

    @nodes.setter
    def nodes(self, nodes):
        self._nodes = list(nodes)

    @edges.setter
    def edges(self, edges):
        self._edges = [tuple(e) for e in edges]

    def add_node(self, node):
        if node not in self._nodes:
            self._nodes.append(node)

    def add_edge(self, source, destination):
        if source not in self._nodes:
            self.add_node(source)
        if destination not in self._nodes:
            self.add_node(destination)
        if (source, destination) not in self._edges:
            self._edges.append((source, destination))

    def get_connections(self):
        connections = {node: [] for node in self._nodes}
        for source, destination in self._edges:
            connections[source].append(destination)
        return connections
    
    def getParentNodes(self) -> List:
        return [edge[0] for edge in self._edges]
    
    def getChildNodes(self) -> List:
        return [edge[-1] for edge in self._edges]
    
    #Necesita una revisión pero por ahora hace lo que necesito
    def prune_to_longest_paths(self):
        connections = self.get_connections()
        roots = [n for n in self._nodes if n not in self.getChildNodes()]

        longest_paths = []

        def dfs(node, path):
            path = path + [node]
            if node not in connections or not connections[node]:
                longest_paths.append(path)
                return
            for child in connections[node]:
                dfs(child, path)

        for root in roots:
            dfs(root, [])

        leaf_to_path = {}
        for path in longest_paths:
            leaf = path[-1]
            if leaf not in leaf_to_path or len(path) > len(leaf_to_path[leaf]):
                leaf_to_path[leaf] = path

        new_nodes = set()
        new_edges = set()
        for path in leaf_to_path.values():
            new_nodes.update(path)
            new_edges.update([(path[i], path[i+1]) for i in range(len(path)-1)])

        self._nodes = list(new_nodes)
        self._edges = list(new_edges)

    def __str__(self):
        return f"Graph(Nodes: {self._nodes},\n Edges: {self._edges})"