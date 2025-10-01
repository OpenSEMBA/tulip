import unittest

from src.Graph import Graph

class TestGraph(unittest.TestCase):

    def setUp(self) -> None:
        self.graph = Graph()

    def tearDown(self):
        del self.graph

    def test_addNode(self) -> None:
        self.graph.add_node('A')
        self.assertIn('A', self.graph.nodes)

        self.graph.add_node('A')
        self.assertEqual(self.graph.nodes.count('A'), 1)

    def test_addEdge(self) -> None:
        self.graph.add_edge('A', 'B')
        self.assertIn(('A', 'B'), self.graph.edges)

        self.assertIn('A', self.graph.nodes)
        self.assertIn('B', self.graph.nodes)

        self.graph.add_edge('A', 'B')
        self.assertEqual(self.graph.edges.count(('A', 'B')), 1)

    def test_settersAndGetters(self) -> None:
        nodes = ['X', 'Y']
        edges = [('X', 'Y')]
        self.graph.nodes = nodes
        self.graph.edges = edges
        self.assertEqual(self.graph.nodes, nodes)
        self.assertEqual(self.graph.edges, edges)

    def test_GetConnections(self) -> None:
        self.graph.add_edge('A', 'B')
        self.graph.add_edge('A', 'C')
        self.graph.add_node('D')  # no connections
        connections = self.graph.get_connections()
        expected = {
            'A': ['B', 'C'],
            'B': [],
            'C': [],
            'D': []
        }
        self.assertEqual(connections, expected)

    def test_str(self) -> None:
        self.graph.add_edge('A', 'B')
        s = str(self.graph)
        self.assertIn('A', s)
        self.assertIn('B', s)
        self.assertIn('Edges', s)

    def testPruneToLongestPaths(self) -> None:
        self.graph.nodes = ['A' ,'B', 'C', 'D', 'E', 'F', 'G']
        self.graph.edges = [
            ('A', 'B'), ('A', 'C'), ('A', 'D'), ('A', 'E'),
            ('B', 'C'), ('B', 'E'),
            ('C', 'E'),
            ('F', 'G')
        ]

        expectedEdges = [
            ('A', 'B'), ('A', 'D'),
            ('B', 'C'),
            ('C', 'E'),
            ('F', 'G')
        ]

        self.graph.prune_to_longest_paths()
        self.assertListEqual(sorted(self.graph.edges), sorted(expectedEdges))

    def testGetRoots(self) -> None:
        self.graph.nodes = ['A' ,'B', 'C', 'D', 'E', 'F', 'G']
        self.graph.edges = [
            ('A', 'B'), ('A', 'D'),
            ('B', 'C'),
            ('C', 'E'),
            ('F', 'G')
        ]
        self.assertListEqual(self.graph.roots, ['A', 'F'])

if __name__ == '__main__':
    unittest.main()