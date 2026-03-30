#include <gtest/gtest.h>

#include "DirectedGraph.h"

class DirectedGraphTest : public ::testing::Test {};

using namespace pulmtln;

TEST_F(DirectedGraphTest, getVertices)
{
	DirectedGraph g;
	g.addVertex(1);
	g.addVertex(5);

	EXPECT_EQ(2, g.getVertices().size());
}

TEST_F(DirectedGraphTest, findCycles_oriented_graph)
{
	// 0->1->2
	//  <---
	DirectedGraph g; 
	g.addEdge(0, 1);
	g.addEdge(1, 2);
	g.addEdge(2, 0);

	ASSERT_EQ(1, g.findCycles().size());
	EXPECT_EQ(3, g.findCycles().front().size());
}

TEST_F(DirectedGraphTest, findCycles)
{
	{
		//5 ---->---- 6
		//^			  v
		//2 ---->---- 3
		//^			  v
		//1 ----<---- 0
		//^			  v	
		//4 ----<---- 7

		DirectedGraph g;
		g.addEdge(4, 1);
		g.addEdge(1, 2);
		g.addEdge(2, 5);
		g.addEdge(5, 6);
		g.addEdge(6, 3);
		g.addEdge(3, 0);
		g.addEdge(0, 7);
		g.addEdge(7, 4);
		g.addEdge(2, 3);
		g.addEdge(0, 1);
		EXPECT_EQ(4, g.findCycles().size());
	}
	{
		//5 ----<---- 6
		//v			  ^
		//2 ----<---- 3
		//v			  ^
		//1 ---->---- 0
		//v			  ^	
		//4 ---->---- 7

		DirectedGraph g;
		g.addEdge(1, 4);
		g.addEdge(2, 1);
		g.addEdge(5, 2);
		g.addEdge(6, 5);
		g.addEdge(3, 6);
		g.addEdge(0, 3);
		g.addEdge(7, 0);
		g.addEdge(4, 7);
		g.addEdge(3, 2);
		g.addEdge(1, 0);
		EXPECT_EQ(4, g.findCycles().size());
	}
}

TEST_F(DirectedGraphTest, split_1)
{
	DirectedGraph g;
	g.addEdge(10, 11);
	g.addVertex(12);

	EXPECT_EQ(2, g.split().size());
}

TEST_F(DirectedGraphTest, split_2)
{
	DirectedGraph g;
	g.addEdge(10, 11);
	g.addEdge(12, 13);

	EXPECT_EQ(2, g.split().size());
}

TEST_F(DirectedGraphTest, domainBoundaries)
{
	//
	//  1-----------2
	//  |///////////|
	//  |//4 --- 5//|
	//  |//|     |//|
	//  |//6 --- 7//|
	//  |///////////|
	//  8-----------9

	DirectedGraph g;
	g.addClosedPath({ 1,2,4 });
	g.addClosedPath({ 2,5,4 });
	g.addClosedPath({ 2,9,5 });
	g.addClosedPath({ 5,9,7 });
	g.addClosedPath({ 7,9,6 });
	g.addClosedPath({ 6,9,8 });
	g.addClosedPath({ 6,8,4 });
	g.addClosedPath({ 4,8,1 });

	auto allBoundaries{ g.getBoundaryGraph() };
	EXPECT_EQ(allBoundaries.edgesSize(), 8);
	EXPECT_EQ(allBoundaries.verticesSize(), 8);

	auto bdrs{ allBoundaries.split()};
	ASSERT_EQ(bdrs.size(), 2);
	ASSERT_EQ(bdrs[0].findCycles().size(), 1);
	EXPECT_EQ(Path({ 1, 2, 9, 8 }), bdrs[0].findCycles()[0]);
	EXPECT_EQ(bdrs[1].findCycles().size(), 1);

}