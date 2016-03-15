#ifndef ALL_POINTS_MAX_FLOW_HEADER_GUARD
#define ALL_POINTS_MAX_FLOW_HEADER_GUARD
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/undirected_dfs.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include "custom_dfs.hpp"
namespace allPointsMaxFlow
{
	template<class inputGraph, typename flowType = double> struct flowPredicate
	{
	public:
		typedef typename inputGraph::edge_descriptor edge_descriptor;
		flowPredicate()
			:graph(NULL), capacities(NULL), residualCapacities(NULL)
		{}
		flowPredicate(const inputGraph& graph, const std::vector<flowType>& capacities, const std::vector<flowType>& residualCapacities)
			: graph(&graph), capacities(&capacities), residualCapacities(&residualCapacities)
		{}
		bool operator()(const edge_descriptor& edge) const
		{
			int index = boost::get(boost::edge_index, *graph, edge);
			return (*residualCapacities)[index] > 0;// || (*residualCapacities)[reverseIndex] < (*capacities)[reverseIndex];
		}
	private:
		const inputGraph* graph;
		const std::vector<flowType>* capacities;
		const std::vector<flowType>* residualCapacities;
	};
	template<class inputGraph, typename flowType = double> struct allPointsMaxFlowScratch
	{
		typedef typename boost::property_map<inputGraph, boost::edge_index_t>::const_type edgeIndexMapType;
		typedef typename boost::property_map<inputGraph, boost::vertex_index_t>::const_type vertexIndexMapType;
		typedef boost::iterator_property_map<typename std::vector<flowType>::const_iterator, edgeIndexMapType> edgeCapacityMapType;
		typedef boost::iterator_property_map<typename std::vector<flowType>::iterator, edgeIndexMapType> edgeResidualCapacityMapType;
		typedef boost::iterator_property_map<typename std::vector<typename inputGraph::edge_descriptor>::iterator, vertexIndexMapType> vertexPredecessorMapType;
		typedef boost::iterator_property_map<typename std::vector<boost::default_color_type>::iterator, vertexIndexMapType> vertexColorMapType;
		typedef boost::iterator_property_map<typename std::vector<boost::default_color_type>::iterator, edgeIndexMapType> edgeColorMapType;
		typedef boost::iterator_property_map<typename std::vector<int>::iterator, vertexIndexMapType> distanceMapType;
		std::vector<flowType> edgeResidualCapacityVector;
		std::vector<typename inputGraph::edge_descriptor> vertexPredecessorVector;
		std::vector<boost::default_color_type> colorVector1;
		std::vector<boost::default_color_type> colorVector2;
		std::vector<int> distanceVector;

		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, boost::property<boost::edge_index_t, int, boost::property<boost::edge_weight_t, flowType> > > starGraphType;
		typedef typename boost::property_map<starGraphType, boost::edge_index_t>::const_type starEdgeIndexMapType;
		typedef typename boost::property_map<starGraphType, boost::vertex_index_t>::const_type starVertexIndexMapType;
		typedef typename boost::property_map<starGraphType, boost::edge_weight_t>::const_type starEdgeWeightMapType;
		typedef boost::iterator_property_map<typename std::vector<boost::default_color_type>::iterator, starVertexIndexMapType> starVertexColorMapType;
		typedef boost::iterator_property_map<typename std::vector<boost::default_color_type>::iterator, starEdgeIndexMapType> starEdgeColorMapType;
		
		typedef boost::filtered_graph<inputGraph, flowPredicate<inputGraph, flowType> > filteredGraph;
		typedef typename boost::property_map<filteredGraph, boost::vertex_index_t>::const_type filteredVertexIndexMapType;
		typedef typename boost::property_map<filteredGraph, boost::edge_index_t>::type filteredEdgeIndexMapType;
		typedef boost::iterator_property_map<typename std::vector<boost::default_color_type>::iterator, filteredVertexIndexMapType> filteredVertexColorMapType;
		typedef boost::iterator_property_map<typename std::vector<boost::default_color_type>::iterator, filteredEdgeIndexMapType> filteredEdgeColorMapType;
	};
	template<typename inputGraph, typename flowType = double> struct minimumEdgeWeightVisitorState
	{
	public:
		minimumEdgeWeightVisitorState(std::vector<flowType>& flowMatrix)
			: flowMatrix(flowMatrix)
		{}
		typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::vertex_descriptor start;
		std::vector<typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::edge_descriptor> edges;
		std::vector<flowType> minimumValues;
		std::vector<flowType>& flowMatrix;
		std::size_t nVertices;
	};
	template<class inputGraph, typename flowType> class minimumEdgeWeightVisitor : public boost::default_dfs_visitor
	{
	public:
		minimumEdgeWeightVisitor(minimumEdgeWeightVisitorState<inputGraph, flowType>& stateRef)
			:stateRef(stateRef)
		{}
		void tree_edge(const typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::edge_descriptor& e, const typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType& star)
		{
			flowType currentEdgeWeight = boost::get(boost::edge_weight, star, e);
			if(stateRef.minimumValues.size() == 0)
			{
				stateRef.minimumValues.push_back(currentEdgeWeight);
				stateRef.edges.push_back(e);
			}
			else if (currentEdgeWeight < *stateRef.minimumValues.rbegin())
			{
				stateRef.minimumValues.push_back(currentEdgeWeight);
				stateRef.edges.push_back(e);
			}
			if(e.m_target == stateRef.start) throw std::runtime_error("");
			stateRef.flowMatrix[e.m_target + stateRef.nVertices * stateRef.start] = stateRef.flowMatrix[stateRef.start + stateRef.nVertices * e.m_target] = *stateRef.minimumValues.rbegin();
		}
		void finish_edge(const typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::edge_descriptor e, const typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType& star)
		{
			if(stateRef.edges.size() == 0) return;
			typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::edge_descriptor top = *stateRef.edges.rbegin();
			if(e == top || (e.m_target == top.m_source && e.m_source == top.m_target))
			{
				stateRef.edges.pop_back();
				stateRef.minimumValues.pop_back();
			}
		}
	private:
		minimumEdgeWeightVisitorState<inputGraph, flowType>& stateRef;
	};
	template<class inputGraph, typename flowType> void extractFlow(const typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType& tree, std::vector<flowType>& flowMatrix, allPointsMaxFlowScratch<inputGraph, flowType>& scratch)
	{
		const std::size_t nVertices = boost::num_vertices(tree);
		const std::size_t nEdges = boost::num_edges(tree);
		scratch.colorVector1.resize(nVertices);
		scratch.colorVector2.resize(nEdges);
		if (nVertices == 2)
		{
			typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::edge_iterator begin, end;
			boost::tie(begin, end) = boost::edges(tree);
			flowMatrix[1] = flowMatrix[2] = boost::get(boost::edge_weight, tree, *begin);
		}
		else if (nVertices < 2)
		{
			return;
		}
		typename allPointsMaxFlowScratch<inputGraph, flowType>::starVertexIndexMapType vertexIndexMap = boost::get(boost::vertex_index, tree);
		typename allPointsMaxFlowScratch<inputGraph, flowType>::starEdgeIndexMapType edgeIndexMap = boost::get(boost::edge_index, tree);

		//Identify the minimum weight edge between any pair of vertices in currentVertices
		typename allPointsMaxFlowScratch<inputGraph, flowType>::starVertexColorMapType vertexColorMap(scratch.colorVector1.begin(), vertexIndexMap);
		typename allPointsMaxFlowScratch<inputGraph, flowType>::starEdgeColorMapType edgeColorMap(scratch.colorVector2.begin(), edgeIndexMap);
		minimumEdgeWeightVisitorState<typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType, flowType> visitorState(flowMatrix);
		visitorState.nVertices = nVertices;
		minimumEdgeWeightVisitor<typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType, flowType> visitor(visitorState);
		typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::vertex_iterator current, end;
		boost::tie(current, end) = boost::vertices(tree);
		for(; current != end; current++)
		{
			if(visitorState.edges.size() > 0 || visitorState.minimumValues.size() > 0) throw std::runtime_error("Internal error");
			visitorState.start = *current;
			std::fill(scratch.colorVector1.begin() + *current, scratch.colorVector1.end(), boost::white_color);
			std::fill(scratch.colorVector2.begin(), scratch.colorVector2.end(), boost::white_color);
			boost::custom_undirected_dfs(tree, visitor, vertexColorMap, edgeColorMap, *current);
		}
	}
	template<class inputGraph, typename flowType> void allPointsMaxFlow(std::vector<flowType>& flowMatrix, const std::vector<flowType>& capacities, const inputGraph& graph, allPointsMaxFlowScratch<inputGraph, flowType>& scratch)
	{
		const std::size_t nVertices = boost::num_vertices(graph);
		const std::size_t nEdges = boost::num_edges(graph);
		//scratch data for max flow algorithm
		scratch.edgeResidualCapacityVector.resize(nEdges);
		scratch.vertexPredecessorVector.resize(nEdges);
		scratch.colorVector1.resize(nVertices);
		scratch.colorVector2.resize(nEdges);
		scratch.distanceVector.resize(nVertices);

		typename allPointsMaxFlowScratch<inputGraph, flowType>::edgeIndexMapType edgeIndexMap = boost::get(boost::edge_index, graph);
		typename allPointsMaxFlowScratch<inputGraph, flowType>::vertexIndexMapType vertexIndexMap = boost::get(boost::vertex_index, graph);
		typename allPointsMaxFlowScratch<inputGraph, flowType>::edgeResidualCapacityMapType residualCapacityMap(scratch.edgeResidualCapacityVector.begin(), edgeIndexMap);
		typename allPointsMaxFlowScratch<inputGraph, flowType>::edgeCapacityMapType edgeCapacityMap(capacities.begin(), edgeIndexMap);
		typename allPointsMaxFlowScratch<inputGraph, flowType>::vertexPredecessorMapType vertexPredecessorMap(scratch.vertexPredecessorVector.begin(), vertexIndexMap);
		typename allPointsMaxFlowScratch<inputGraph, flowType>::vertexColorMapType colorMap(scratch.colorVector1.begin(), vertexIndexMap);
		typename allPointsMaxFlowScratch<inputGraph, flowType>::distanceMapType distanceMap(scratch.distanceVector.begin(), vertexIndexMap);

		//construct star graph
		typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType star(nVertices);
		for (std::size_t i = 1; i < nVertices; i++)
		{
			boost::add_edge(i, 0, (int)i-1, star);
		}
		for (std::size_t s = 1; s < nVertices; s++)
		{
			typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::edge_descriptor uniqueEdge = *boost::out_edges(s, star).first;
			typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::vertex_descriptor t = uniqueEdge.m_target;
			flowType maxFlow = boost::boykov_kolmogorov_max_flow(graph, s, t, boost::residual_capacity_map(residualCapacityMap).capacity_map(edgeCapacityMap).predecessor_map(vertexPredecessorMap).color_map(colorMap).distance_map(distanceMap));
			boost::put(boost::edge_weight, star, uniqueEdge, maxFlow);

			//get out the two components of graph which are seperated by the mincut we just found
			typename allPointsMaxFlowScratch<inputGraph, flowType>::filteredGraph filtered(graph, flowPredicate<inputGraph, flowType>(graph, capacities, scratch.edgeResidualCapacityVector));
			typename allPointsMaxFlowScratch<inputGraph, flowType>::filteredVertexIndexMapType filteredVertexMap = boost::get(boost::vertex_index, filtered);
			typename allPointsMaxFlowScratch<inputGraph, flowType>::filteredEdgeIndexMapType filteredEdgeMap = boost::get(boost::edge_index, filtered);
			typename allPointsMaxFlowScratch<inputGraph, flowType>::filteredVertexColorMapType filteredVertexColorMap(scratch.colorVector1.begin(), filteredVertexMap);
			typename allPointsMaxFlowScratch<inputGraph, flowType>::filteredEdgeColorMapType filteredEdgeColorMap(scratch.colorVector2.begin(), filteredEdgeMap);
			std::fill(scratch.colorVector1.begin(), scratch.colorVector1.end(), boost::white_color);
			std::fill(scratch.colorVector2.begin(), scratch.colorVector2.end(), boost::white_color);
			boost::dfs_visitor<boost::null_visitor> visitor = boost::make_dfs_visitor(boost::null_visitor());
			boost::detail::undir_dfv_impl(filtered, s, visitor, filteredVertexColorMap, filteredEdgeColorMap);
			//at this point the two components are white and black
			for (std::size_t i = s + 1; i < nVertices; i++)
			{
				bool isNeighbour = false;
				typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::out_edge_iterator current, end;
				boost::tie(current, end) = boost::out_edges(i, star);
				while (current != end)
				{
					if (current->m_target == t) 
					{
						isNeighbour = true;
						break;
					}
					current++;
				}
				if (scratch.colorVector1[i] == boost::black_color && isNeighbour)
				{
					int oldEdgeIndex = boost::get(boost::edge_index, star, *current);
					boost::remove_edge(*current, star);
					boost::add_edge(i, s, oldEdgeIndex, star);
				}
			}
		}
		std::fill(flowMatrix.begin(), flowMatrix.end(), std::numeric_limits<double>::infinity());
		extractFlow<inputGraph, flowType>(star, flowMatrix, scratch);
	}
}
#endif
