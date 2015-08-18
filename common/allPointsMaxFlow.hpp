#ifndef ALL_POINTS_MAX_FLOW_HEADER_GUARD
#define ALL_POINTS_MAX_FLOW_HEADER_GUARD
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/iterator/counting_iterator.hpp>
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
		typedef boost::iterator_property_map<typename std::vector<boost::default_color_type>::iterator, vertexIndexMapType> colorMapType;
		typedef boost::iterator_property_map<typename std::vector<int>::iterator, vertexIndexMapType> distanceMapType;
		std::vector<flowType> edgeResidualCapacityVector;
		std::vector<typename inputGraph::edge_descriptor> vertexPredecessorVector;
		std::vector<boost::default_color_type> colorVector;
		std::vector<int> distanceVector;

		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, boost::property<boost::edge_index_t, int, boost::property<boost::edge_weight_t, flowType> > > starGraphType;
		typedef typename boost::property_map<starGraphType, boost::edge_index_t>::type starEdgeIndexMapType;
		typedef typename boost::property_map<starGraphType, boost::edge_weight_t>::type starEdgeWeightMapType;
		std::vector<std::vector<int> > scratchIntVectors;
		
		std::vector<bool> minimumEdgeSearchVertices;

		typedef boost::filtered_graph<inputGraph, flowPredicate<inputGraph, flowType> > filteredGraph;
		typedef typename boost::property_map<filteredGraph, boost::vertex_index_t>::const_type filteredVertexIndexMapType;
		typedef boost::iterator_property_map<typename std::vector<boost::default_color_type>::iterator, filteredVertexIndexMapType> filteredColorMapType;
	};
	template<typename flowType = double> struct minimumEdgeWeightVisitorState
	{
	public:
		minimumEdgeWeightVisitorState(std::vector<bool>& minimumEdgeSearchVertices)
			:minimumEdgeSearchVertices(minimumEdgeSearchVertices), minimumEdgeIndex(-1), minimumWeight(std::numeric_limits<flowType>::max()), sourceVertex(-1), targetVertex(-1)
		{}
		std::vector<bool>& minimumEdgeSearchVertices;
		int minimumEdgeIndex;
		flowType minimumWeight;
		int sourceVertex;
		int targetVertex;
	};
	template<class inputGraph, typename flowType> class minimumEdgeWeightVisitor : public boost::default_dfs_visitor
	{
	public:
		minimumEdgeWeightVisitor(minimumEdgeWeightVisitorState<flowType>& stateRef)
			:stateRef(stateRef)
		{}
		void tree_edge(const typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::edge_descriptor& e, const typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType& star)
		{
			flowType currentEdgeWeight = boost::get(boost::edge_weight, star, e);
			if (currentEdgeWeight < stateRef.minimumWeight && stateRef.minimumEdgeSearchVertices[e.m_source] && stateRef.minimumEdgeSearchVertices[e.m_target])
			{
				stateRef.minimumWeight = currentEdgeWeight;
				stateRef.minimumEdgeIndex = boost::get(boost::edge_index, star, e);
				stateRef.sourceVertex = (int)e.m_source;
				stateRef.targetVertex = (int)e.m_target;
			}
		}
	private:
		minimumEdgeWeightVisitorState<flowType>& stateRef;
	};
	template<class inputGraph, typename flowType> void extractFlow(const typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType& tree, std::vector<flowType>& flowMatrix, std::vector<int>& currentVertices, allPointsMaxFlowScratch<inputGraph, flowType>& scratch)
	{
		const std::size_t nVertices = scratch.colorVector.size();
		if (currentVertices.size() == 2)
		{
			std::size_t vertexOne = currentVertices[0], vertexTwo = currentVertices[1];
			typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::out_edge_iterator current, end;
			boost::tie(current, end) = boost::out_edges(vertexOne, tree);
			while (current->m_target != vertexTwo)
			{
				current++;
			}
			flowMatrix[vertexOne + nVertices * vertexTwo] = flowMatrix[vertexTwo + nVertices * vertexOne] = boost::get(boost::edge_weight, tree, *current);
		}
		else if (currentVertices.size() == 1)
		{
			return;
		}
		else if (currentVertices.size() == 0)
		{
			throw std::runtime_error("This branch should be impossible");
		}
		typename allPointsMaxFlowScratch<inputGraph, flowType>::vertexIndexMapType vertexIndexMap = boost::get(boost::vertex_index, tree);

		//Identify the minimum weight edge between any pair of vertices in currentVertices
		std::fill(scratch.colorVector.begin(), scratch.colorVector.end(), boost::black_color);
		scratch.minimumEdgeSearchVertices.resize(nVertices);
		std::fill(scratch.minimumEdgeSearchVertices.begin(), scratch.minimumEdgeSearchVertices.end(), false);
		for (std::size_t i = 0; i < currentVertices.size(); i++)
		{
			scratch.colorVector[currentVertices[i]] = boost::white_color;
			scratch.minimumEdgeSearchVertices[currentVertices[i]] = true;
		}
		typename allPointsMaxFlowScratch<inputGraph, flowType>::colorMapType colorMap(scratch.colorVector.begin(), vertexIndexMap);
		minimumEdgeWeightVisitorState<flowType> visitorState(scratch.minimumEdgeSearchVertices);
		minimumEdgeWeightVisitor<inputGraph, flowType> visitor(visitorState);
		boost::depth_first_search(tree, visitor, colorMap, currentVertices[0]);

		//get out the two vertices of that edge
		int targetVertex = visitorState.targetVertex;
		int sourceVertex = visitorState.sourceVertex;
		flowType minimumWeight = visitorState.minimumWeight;
		//Identify all vertices on one side or the other.
		std::fill(scratch.colorVector.begin(), scratch.colorVector.end(), boost::white_color);
		scratch.colorVector[sourceVertex] = scratch.colorVector[targetVertex] = boost::black_color;
		boost::detail::depth_first_visit_impl(tree, sourceVertex, visitor, colorMap, boost::detail::nontruth2());
		scratch.colorVector[targetVertex] = boost::white_color;

		//get out vertices on both sides of the bottleneck
		std::vector<int> sideOne, sideTwo;
		if (scratch.scratchIntVectors.size() > 0)
		{
			sideOne = std::move(*scratch.scratchIntVectors.rbegin());
			scratch.scratchIntVectors.pop_back();
		}
		if (scratch.scratchIntVectors.size() > 0)
		{
			sideTwo = std::move(*scratch.scratchIntVectors.rbegin());
			scratch.scratchIntVectors.pop_back();
		}
		sideOne.clear();
		sideTwo.clear();

		for (std::size_t i = 0; i < currentVertices.size(); i++)
		{
			if (scratch.colorVector[currentVertices[i]] == boost::white_color)
			{
				sideOne.push_back(currentVertices[i]);
			}
			else sideTwo.push_back(currentVertices[i]);
			for (std::size_t j = i + 1; j < currentVertices.size(); j++)
			{
				if (scratch.colorVector[currentVertices[i]] != scratch.colorVector[currentVertices[j]])
				{
					flowMatrix[currentVertices[i] + currentVertices[j] * nVertices] = flowMatrix[currentVertices[j] + currentVertices[i]*nVertices] = minimumWeight;
				}
			}
		}
		extractFlow(tree, flowMatrix, sideOne, scratch);
		extractFlow(tree, flowMatrix, sideTwo, scratch);
		scratch.scratchIntVectors.push_back(std::move(sideOne));
		scratch.scratchIntVectors.push_back(std::move(sideTwo));
	}
	template<class inputGraph, typename flowType> void allPointsMaxFlow(std::vector<flowType>& flowMatrix, const std::vector<flowType>& capacities, const inputGraph& graph, allPointsMaxFlowScratch<inputGraph, flowType>& scratch)
	{
		const std::size_t nVertices = boost::num_vertices(graph);
		const std::size_t nEdges = boost::num_edges(graph);
		//scratch data for max flow algorithm
		scratch.edgeResidualCapacityVector.resize(nEdges);
		scratch.vertexPredecessorVector.resize(nEdges);
		scratch.colorVector.resize(nVertices);
		scratch.distanceVector.resize(nVertices);

		typename allPointsMaxFlowScratch<inputGraph, flowType>::edgeIndexMapType edgeIndexMap = boost::get(boost::edge_index, graph);
		typename allPointsMaxFlowScratch<inputGraph, flowType>::vertexIndexMapType vertexIndexMap = boost::get(boost::vertex_index, graph);
		typename allPointsMaxFlowScratch<inputGraph, flowType>::edgeResidualCapacityMapType residualCapacityMap(scratch.edgeResidualCapacityVector.begin(), edgeIndexMap);
		typename allPointsMaxFlowScratch<inputGraph, flowType>::edgeCapacityMapType edgeCapacityMap(capacities.begin(), edgeIndexMap);
		typename allPointsMaxFlowScratch<inputGraph, flowType>::vertexPredecessorMapType vertexPredecessorMap(scratch.vertexPredecessorVector.begin(), vertexIndexMap);
		typename allPointsMaxFlowScratch<inputGraph, flowType>::colorMapType colorMap(scratch.colorVector.begin(), vertexIndexMap);
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
			typename allPointsMaxFlowScratch<inputGraph, flowType>::filteredColorMapType filteredColorMap(scratch.colorVector.begin(), filteredVertexMap);
			std::fill(scratch.colorVector.begin(), scratch.colorVector.end(), boost::white_color);
			boost::dfs_visitor<boost::null_visitor> visitor = boost::make_dfs_visitor(boost::null_visitor());
			boost::detail::depth_first_visit_impl(filtered, s, visitor, filteredColorMap, boost::detail::nontruth2());
			//at this point the two components are white and black
			for (std::size_t i = s + 1; i < nVertices; i++)
			{
				bool isNeighbour = false;
				typename allPointsMaxFlowScratch<inputGraph, flowType>::starGraphType::out_edge_iterator current, end;
				boost::tie(current, end) = boost::out_edges(i, star);
				while (current != end)
				{
					if (current->m_target == t) isNeighbour = true;
					current++;
				}
				if (scratch.colorVector[i] == boost::black_color && isNeighbour)
				{
					boost::remove_edge(i, t, star);
					boost::add_edge(i, s, star);
				}
			}
		}
		std::vector<int> allVertices;
		if (scratch.scratchIntVectors.size() > 0)
		{
			allVertices = std::move(*scratch.scratchIntVectors.rbegin());
			scratch.scratchIntVectors.pop_back();
		}
		allVertices.clear();
		allVertices.insert(allVertices.begin(), boost::counting_iterator<int>(0), boost::counting_iterator<int>((int)nVertices));

		extractFlow<inputGraph, flowType>(star, flowMatrix, allVertices, scratch);
		
		scratch.scratchIntVectors.push_back(std::move(allVertices));
	}
}
#endif
